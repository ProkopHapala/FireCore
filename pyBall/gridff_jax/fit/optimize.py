"""Parameter fitting for the hybrid student model with JAX/Optax or SciPy fallback."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from ..array_api import HAS_JAX, HAS_OPTAX, as_numpy
from ..config import TrainingConfig
from ..hybrid_energy.model import HybridGridFFModel, HybridParameters, default_hybrid_parameters
from ..interfaces import DensityData, PoseBatch, TeacherResult

if HAS_JAX:
    import jax
    import jax.numpy as jnp

if HAS_OPTAX:
    import optax


@dataclass
class FitResult:
    params: HybridParameters
    history: list[dict] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    optimizer_result: dict | None = None


class ParameterLayout:
    def __init__(self, elements, base: HybridParameters, training: TrainingConfig, z_ref: float, use_req_plq: bool = False):
        self.elements = tuple(elements)
        self.training = training
        self.z_ref = float(z_ref)
        self.base = base
        self.use_req_plq = bool(use_req_plq)
        self.names = []
        self.bounds = []
        self.physical_bounds = []
        self.group_slices = {}
        self.scalar_indices = {}
        self.transforms = {}
        self._use_qeq = bool(training.fit_chi or training.fit_hardness or training.fit_image_plane)

        if self.use_req_plq:
            if training.fit_req_radius_offset:
                self._append_group("req_radius_offset", (-0.50, 1.20))
            if training.fit_req_energy_scale:
                self._append_group("req_energy_scale", (np.log(0.1), np.log(30.0)), physical_bounds=(0.1, 30.0), transform="exp")
        else:
            self._append_group("pauli", (0.05, 5.0))
            self._append_group("london", (0.05, 5.0))
        if training.fit_c6_coeff:
            self._append_group("c6_coeff", (-500.0, 500.0))
        # DISABLED(2026-03): density-derived channels confirmed ineffective in ablation scan.
        # if training.fit_density_lap:
        #     self._append_group("density_lap", (-5.0, 5.0))
        # if training.fit_density_grad:
        #     self._append_group("density_grad", (-5.0, 5.0))
        # if training.fit_london_alt:
        #     self._append_group("london_alt", (-5.0, 5.0))
        if training.fit_static_charge:
            self._append_group("static_charge", (-0.95, 2.0))
        if training.fit_reactive:
            self._append_group("reactive", (0.0, 5.0))
        if training.fit_chi:
            self._append_group("chi", (-15.0, 1.0))
        if training.fit_hardness:
            self._append_group("hardness", (0.5, 25.0))
        if training.fit_sample_shift_z:
            self._append_scalar("sample_shift_z", (-1.5, 1.5))
        if training.fit_coulomb_shift_z:
            self._append_scalar("coulomb_shift_z", (-1.5, 1.5))
        if self._use_qeq:
            self._append_scalar("image_scale", (0.0, 3.0))
            self._append_scalar("reservoir_chi", (-15.0, 1.0))
            self._append_scalar("reservoir_hardness", (0.5, 25.0))
        if training.fit_image_plane:
            self._append_scalar("image_plane", (z_ref - 0.5, z_ref + 3.0))

        self.lower = np.asarray([item[0] for item in self.bounds], dtype=float)
        self.upper = np.asarray([item[1] for item in self.bounds], dtype=float)

    def _append_group(self, prefix: str, bounds, physical_bounds=None, transform: str = "identity"):
        start = len(self.names)
        for element in self.elements:
            self.names.append((prefix, element))
            self.bounds.append(bounds)
            self.physical_bounds.append(bounds if physical_bounds is None else physical_bounds)
        self.group_slices[prefix] = slice(start, len(self.names))
        self.transforms[prefix] = transform

    def _append_scalar(self, prefix: str, bounds):
        self.scalar_indices[prefix] = len(self.names)
        self.names.append((prefix, None))
        self.bounds.append(bounds)
        self.physical_bounds.append(bounds)
        self.transforms[prefix] = "identity"

    def _to_physical(self, prefix: str, value, xp):
        transform = self.transforms.get(prefix, "identity")
        if transform == "exp":
            return xp.exp(value)
        return value

    def _to_internal(self, prefix: str, value):
        transform = self.transforms.get(prefix, "identity")
        if transform == "exp":
            return float(np.log(max(float(value), 1.0e-12)))
        return float(value)

    def _base_tree(self, use_jax: bool):
        xp = jnp if use_jax and HAS_JAX else np
        return {
            "pauli": xp.asarray([self.base.pauli[element] for element in self.elements], dtype=xp.float64),
            "london": xp.asarray([self.base.london[element] for element in self.elements], dtype=xp.float64),
            "reactive": xp.asarray([self.base.reactive[element] for element in self.elements], dtype=xp.float64),
            "static_charge": xp.asarray([self.base.static_charge[element] for element in self.elements], dtype=xp.float64),
            "c6_coeff": xp.asarray([self.base.c6_coeff.get(element, 1.0) for element in self.elements], dtype=xp.float64),
            "density_lap": xp.asarray([self.base.density_lap.get(element, 0.0) for element in self.elements], dtype=xp.float64),
            "density_grad": xp.asarray([self.base.density_grad.get(element, 0.0) for element in self.elements], dtype=xp.float64),
            "london_alt": xp.asarray([self.base.london_alt.get(element, 0.0) for element in self.elements], dtype=xp.float64),
            "req_radius_offset": xp.asarray([self.base.req_radius_offset[element] for element in self.elements], dtype=xp.float64),
            "req_energy_scale": xp.asarray([self.base.req_energy_scale[element] for element in self.elements], dtype=xp.float64),
            "chi": xp.asarray([self.base.chi[element] for element in self.elements], dtype=xp.float64),
            "hardness": xp.asarray([self.base.hardness[element] for element in self.elements], dtype=xp.float64),
            "image_scale": xp.asarray(self.base.image_scale, dtype=xp.float64),
            "image_plane": xp.asarray(self.base.image_plane, dtype=xp.float64),
            "image_damping": xp.asarray(self.base.image_damping, dtype=xp.float64),
            "sample_shift_z": xp.asarray(self.base.sample_shift_z, dtype=xp.float64),
            "coulomb_shift_z": xp.asarray(self.base.coulomb_shift_z, dtype=xp.float64),
            "reservoir_chi": xp.asarray(self.base.reservoir_chi, dtype=xp.float64),
            "reservoir_hardness": xp.asarray(self.base.reservoir_hardness, dtype=xp.float64),
            "total_charge": xp.asarray(self.base.total_charge, dtype=xp.float64),
            "energy_offset": xp.asarray(self.base.energy_offset, dtype=xp.float64),
        }

    def vector_to_tree(self, vector, use_jax: bool):
        xp = jnp if use_jax and HAS_JAX else np
        vector = xp.asarray(vector, dtype=xp.float64)
        tree = self._base_tree(use_jax=use_jax)
        for prefix in ("pauli", "london", "c6_coeff", "density_lap", "density_grad", "london_alt", "static_charge", "reactive", "chi", "hardness", "req_radius_offset", "req_energy_scale"):
            if prefix in self.group_slices:
                tree[prefix] = self._to_physical(prefix, vector[self.group_slices[prefix]], xp)
        for prefix in ("sample_shift_z", "coulomb_shift_z", "image_scale", "reservoir_chi", "reservoir_hardness", "image_plane"):
            if prefix in self.scalar_indices:
                tree[prefix] = self._to_physical(prefix, vector[self.scalar_indices[prefix]], xp)
        return tree

    def pack(self, params: HybridParameters):
        values = []
        for prefix, element in self.names:
            if prefix in ("pauli", "london", "c6_coeff", "density_lap", "density_grad", "london_alt", "static_charge", "reactive", "chi", "hardness", "req_radius_offset", "req_energy_scale"):
                values.append(self._to_internal(prefix, getattr(params, prefix)[element]))
            else:
                values.append(self._to_internal(prefix, getattr(params, prefix)))
        return np.asarray(values, dtype=float)

    def unpack(self, vector):
        vector = np.asarray(vector, dtype=float)
        params = HybridParameters(
            pauli=dict(self.base.pauli),
            london=dict(self.base.london),
            reactive=dict(self.base.reactive),
            static_charge=dict(self.base.static_charge),
            c6_coeff=dict(self.base.c6_coeff) if self.base.c6_coeff else {},
            density_lap=dict(self.base.density_lap) if self.base.density_lap else {},
            density_grad=dict(self.base.density_grad) if self.base.density_grad else {},
            london_alt=dict(self.base.london_alt) if self.base.london_alt else {},
            req_radius_offset=dict(self.base.req_radius_offset),
            req_energy_scale=dict(self.base.req_energy_scale),
            chi=dict(self.base.chi),
            hardness=dict(self.base.hardness),
            image_scale=self.base.image_scale,
            image_plane=self.base.image_plane,
            image_damping=self.base.image_damping,
            sample_shift_z=self.base.sample_shift_z,
            coulomb_shift_z=self.base.coulomb_shift_z,
            reservoir_chi=self.base.reservoir_chi,
            reservoir_hardness=self.base.reservoir_hardness,
            total_charge=self.base.total_charge,
            energy_offset=self.base.energy_offset,
        )
        for value, (prefix, element) in zip(vector, self.names):
            physical_value = self._to_physical(prefix, float(value), np)
            if prefix in ("pauli", "london", "c6_coeff", "density_lap", "density_grad", "london_alt", "static_charge", "reactive", "chi", "hardness", "req_radius_offset", "req_energy_scale"):
                getattr(params, prefix)[element] = float(physical_value)
            else:
                setattr(params, prefix, float(physical_value))
        return params

    def regularization(self, vector, use_jax: bool):
        xp = jnp if use_jax and HAS_JAX else np
        vector = xp.asarray(vector, dtype=xp.float64)
        reg = xp.asarray(0.0, dtype=xp.float64)
        if "req_radius_offset" in self.group_slices and self.training.req_radius_regularization > 0.0:
            radius_offset = vector[self.group_slices["req_radius_offset"]]
            reg = reg + xp.asarray(self.training.req_radius_regularization, dtype=xp.float64) * xp.mean(radius_offset * radius_offset)
        if "req_energy_scale" in self.group_slices and self.training.req_energy_regularization > 0.0:
            log_scale = vector[self.group_slices["req_energy_scale"]]
            reg = reg + xp.asarray(self.training.req_energy_regularization, dtype=xp.float64) * xp.mean(log_scale * log_scale)
        if "c6_coeff" in self.group_slices and self.training.c6_regularization > 0.0:
            c6_vals = vector[self.group_slices["c6_coeff"]]
            deviation = c6_vals - xp.asarray(1.0, dtype=xp.float64)
            reg = reg + xp.asarray(self.training.c6_regularization, dtype=xp.float64) * xp.mean(deviation * deviation)
        # DISABLED(2026-03): regularization for disabled channels.
        # if "density_lap" in self.group_slices and self.training.density_lap_regularization > 0.0:
        #     lap_vals = vector[self.group_slices["density_lap"]]
        #     reg = reg + xp.asarray(self.training.density_lap_regularization, dtype=xp.float64) * xp.mean(lap_vals * lap_vals)
        # if "density_grad" in self.group_slices and self.training.density_grad_regularization > 0.0:
        #     grad_vals = vector[self.group_slices["density_grad"]]
        #     reg = reg + xp.asarray(self.training.density_grad_regularization, dtype=xp.float64) * xp.mean(grad_vals * grad_vals)
        # if "london_alt" in self.group_slices and self.training.london_alt_regularization > 0.0:
        #     lalt_vals = vector[self.group_slices["london_alt"]]
        #     reg = reg + xp.asarray(self.training.london_alt_regularization, dtype=xp.float64) * xp.mean(lalt_vals * lalt_vals)
        return reg

    def constraint_report(self, vector):
        vector = np.asarray(vector, dtype=float)
        threshold = float(self.training.constraint_bound_fraction)
        report = {}
        constraint_limited = False
        for index, (prefix, element) in enumerate(self.names):
            physical_value = float(self._to_physical(prefix, vector[index], np))
            lower, upper = self.physical_bounds[index]
            span = max(float(upper) - float(lower), 1.0e-12)
            lower_frac = (physical_value - float(lower)) / span
            upper_frac = (float(upper) - physical_value) / span
            bound_fraction = min(lower_frac, upper_frac)
            near_bound = bool(bound_fraction <= threshold)
            constraint_limited = constraint_limited or near_bound
            key = prefix if element is None else f"{prefix}:{element}"
            report[key] = {
                "value": physical_value,
                "lower": float(lower),
                "upper": float(upper),
                "distance_fraction": float(bound_fraction),
                "near_bound": near_bound,
            }
        return {
            "constraint_limited": constraint_limited,
            "threshold_fraction": threshold,
            "parameters": report,
        }


def _split_indices(n_samples, training: TrainingConfig):
    n_test = int(round(training.test_split * n_samples))
    n_val = int(round(training.validation_split * n_samples))
    n_train = max(n_samples - n_val - n_test, 1)
    rng = np.random.default_rng(training.seed + n_samples)
    indices = rng.permutation(n_samples)
    return {
        "train": indices[:n_train],
        "val": indices[n_train : n_train + n_val],
        "test": indices[n_train + n_val :],
    }


def _select_pose_batch(poses: PoseBatch, indices):
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=poses.positions[indices],
        pose_params=poses.pose_params[indices],
        site_labels=[poses.site_labels[i] for i in indices],
        metadata=dict(poses.metadata),
    )


def _select_teacher_result(teacher: TeacherResult, indices):
    return TeacherResult(
        energies=teacher.energies[indices],
        forces=teacher.forces[indices],
        metadata=dict(teacher.metadata),
    )


def _loss_on_split_numpy(model, items, params, force_weight: float, energy_weight: float, regularization: float = 0.0):
    total_loss = 0.0
    total_weight = 0.0
    compute_forces = force_weight > 0.0
    for poses, teacher in items:
        prediction = model.evaluate_batch(poses, params=params, compute_forces=compute_forces)
        energy_residual = prediction.energies - teacher.energies
        item_loss = energy_weight * np.mean(energy_residual * energy_residual)
        if compute_forces:
            force_residual = prediction.forces - teacher.forces
            item_loss += force_weight * np.mean(force_residual * force_residual)
        total_loss += item_loss * len(poses.positions)
        total_weight += len(poses.positions)
    if total_weight <= 0.0:
        return float(regularization)
    return total_loss / total_weight + float(regularization)


def _jaxify_split(model, items):
    packed = []
    for poses, teacher in items:
        fixed_charges = np.asarray(poses.adsorbate.charges, dtype=float) if poses.adsorbate.charges is not None else np.zeros(poses.adsorbate.natoms, dtype=float)
        charge_mode = 1.0 if poses.adsorbate.charges is not None and np.max(np.abs(poses.adsorbate.charges)) > 0.0 else 0.0
        packed.append(
            {
                "positions": jnp.asarray(np.asarray(poses.positions, dtype=float), dtype=jnp.float64),
                "species_idx": jnp.asarray(model.make_species_indices(poses.adsorbate.symbols), dtype=jnp.int32),
                "fixed_charges": jnp.asarray(fixed_charges, dtype=jnp.float64),
                "charge_mode": jnp.asarray(charge_mode, dtype=jnp.float64),
                "energies": jnp.asarray(np.asarray(teacher.energies, dtype=float), dtype=jnp.float64),
                "forces": jnp.asarray(np.asarray(teacher.forces, dtype=float), dtype=jnp.float64),
                "n_samples": float(len(poses.positions)),
            }
        )
    return tuple(packed)


def _loss_on_split_jax(model, layout, vector, items, force_weight: float, energy_weight: float):
    params_tree = layout.vector_to_tree(vector, use_jax=True)
    total_loss = jnp.array(0.0, dtype=jnp.float64)
    total_weight = jnp.array(0.0, dtype=jnp.float64)
    for item in items:
        item_loss = model.loss_for_batch_jax(
            positions_batch=item["positions"],
            species_idx=item["species_idx"],
            fixed_charges=item["fixed_charges"],
            charge_mode=item["charge_mode"],
            params_tree=params_tree,
            energies_true=item["energies"],
            forces_true=item["forces"],
            force_weight=force_weight,
            energy_weight=energy_weight,
        )
        total_loss = total_loss + item_loss * item["n_samples"]
        total_weight = total_weight + item["n_samples"]
    return total_loss / jnp.maximum(total_weight, 1.0) + layout.regularization(vector, use_jax=True)


def _prepare_split_data(datasets, training: TrainingConfig):
    split_data = {name: [] for name in ("train", "val", "test")}
    split_indices = {}
    for poses, teacher in datasets:
        splits = _split_indices(len(poses.positions), training)
        split_indices[poses.adsorbate.name] = {key: indices.tolist() for key, indices in splits.items()}
        for split_name, indices in splits.items():
            if len(indices) == 0:
                continue
            split_data[split_name].append((_select_pose_batch(poses, indices), _select_teacher_result(teacher, indices)))
    return split_data, split_indices


def _fit_with_scipy(density, datasets, model, training, layout, split_data, split_indices):
    history = []
    x0 = layout.pack(layout.base)

    def objective(vector):
        params = layout.unpack(vector)
        return _loss_on_split_numpy(
            model,
            split_data["train"],
            params,
            force_weight=training.force_weight,
            energy_weight=training.energy_weight,
            regularization=float(layout.regularization(vector, use_jax=False)),
        )

    def callback(vector):
        params = layout.unpack(vector)
        train_loss = _loss_on_split_numpy(
            model,
            split_data["train"],
            params,
            training.force_weight,
            training.energy_weight,
            regularization=float(layout.regularization(vector, use_jax=False)),
        )
        val_loss = train_loss
        if split_data["val"]:
            val_loss = _loss_on_split_numpy(
                model,
                split_data["val"],
                params,
                training.force_weight,
                training.energy_weight,
                regularization=float(layout.regularization(vector, use_jax=False)),
            )
        history.append({"step": len(history), "train_loss": float(train_loss), "val_loss": float(val_loss)})

    result = optimize.minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=layout.bounds,
        callback=callback,
        options={"maxiter": training.max_steps, "maxfun": max(200, training.max_steps * 50)},
    )
    final_params = layout.unpack(result.x)
    constraint_report = layout.constraint_report(result.x)
    metrics = {
        "train_loss": float(_loss_on_split_numpy(model, split_data["train"], final_params, training.force_weight, training.energy_weight, regularization=float(layout.regularization(result.x, use_jax=False)))),
        "val_loss": float(_loss_on_split_numpy(model, split_data["val"], final_params, training.force_weight, training.energy_weight, regularization=float(layout.regularization(result.x, use_jax=False)))) if split_data["val"] else None,
        "test_loss": float(_loss_on_split_numpy(model, split_data["test"], final_params, training.force_weight, training.energy_weight, regularization=float(layout.regularization(result.x, use_jax=False)))) if split_data["test"] else None,
        "regularization": float(layout.regularization(result.x, use_jax=False)),
        "split_indices": split_indices,
        "backend": "scipy_lbfgsb",
        "constraint_report": constraint_report,
    }
    return FitResult(
        params=final_params,
        history=history,
        metrics=metrics,
        optimizer_result={"success": bool(result.success), "message": result.message, "nit": int(result.nit)},
    )


def _fit_with_optax(density, datasets, model, training, layout, split_data, split_indices):
    train_items = _jaxify_split(model, split_data["train"])
    val_items = _jaxify_split(model, split_data["val"])
    test_items = _jaxify_split(model, split_data["test"])

    lower = jnp.asarray(layout.lower, dtype=jnp.float64)
    upper = jnp.asarray(layout.upper, dtype=jnp.float64)
    vector = jnp.asarray(layout.pack(layout.base), dtype=jnp.float64)
    best_vector = vector
    optimizer = optax.adam(training.learning_rate)
    opt_state = optimizer.init(vector)
    patience_limit = int(getattr(training, "early_stopping_patience", 40))
    patience = 0
    best_val = np.inf
    history = []

    def train_loss_raw(current_vector):
        return _loss_on_split_jax(model, layout, current_vector, train_items, training.force_weight, training.energy_weight)

    def eval_loss_raw(current_vector, items):
        return _loss_on_split_jax(model, layout, current_vector, items, training.force_weight, training.energy_weight)

    train_loss_fn = jax.jit(train_loss_raw)
    if val_items:
        val_loss_fn = jax.jit(lambda current_vector: eval_loss_raw(current_vector, val_items))
    else:
        val_loss_fn = None
    if test_items:
        test_loss_fn = jax.jit(lambda current_vector: eval_loss_raw(current_vector, test_items))
    else:
        test_loss_fn = None

    @jax.jit
    def step(current_vector, current_state):
        loss, grads = jax.value_and_grad(train_loss_raw)(current_vector)
        updates, current_state = optimizer.update(grads, current_state, current_vector)
        current_vector = optax.apply_updates(current_vector, updates)
        current_vector = jnp.clip(current_vector, lower, upper)
        return current_vector, current_state, loss

    for step_index in range(int(training.max_steps)):
        vector, opt_state, train_loss = step(vector, opt_state)
        train_value = float(train_loss)
        if val_loss_fn is None:
            val_value = train_value
        else:
            val_value = float(val_loss_fn(vector))
        history.append({"step": step_index, "train_loss": train_value, "val_loss": val_value})
        if val_value < best_val - 1.0e-10:
            best_val = val_value
            best_vector = vector
            patience = 0
        else:
            patience += 1
            if patience >= patience_limit:
                break

    final_vector = best_vector
    final_params = layout.unpack(as_numpy(final_vector))
    final_vector_np = as_numpy(final_vector)
    constraint_report = layout.constraint_report(final_vector_np)
    metrics = {
        "train_loss": float(train_loss_fn(final_vector)),
        "val_loss": float(val_loss_fn(final_vector)) if val_loss_fn is not None else None,
        "test_loss": float(test_loss_fn(final_vector)) if test_loss_fn is not None else None,
        "regularization": float(layout.regularization(final_vector_np, use_jax=False)),
        "split_indices": split_indices,
        "backend": "optax_adam_jax",
        "constraint_report": constraint_report,
    }
    return FitResult(
        params=final_params,
        history=history,
        metrics=metrics,
        optimizer_result={
            "success": True,
            "message": "completed",
            "nit": len(history),
            "best_val_loss": None if np.isinf(best_val) else float(best_val),
            "early_stopped": len(history) < int(training.max_steps),
        },
    )


def fit_hybrid_parameters(
    density: DensityData,
    datasets: list[tuple[PoseBatch, TeacherResult]],
    model: HybridGridFFModel | None = None,
    training: TrainingConfig | None = None,
    initial_params: HybridParameters | None = None,
):
    training = training or TrainingConfig()
    elements = tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols}))
    z_image = datasets[0][0].metadata["z_image"]
    base_params = default_hybrid_parameters(elements, z_image=z_image)
    if initial_params is not None:
        for element in elements:
            if element in initial_params.pauli:
                base_params.pauli[element] = float(initial_params.pauli[element])
            if element in initial_params.london:
                base_params.london[element] = float(initial_params.london[element])
            if element in initial_params.reactive:
                base_params.reactive[element] = float(initial_params.reactive[element])
            if element in initial_params.static_charge:
                base_params.static_charge[element] = float(initial_params.static_charge[element])
            if element in initial_params.c6_coeff:
                base_params.c6_coeff[element] = float(initial_params.c6_coeff[element])
            if element in initial_params.density_lap:
                base_params.density_lap[element] = float(initial_params.density_lap[element])
            if element in initial_params.density_grad:
                base_params.density_grad[element] = float(initial_params.density_grad[element])
            if element in initial_params.london_alt:
                base_params.london_alt[element] = float(initial_params.london_alt[element])
            if element in initial_params.req_radius_offset:
                base_params.req_radius_offset[element] = float(initial_params.req_radius_offset[element])
            if element in initial_params.req_energy_scale:
                base_params.req_energy_scale[element] = float(initial_params.req_energy_scale[element])
            if element in initial_params.chi:
                base_params.chi[element] = float(initial_params.chi[element])
            if element in initial_params.hardness:
                base_params.hardness[element] = float(initial_params.hardness[element])
        base_params.image_scale = float(initial_params.image_scale)
        base_params.image_plane = float(initial_params.image_plane)
        base_params.image_damping = float(initial_params.image_damping)
        base_params.sample_shift_z = float(initial_params.sample_shift_z)
        base_params.coulomb_shift_z = float(initial_params.coulomb_shift_z)
        base_params.reservoir_chi = float(initial_params.reservoir_chi)
        base_params.reservoir_hardness = float(initial_params.reservoir_hardness)
        base_params.total_charge = float(initial_params.total_charge)
        base_params.energy_offset = float(initial_params.energy_offset)
    model = model or HybridGridFFModel(density=density, reactive_elements=elements, prefer_jax=True)
    layout = ParameterLayout(
        elements,
        base_params,
        training=training,
        z_ref=datasets[0][0].metadata["z_ref"],
        use_req_plq=model.hybrid_config.use_req_plq,
    )
    split_data, split_indices = _prepare_split_data(datasets, training)

    if model.use_jax and HAS_OPTAX:
        return _fit_with_optax(density, datasets, model, training, layout, split_data, split_indices)
    return _fit_with_scipy(density, datasets, model, training, layout, split_data, split_indices)


def _fit_c6_lstsq(
    model: HybridGridFFModel,
    datasets: list[tuple[PoseBatch, TeacherResult]],
    stage1_params: HybridParameters,
    use_raw_r6: bool = False,
    verbose: bool = True,
) -> FitResult:
    """Fit C₆ coefficients and energy offset via linear least-squares.

    The C₆ correction is linear in the coefficients:
        E_correction = Σ_el c6_coeff_el × basis_el + offset

    lstsq gives the exact MSE-optimal solution in one step, with unconstrained
    coefficients (can be negative) and an energy offset term.

    The design matrix is built by sampling whatever dispersion grid is installed
    in the model (TS-damped by default, or raw 1/r⁶ after install_raw_r6_grid()).

    Args:
        use_raw_r6: Label flag for logging; the actual grid swap must be done
            before calling this function via model.install_raw_r6_grid().
    """
    from copy import deepcopy

    basis_label = "raw_r6" if use_raw_r6 else "ts_damped"

    # Gather all training poses and teacher energies
    all_poses = []
    all_teacher_e = []
    for poses, teacher in datasets:
        for k in range(len(poses.positions)):
            all_poses.append(poses.positions[k])
            all_teacher_e.append(float(teacher.energies[k]))
    all_teacher_e = np.asarray(all_teacher_e, dtype=float)
    n_poses = len(all_poses)

    # Evaluate baseline energies with c6_coeff=0, energy_offset=0
    baseline_energies = np.zeros(n_poses, dtype=float)
    idx = 0
    for poses, _ in datasets:
        result = model.evaluate_batch(poses, params=stage1_params, compute_forces=False)
        for k in range(len(poses.positions)):
            baseline_energies[idx] = float(result.energies[k])
            idx += 1

    # Build design matrix
    elements = model.parameter_elements
    n_el = len(elements)
    X = np.zeros((n_poses, n_el), dtype=float)
    idx = 0
    if verbose:
        print(f"[lstsq] Building {basis_label} design matrix for {n_poses} poses ...", flush=True)
    # Always use compute_dispersion_design_matrix — it samples whatever grid
    # is installed (TS-damped or raw r6 via install_raw_r6_grid()).
    # This ensures lstsq basis matches evaluation basis exactly.
    for poses, _ in datasets:
        X_batch = model.compute_dispersion_design_matrix(poses, stage1_params)
        n_batch = len(poses.positions)
        X[idx:idx + n_batch] = X_batch
        idx += n_batch

    # Residual: what the baseline PLQ misses
    residual = all_teacher_e - baseline_energies

    # Augment with offset column: [X | ones]
    X_aug = np.column_stack([X, np.ones(n_poses)])

    # Solve: residual ≈ X_aug @ D
    D, _, _, _ = np.linalg.lstsq(X_aug, residual, rcond=None)

    # Extract coefficients
    c6_coeff_values = D[:n_el]
    energy_offset = float(D[n_el])

    if verbose:
        print(f"[lstsq/{basis_label}] C₆ coefficients: {dict(zip(elements, c6_coeff_values))}", flush=True)
        print(f"[lstsq/{basis_label}] energy offset: {energy_offset * 1000:.1f} meV", flush=True)
        corrected = baseline_energies + X @ c6_coeff_values + energy_offset
        rmse_before = np.sqrt(np.mean((baseline_energies - all_teacher_e) ** 2)) * 1000
        rmse_after = np.sqrt(np.mean((corrected - all_teacher_e) ** 2)) * 1000
        print(f"[lstsq/{basis_label}] train RMSE: {rmse_before:.1f} → {rmse_after:.1f} meV", flush=True)

    # Build result params
    result_params = deepcopy(stage1_params)
    for i, el in enumerate(elements):
        result_params.c6_coeff[el] = float(c6_coeff_values[i])
    result_params.energy_offset = energy_offset

    corrected = baseline_energies + X @ c6_coeff_values + energy_offset
    train_loss_after = float(np.mean((corrected - all_teacher_e) ** 2))

    return FitResult(
        params=result_params,
        history=[{"step": 0, "train_loss": train_loss_after, "val_loss": train_loss_after, "stage": 2, "method": f"lstsq_{basis_label}"}],
        metrics={
            "train_loss": train_loss_after,
            "lstsq_residual_norm": float(np.linalg.norm(residual - X_aug @ D)),
            "backend": f"lstsq_{basis_label}",
        },
        optimizer_result={"success": True, "method": f"lstsq_{basis_label}", "n_params": n_el + 1},
    )


def fit_two_stage_c6(
    density: DensityData,
    datasets: list[tuple[PoseBatch, TeacherResult]],
    model: HybridGridFFModel | None = None,
    training: TrainingConfig | None = None,
    initial_params: HybridParameters | None = None,
    lstsq_datasets: list[tuple[PoseBatch, TeacherResult]] | None = None,
    use_raw_r6: bool = False,
    verbose: bool = True,
) -> FitResult:
    """Two-stage C₆ optimization: fit REQ first, then freeze REQ and fit C₆.

    Args:
        lstsq_datasets: Optional larger dataset for the C₆ lstsq fit (Stage 2).
            If None, uses the same ``datasets`` as Stage 1. Since Stage 2 is a
            5-parameter linear fit, using all available data (including holdout)
            reduces overfitting risk and maximises accuracy.
        use_raw_r6: If True, Stage 2 uses raw -Σ 1/|r-r_j|⁶ basis (no TS damping,
            no c6_pair factor) matching the proof-of-concept. Default False uses
            TS-damped dispersion grids.

    Physics rationale: density-based PLQ captures short-range Pauli/exchange
    (exp(-κz) decay) but cannot reproduce the long-range 1/r⁶ vdW-surf
    dispersion tail present in MAD-SURF teacher data (trained on PBE+vdW^surf).
    The C₆ grid channel adds the missing 1/r⁶ shape. Joint REQ+C₆ optimization
    fails due to parameter coupling (202 meV vs 37 meV baseline). Decoupling
    into stages avoids bad minima.

    Stage 1: Fit REQ with c6_coeff=0.0 (zero C₆ contribution = PLQ-only baseline)
    Stage 2: Freeze all REQ, fit only c6_coeff from initial value 0.0
    Stage 3 (optional): Gentle joint refinement with very small learning rate
    """
    from copy import deepcopy
    training = training or TrainingConfig()

    # Build initial params with c6_coeff zeroed out so Stage 1 is PLQ-only.
    # default_hybrid_parameters sets c6_coeff=1.0; we must override to 0.0.
    elements = tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols}))
    z_image = datasets[0][0].metadata["z_image"]
    stage1_init = initial_params if initial_params is not None else default_hybrid_parameters(elements, z_image=z_image)
    stage1_init = deepcopy(stage1_init)
    for el in list(stage1_init.c6_coeff.keys()):
        stage1_init.c6_coeff[el] = 0.0

    # ---- STAGE 1: REQ-only (density PLQ, c6_coeff frozen at 0.0) ----
    stage1_training = deepcopy(training)
    stage1_training.fit_c6_coeff = False
    if verbose:
        print("[two-stage] Stage 1: fitting REQ parameters (C₆ disabled)...", flush=True)

    stage1_result = fit_hybrid_parameters(
        density=density,
        datasets=datasets,
        model=model,
        training=stage1_training,
        initial_params=stage1_init,
    )
    if verbose:
        s1_loss = stage1_result.metrics.get("train_loss", "?")
        print(f"[two-stage] Stage 1 done: train_loss={s1_loss}", flush=True)

    # ---- STAGE 2: C₆ via lstsq (exact linear solution) ----
    if model is None:
        raise ValueError("fit_two_stage_c6 requires an explicit model (cannot be None)")
    if verbose:
        print("[two-stage] Stage 2: fitting C₆ coefficients via lstsq (REQ frozen)...", flush=True)

    stage2_result = _fit_c6_lstsq(
        model=model,
        datasets=lstsq_datasets if lstsq_datasets is not None else datasets,
        stage1_params=stage1_result.params,
        use_raw_r6=use_raw_r6,
        verbose=verbose,
    )
    if verbose:
        print("[two-stage] Stage 2 done.", flush=True)

    final_result = stage2_result

    # ---- STAGE 3 (optional): gentle joint refinement ----
    if training.stage3_refine:
        stage3_training = deepcopy(training)
        stage3_training.fit_c6_coeff = True
        stage3_training.max_steps = training.stage3_max_steps
        stage3_training.learning_rate = training.stage3_learning_rate
        if verbose:
            print("[two-stage] Stage 3: gentle joint refinement...", flush=True)

        final_result = fit_hybrid_parameters(
            density=density,
            datasets=datasets,
            model=model,
            training=stage3_training,
            initial_params=stage2_result.params,
        )
        if verbose:
            s3_loss = final_result.metrics.get("train_loss", "?")
            print(f"[two-stage] Stage 3 done: train_loss={s3_loss}", flush=True)

    # Combine histories from all stages
    combined_history = []
    offset = 0
    all_results = [stage1_result, stage2_result]
    if training.stage3_refine:
        all_results.append(final_result)
    for stage_num, result in enumerate(all_results, start=1):
        for entry in result.history:
            combined_entry = dict(entry)
            combined_entry["stage"] = stage_num
            combined_entry["global_step"] = offset + entry.get("step", 0)
            combined_history.append(combined_entry)
        offset += len(result.history)

    return FitResult(
        params=final_result.params,
        history=combined_history,
        metrics={
            **final_result.metrics,
            "stage1_train_loss": stage1_result.metrics.get("train_loss"),
            "stage2_train_loss": stage2_result.metrics.get("train_loss"),
            "two_stage": True,
        },
        optimizer_result=final_result.optimizer_result,
    )
