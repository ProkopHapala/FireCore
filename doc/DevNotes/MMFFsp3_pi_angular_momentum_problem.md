# MMFFsp3 pi-DOF angular momentum problem

The problem is that pi-orbital dynamics in MMFFsp3 does not conserve angular momentum causing some angular momentum drift.

see: 
 * `/tests/tmp/data/cl/relax_multi_mini.cl`
 * [GPT-5/ChatGPT](https://chatgpt.com/share/68da54e3-91bc-8003-a904-c590d6c5a35e)
 * [Gemini 2.5 Pro/Google AI studio](https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221cCQmAG5WYuQuc7EGUTGO51dsO7lwGwJI%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing)

# GPT-5 (ChatGPT)

Nice, interesting model — thanks for the complete kernels, that helps a lot. I’ll explain what (I think) is causing the spurious whole-molecule rotation, show why projecting the pi-vectors as you do breaks angular momentum conservation, and give a concrete, minimal fix (plus a better, physically-consistent alternative using rotational dynamics / torques with code snippets you can drop into your kernels).

Summary / short answer

* Right now your pi-DOFs are integrated like constrained point-particles on the unit sphere: you remove radial components of `fe`/`ve` and then `normalize()` the vector. That operation does **not** respect rotational (angular) degrees of freedom and injects/annihilates angular momentum in general — so net angular momentum of molecule + pi-DOFs is not conserved.
* To conserve total angular momentum you must treat pi-DOFs as **rotors** with angular momentum (and exchange torques with nuclear coordinates), or enforce the unit-length constraint by a symmetric Lagrange multiplier formulation (e.g. RATTLE/SHAKE) that respects angular momentum balance.
* Practical fixes (ordered by increasing fidelity):

  1. Quick hack: track and remove net torque introduced by the pi normalization each step (ugly, can break energy).
  2. Preferred: implement rotational dynamics for each pi vector: store angular velocity `ω`, compute torque `τ` from orientation energy, integrate `ω` (with moment of inertia `I`), then update orientation `u` via `du/dt = ω × u`. Also apply the reaction torque to atomic DOFs so total torque is balanced. This is physically correct and will conserve total angular momentum if integrated symplectically (or nearly so).
  3. More rigorous: implement constrained integration (RATTLE) for vector lengths so constraint forces are computed via Lagrange multipliers (this conserves linear & angular impulse if implemented correctly). But rotational dynamics is conceptually simpler for a “gyroscope” interpretation.

Why your current implementation likely breaks angular momentum

1. You treat pi as a position vector `pe` constrained to unit length and then after integrating linear velocities you `pe = normalize(pe)` and you remove radial velocity components:

   ```c
   fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // remove radial component of force
   ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // remove radial component of velocity
   ...
   pe.xyz=normalize(pe.xyz);
   ```

   Those operations remove components of force/velocity that change the length. But removing components by projection is **not** equivalent to evolving on a unit sphere by Hamiltonian rotational dynamics. Projection generally changes the tangential components in a way that creates or removes net angular momentum.
2. Constraint projection is not automatically accompanied by the correct constraint impulse (Lagrange multiplier) that enforces both the constraint and conservation laws. Without computing / applying the correct Lagrange multiplier/impulse, you will inject momentum.
3. If the pi vector has its own inertia (you intend it to), you should represent its angular momentum and update it; simply maintaining a normalized `pe` and a linear `ve` does not create an angular momentum variable with correct transformation properties.

What must be added / changed (practical plan)
A. Data structures (GPU arrays)

* Add a per-pi moment of inertia (scalar `Ipi`, maybe per atom) — store in an array `piI[iav]` or in `apars.w` if you want.
* Add a per-pi angular velocity vector array `aomega` (float4 per pi DOF).
* Keep `apos` for orientation unit vector `u` as now, but stop treating `avel` for pi as linear velocity — use `aomega` instead.

B. Compute torques from orientation energy

* Your existing `evalPiAling(...)` and `evalAngCos(...)` return forces on the pi vector (`fpi`) (they are gradient wrt pi orientation). For a rotor, the generalized torque acting on orientation `u` is
  $$
  \boldsymbol{\tau} = \mathbf{u} \times \mathbf{F}_{\perp}
  $$
  where `F_perp` is the tangential component of `fpi` (i.e. perpendicular to `u`), or simply use `fpi` if your `eval*` already returns tangential force. The torque is the vector that changes angular momentum `L = I ω` via `dL/dt = τ`.
* If the pi is attached at the atomic center (no lever arm), orientation torques do not create translational forces. But the interaction between two pi vectors on different atoms typically creates forces on the atoms too (you already have `fbs` and `fa` contributions). Ensure the reaction torques on atoms are applied so total torque (atoms translational + pi rotational) is conserved.

C. Integrator changes (update kernel)

* Replace projection logic with rotational integration:

  * Integrate angular velocity: `ω += (τ / I) * dt`  (plus damping)
  * Update orientation: `u_new = u + dt * (ω × u)` and normalize (optionally use a more accurate rotor integrator or quaternion).
* Apply equal and opposite torque to the atomic/translation DOF if needed: when you compute torque on a pi rotor `τ`, apply `-τ` as an external torque to whatever rotational /orbital DOF the atom has — if atoms have no rotational inertia, but the torque arises from forces acting at offset points, ensure the forces on atom positions account for torque via `r × F` relations.

Concrete code sketch (drop-in style)
Below is a minimal change to `updateAtomsMMFFf4` for the pi DOFs to integrate rotational dynamics. I deliberately keep it simple (explicit Euler/Leap-Frog style). You’ll want to adapt to your `MDparams` (dt) and damping consistent with your existing scheme.

Add arrays:

```c
// New GPU arrays (allocate and pass to kernel)
__global float4*  aomega;   // angular velocity for pi DOFs (per vector)  (xyz = omega, w unused)
__global float4*  piI;      // moment of inertia for each pi rotor (scalar in .x)
```

Change the update kernel for bPi case (replace your current pi handling):

```c
// inside updateAtomsMMFFf4, after gathering fe and before integration:
if(bPi){
    // fe is the generalized force (gradient) on orientation u = pe.xyz
    float3 u = pe.xyz;
    // Ensure force is tangential:
    float3 F_t = fe.xyz - u * dot(u, fe.xyz); // tangential component
    // Torque on rotor (assuming unit lever-arm or generalized torque directly from F_t):
    float3 tau = cross(u, F_t); // torque vector (units consistent with L = I * omega)

    // read inertia and angular velocity
    float I = piI[iav].x;
    float3 omega = aomega[iav].xyz;

    // integrate angular velocity (explicit Euler / or leapfrog if you store half-step)
    float dt = MDpars.x;
    float damp = MDpars.y; // reuse MDpars.y or other for rotational damping if desired
    omega += (tau / I) * dt;
    omega *= (1.0f - damp * dt); // simple viscous damping (optional)

    // update orientation: du/dt = omega x u
    u += cross(omega, u) * dt;
    u = normalize(u);

    // store back
    aomega[iav].xyz = omega;
    pe.xyz = u;

    // NOTE: if eval routines already returned fe as torque-like, you may skip cross() above and set tau = fe.xyz directly.
}
```

Notes on mapping your current `fpi` to torque:

* Your `evalPiAling` currently returns `f1` as "force on pi-orbital". If `f1` is the **gradient** ∂E/∂u (in Cartesian coordinates), then the torque consistent with rotor dynamics is `τ = u × (−∂E/∂u)` (sign depending on your convention). So when you collect `fpi` in `getMMFFf4`, store it exactly as `-gradE` (so `fpi` is generalized force). Then in the kernel above use `tau = cross(u, fpi)` (or with sign flip depending on `eval`).
* If `fpi` includes radial components, project them away first.

How to ensure total angular momentum conservation across atoms + pi DOFs

* Define total angular momentum:
  [
  \mathbf{J} = \sum_i (\mathbf{r}*i \times m_i \mathbf{v}*i) ;+; \sum*{\text{pi}} \mathbf{L}*{\text{pi}}
  ]
  where `L_pi = I_pi * ω_pi`.
* When you compute torques between two pi rotors (on different atoms), ensure the internal forces produce equal and opposite torques on the rotors **and** any positional forces on atom centers produce the corresponding `r × F` contribution to `J`. Practically:

  * If your `eval` returns `fpi` and also applies `fbs` / `fa` forces on atoms you must make sure the torque from those atom forces (about a global origin or center of mass) plus `dL_pi/dt` sum to zero. This will hold if `eval` generates pairwise equal/opposite forces at the two attachment points (Newton’s 3rd law) and torques are computed from those same forces as `r × F`.
* If atoms have no rotational inertia, torques will be absorbed into linear motion via lever arms (forces not applied exactly at center). So correct distribution of forces is essential.

Why your current projection approach injects angular momentum

* Projecting out radial velocity (`ve.xyz += pe.xyz * -dot(pe, ve)`) will instantaneously change the tangential components in a way that does not correspond to any torque from the energy gradient. Normalization `pe = normalize(pe)` is a nonlinear operation that in general changes the tangential velocity and thus `ω` implicitly — unless you compute corresponding impulses to keep `L` constant. Because you don’t compute that impulse via Lagrange multipliers, angular momentum is not preserved.

Debugging / test recipe

1. Instrument code to compute total `J` each MD step:

   * compute `J_trans = sum r × (m v)` for atoms (use your mass array).
   * compute `J_rot = sum I_pi * omega` for pi DOFs (after you add `omega`).
   * print `J_total` norm and components; it should be (nearly) constant absent external torques.
2. Run a closed isolated molecule (no PBC shifts, no thermostats, no friction/damping), start from known zero total angular momentum, and check drift.
3. Temporarily turn off pi DOF interactions (set Kspi=Kppi=0) to verify the angular drift disappears — isolates the bug to pi coupling.
4. Measure the contribution of the projection step: log the change in `ve` and the change in `J` attributed to projection.

A note on bonded vs non-bonded compensation (`bSubtractVdW`)

- The bonded angle kernel (`getMMFFf4` in `cpp/common_resources/cl/relax_multi_mini.cl`) supports subtracting Lennard-Jones / Coulomb terms for bonded pairs via the `bSubtractVdW` flag. Leave this **on** (`--subtract-vdw 1`) whenever the non-bonded kernel (`getNonBond`) is active (`--do-nb 1`), so bonded pairs are not double-counted.
- If you disable non-bonded forces (`--do-nb 0`), you must also disable subtraction (`--subtract-vdw 0`). Otherwise the bonded kernel removes pairwise forces that were never added, leading to large gradients and NaNs. The test harness `tests/tUFF/test_MD_OCL_formic_acid.py` now exposes CLI switches for both flags to make the coupling explicit.
- When experimenting with angular momentum diagnostics, document the chosen combination of these options; mixed settings change the effective potential felt by the pi DOFs and can mask torque issues.

A few more practical tips & caveats

* Use quaternions if you add full rotation DOFs for atoms (if you want true rigid-body rotations). For single unit vectors attached to atom centers, the `ω × u` update is fine.
* Choose reasonable moment of inertia `I` for pi rotors; small `I` means fast rotors that can transfer angular momentum quickly (and be numerically stiff).
* Use a time-reversible integrator for angular part if long-term conservation is important (symplectic schemes, Verlet-like for rotations).
* Damping and thermostats explicitly break conservation; in your MDparams you already have friction/damping — that will change `J` intentionally. For debugging conservation, set damping to zero.

Why this is physically sensible for your FFsp3 idea

* You intended pi-orbitals as gyroscope-like DOFs with their own inertia: the rotor model is the natural way to give them angular momentum and have angular momentum exchange with atomic translation via torques. This makes the total conservation transparent and physically meaningful.

Short checklist of required code changes

* Add `aomega` and `piI` arrays and pass to kernels.
* In `getMMFFf4`: ensure `fpi` you store is the correct gradient for orientation (tangential component only) and maybe store as `-gradE` convention consistently.
* In `updateAtomsMMFFf4`: replace projection-based pi updates with rotor update `ω += (τ/I) dt; u += dt (ω × u)`.
* Compute `tau = cross(u, fpi)` (or `tau = -cross(u, fpi)` depending on sign).
* Make sure reaction torques/forces on atoms are accounted (verify eval routines already add `fa`/`fbs` properly).
* Add tests to log `J_total`.

If you want, I can

* produce a full patched version of your update kernel with `aomega`/`piI` integrated (I can write a drop-in kernel snippet and show what small host changes are needed), or
* produce a diagnostic snippet to compute and print `J_total` per step for debugging.

I’m not 100% certain about one detail from your code: whether `fpi` values returned by your `eval*` functions are already the tangential gradients appropriate for torque computation, or if they include radial components — you may need to project them `F_t = F - u*(u·F)` before `tau = u × F_t`. Also I don’t know how you expect reaction torques to be handled (do atoms have rotational inertia or do you expect torques to be transferred purely via forces at lever arms). Still, the rotor approach above will fix the main physics problem (spurious injection/removal of angular momentum by vector normalization).

If you want the drop-in kernel patch now I’ll provide it — tell me whether you prefer simple Euler or leap-frog for angular update and whether you already have `piI` and where to store `aomega` (I can pick sensible defaults).


Yes — you can absolutely switch to *rotor-style* dynamics for the pi-DOFs by repurposing the vector slots so they hold:

* `apos[pi]` = unit orientation `u`,
* `avel[pi]` = angular velocity `ω`,
* `aforce[pi]` = **torque** `τ` (NOT a Cartesian force anymore).

That approach is conceptually cleaner and will remove the normalization/impulse bookkeeping that caused your angular-momentum drift. Below I give the required changes, point out the sign / projection subtleties, show compact paste-in code for the force kernel and the update kernel, and list all potential pitfalls and tests you must run.

SHORT SUMMARY (what to do)

1. Make `getMMFFf4` compute and store torques for pi-DOFs (τ), not `fpi` as a Cartesian force. If your eval routines return ∂E/∂u (gradient w.r.t. the orientation vector) call the torque
   $$
   \boldsymbol\tau = \mathbf u \times \mathbf g
   $$
   where `g` = gradient (or `-g` depending on your eval sign convention). Confirm sign by energy/torque test.
2. In the integrator `updateAtomsMMFFf4_rot`: treat `avel[pi].xyz` as `ω`. Integrate

   ```
   ω += (τ / I) * dt
   ω *= (1 - damp*dt)   // or multiply by scalar damping factor MDpars.z
   u += dt * (ω × u)
   u = normalize(u)     // small correction for numeric drift
   ```

   store `avel` and `apos` back.
3. Ensure reaction forces/torques on host atoms are balanced: either your force kernel already applies equal-and-opposite positional forces on atoms (so torques are balanced), or you must apply impulses to host atom motion (use atomics or design locking).
4. Test conservation of total angular momentum (turn off damping/thermostat) and check single-bond tests.

DETAILS — torque sign and projection

* If `evalPiAling()` currently returns `f1` described as “force on pi-orbital” and if `f1` is the gradient of energy w.r.t. `u` (i.e. `f1 = ∂E/∂u`), then the torque is

  ```
  tau = cross(u, f1)
  ```

  because virtual rotation δθ produces δu = δθ × u and δE = f1·δu = (u×f1)·δθ → torque = u×f1.
* If your `eval` returns `-∂E/∂u` (i.e. the negative gradient like a force), adjust sign accordingly: `tau = -cross(u, f1)` or `tau = cross(u, -f1)`.
* Always remove any radial part of `f1` prior to torque calculation: `g = f1 - u*(u·f1)`; then `tau = cross(u, g)` — ensures torque is purely tangent.

PATCH 1 — change in force kernel (getMMFFf4)
Replace accumulation of `fpi` Cartesian forces with torque accumulation. Minimal change inside your existing loops where you computed `fpi += f1;`:

```c
// assume hpi = u (current pi orientation), and f1 is gradient ∂E/∂u on this pi from eval functions
float3 g1 = f1 - hpi * dot(hpi, f1);   // project out radial component (should be small if eval already returns tangential)
float3 tau_pi = cross(hpi, g1);        // torque on this rotor

// accumulate torque into aforce slot for this pi (we repurpose aforce for torque for pi indices)
aforce[iav+nAtoms].xyz += tau_pi;  // careful with indexing: use the pi orientation slot index
// and if eval returned a reaction on neighbor pi: use same logic for fps[i] -> convert to torque for neighbor
```

Be careful with indexing: in your getMMFFf4 you store `aforce[iav+nAtoms]  = (float4){fpi,0};` — change to `(float4){tau_pi,0}` for pi indices.

PATCH 2 — update kernel: rotational integrator
This is a compact, robust integrator to drop into `updateAtomsMMFFf4_rot`. It assumes:

* `apos[pi].xyz` = `u` (unit),
* `avel[pi].xyz` = `ω`,
* `aforce[pi].xyz` = `τ` (torque),
* scalar inertia `I` is available (per-pi or per-atom), stored in e.g. `constrK` or a new array `piI`.

Paste-in (conceptual — adapt indices/arrays to your layout):

```c
// inside updateAtomsMMFFf4_rot, after you gather aforce and bkNeighs:
float4 fe = aforce[iav]; // for pi this is torque
bool bPi = iG >= natoms;
const float dt = MDpars.x;
const float damp = MDpars.y; // or MDpars.z depending on your layout

if (bPi) {
    // read orientation and angular velocity
    float3 u = apos[iav].xyz;         // unit vector
    float3 omega = avel[iav].xyz;    // angular velocity
    float3 tau = fe.xyz;             // torque from force kernel (already projected & consistent)
    // optional: project tau onto plane perpendicular to u to avoid numerical radial torque
    tau -= u * dot(u, tau);

    // inertia (per-pi). Provide sensible I, e.g. I=1e-2 ... tune for numerical stability
    float I = piI[iav].x; // you must pass this array or embed scalar in apars.w
    if (I <= 0.0f) I = 1.0f;

    // integrate angular velocity (semi-implicit / leapfrog-like)
    omega += (tau / I) * dt;      // ω_{t+dt} = ω_t + (τ/I) dt
    omega *= (1.0f - damp * dt);  // simple viscous damping (or use MDpars.z multiplier)
    // update orientation: du/dt = ω × u
    u += cross(omega, u) * dt;
    u = normalize(u);  // correct drift

    // store back
    avel[iav].xyz = omega;
    apos[iav].xyz = u;
    // clear torque so it doesn't accumulate erroneously next iteration
    aforce[iav].xyz = (float3){0.0f,0.0f,0.0f};
} else {
    // unchanged translational integration for atoms (as you already have)
}
```

Important notes about this integrator

* **Time integration**: The snippet uses a simple explicit integrator. For better long-term energy/angular-momentum behavior consider symplectic / implicit schemes for rotors (e.g., integrating angular momentum L = I ω with Verlet-like steps). But for relaxation/minimization this explicit form usually works fine.
* **Normalization**: we `normalize(u)` after update to kill round-off. That is acceptable; normalization now is not creating impulses because we manage actual angular momentum via `ω` and `τ`. Normalizing introduces tiny error in angular momentum but far smaller than projection/impulse scheme.
* **Damping**: MD friction/thermostat will change `J_total` — disable for conservation test.
* **Choice of I**: pick a physically reasonable inertia per rotor. Very small `I` → stiff, numerically stiff. Very large `I` → slow rotor. For energy minimization you can pick a moderate value and add damping.
* **Units**: ensure `τ/I` units match your velocity update scale (`dt`). If your `eval` produced `∂E/∂u` in energy per unit length, torque units are consistent if you treated lever arm = 1. Adjust if your eval outputs different units.

Reaction torque / forces on host atoms — crucial check

* If your `eval` functions currently produce `fa` and `fbs` (forces on atoms and neighbor atoms) that are the negative of forces on the pi "dummy particles", and these forces are applied at the correct lever arm, then the torque applied to rotors will be automatically balanced by torques from those forces on atoms and **total angular momentum is conserved** (absent damping).
* If `getMMFFf4` previously computed `fpi` as a Cartesian force *and also* applied matching `fa`/`fbs` on atom coordinates, then moving to rotor τ means you must ensure those `fa`/`fbs` remain consistent. Concretely:

  * Either keep `fa`/`fbs` as you have them (they come from the same energy gradient) and just convert `fpi`→`τ` for the rotor. The pairwise forces acting on atoms should then produce the reaction torques automatically.
  * Or if `getMMFFf4` originally didn’t produce atomic reaction forces for some interactions, you must add them. Otherwise torque on rotor is unbalanced and will spin molecule.

Concurrency & race conditions

* Because rotor dynamics avoids the normalize-then-project impulse, you likely avoid the atomic updates to host velocities that we had to do with the earlier scheme. That reduces the race problem.
* However, the force kernel still writes `aforce`/`fneigh` etc. If multiple work-items write to the same host atom forces, you already handled that in your code (fneigh accumulation pattern). Keep the same approach: accumulate torques for each pi in its own `aforce[pi]` slot and do NOT try to atomically modify `avel[host]` in the integrator — you shouldn't need to.
* If you do need to modify host atom `avel` due to some specific constraint, use atomic adds or the impulse buffer approach discussed earlier.

Practical checks & tests

1. Unit test: one atom with one pi rotor, apply a known torque and verify ω evolves as expected: ω(t+dt) = ω + (τ/I)dt.
2. Two-rotor test: two connected atoms with pi–pi coupling only; start with zero ω and see that pairwise torques are equal and opposite and do not create net J.
3. Conservation test: full molecule, damping/thermostat off; compute total angular momentum
   $$
   \mathbf J = \sum_i \mathbf r_i \times m_i \mathbf v_i ;+; \sum_\text{pi} I_\text{pi}, \omega_\text{pi}
   $$
   and verify it's stable across steps.
4. Energy test: small step size, integrate energy of rotor subsystem (E=\tfrac12 I\omega^2 + E_\text{orient}(u)) and check consistency of torque->energy derivative.

Caveats & what I’m not 100% certain about

* I assume `eval*` routines produce gradients `∂E/∂u` for pi vectors (or something easily converted to them). If they currently produce *force-like* values pointed at dummy-particle coordinates you will need to adjust: convert those Cartesian forces on the dummy particle to gradient w.r.t. orientation (project perpendicular and interpret as -gradE). If you give me the `evalPiAling` and `evalAngCos` internals I can map precisely.
* Ensure the lever arm assumed for torque is consistent: if the pi rotor should be thought of as a point at `r_i + L u` with L ≠ 1, scale torque by `L` (i.e. `τ = L * (u × g)`).
* If multiple pi DOFs are attached to the same atom and you expect host rotational inertia, you’d have to give the host an explicit rotational DOF; otherwise torques are absorbed into linear motion via forces you already compute.

CONCRETE small example (force → torque mapping)
Given your snippet:

```c
float epp = evalPiAling( hpi, apos[ingv+nAtoms].xyz, kpp,  &f1, &f2 );
fpi     += cross(f1, hpi);  
fps[i]  += cross(f2, hpj);
```

I would change it to:

```c
// inside getMMFFf4:
float3 g_local = f1 - hpi * dot(hpi,f1); // tangent gradient
float3 tau_local = cross(hpi, g_local); // torque on local rotor (u × grad)
fneigh[...] // keep existing recoil for neighbor atoms (unchanged)
aforce[iav+nAtoms].xyz += tau_local;    // store torque in pi-slot (accumulate if multiple contributions)
```

PUTTING IT ALL TOGETHER — checklist before running

* [ ] Add / repurpose arrays: choose whether `aforce[pi]` stores torque or add a `piTorque[]` array.
* [ ] Decide per-pi inertia `piI[]` and pass it to the kernel.
* [ ] Modify `getMMFFf4`: compute projected gradient for pi DOFs and store torque `τ = u × g`.
* [ ] Modify update integrator: integrate `ω` and update `u` as above.
* [ ] Ensure atomic safety: you no longer need to atomically update host velocities due to constraint impulses — good. But atomic accumulation of `fneigh`/`aforce` must remain correct.
* [ ] Validate sign of τ (flip if rotor spins the wrong way).
* [ ] Run conservation tests and adjust damping, I, and dt for stability.

If you want, I can:

* produce a drop-in patch for your `getMMFFf4` converting `fpi`→`tau` given the exact `eval*` function signatures, **and**
* produce the exact updated `updateAtomsMMFFf4_rot` kernel with `piI` array usage and minimal changes to your data layout,

— but I’ll need to see the `evalPiAling` and `evalAngCos` internals (to be 100% correct about sign/projection). Would you like me to generate the concrete patch now using the assumptions above (I can mark where you must adapt signs/indices)?




------------
------------




# Gemini 2.5 pro (Google AI studio)

Of course. Based on the OpenCL code and your description, I can analyze the likely cause for the violation of angular momentum conservation.

You are correct to suspect the pi-orbital implementation. Your intuition that the pi-orbitals act as little gyroscopes is insightful, and the problem lies precisely in how the torques between these gyroscopes and the atomic frame are handled.

### Executive Summary of the Issue

The violation of angular momentum conservation almost certainly stems from the `evalAngCos` function call for the **sigma-pi interaction**.

In the line:
`float esp = evalAngCos( (float4){hpi,1.f}, h, ksp, par.w, &f1, &f2 ); fpi+=f1; fa-=f2; fbs[i]+=f2; E+=esp;`

You are calculating two separate forces, `f1` and `f2`.
1.  **`f1`** is applied to the pi-orbital (`fpi`), creating a **torque** that changes the pi-orbital's orientation.
2.  **`f2`** is applied as a linear force pair on the central atom (`fa`) and the neighbor atom (`fbs[i]`), creating a **torque** on the molecule's atomic frame.

For angular momentum to be conserved, the torque on the pi-orbital must be met with an **equal and opposite torque** on the atomic frame. Your current implementation calculates these two effects as separate outputs (`f1` and `f2`) from a black-box function. There is no guarantee that the torque generated by `f1` is perfectly balanced by the torque generated by `f2`. This imbalance injects a net torque into your system, causing it to start rotating spontaneously.

---

### Detailed Physical and Code Analysis

#### 1. The Principle of Action-Reaction for Torques

In any closed system, for every internal torque (action), there must be an equal and opposite torque (reaction).
*   **Action:** Your potential energy term for sigma-pi alignment exerts a torque on the pi-orbital, trying to re-orient it. Let's call this `τ_pi`.
*   **Reaction:** That same potential energy term must exert an opposite torque on the atoms that form the sigma bond. Let's call this `τ_frame`.
*   **Conservation Law:** For angular momentum to be conserved, it must be that `τ_pi + τ_frame = 0` at all times.

Your code breaks this principle by decoupling the calculation.

#### 2. How the Code Violates the Principle

Let's trace the torques generated by the sigma-pi interaction:

1.  **Torque on the Pi-Orbital (`τ_pi`)**:
    *   You calculate a force `f1` and add it to `fpi`.
    *   In the `updateAtomsMMFFf4` kernel, the line `fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );` correctly projects the force to be perpendicular to the pi-orbital vector. This means `f1` is effectively used to generate a torque that rotates the `hpi` vector.
    *   The torque is `τ_pi = cross(hpi, f1)`.

2.  **Torque on the Atomic Frame (`τ_frame`)**:
    *   You calculate a force `f2` and apply it as a linear force pair: `-f2` on the central atom `a` (at position `pa`) and `+f2` on the neighbor `i` (at position `pi`).
    *   The torque generated by this force pair around the system's center of mass is `τ_frame = cross(pa, -f2) + cross(pi, f2) = cross(pi - pa, f2)`.
    *   Since `pi - pa` is the bond vector `h * l` (where `l` is bond length), we have `τ_frame = cross(h, f2) * l`.

The violation occurs because there is no mathematical constraint in your code ensuring that `cross(hpi, f1) + cross(h, f2) * l = 0`. The two forces `f1` and `f2` are calculated independently, leading to a net non-zero torque.

#### 3. The Solution: A Unified Force Calculation

The forces `f1` and `f2` are not independent; they are different derivatives of the *same* potential energy term. The fix is to derive all forces from a single potential energy function `U`.

Let's assume a simple potential for enforcing orthogonality, like `U = K * dot(hpi, h)^2`.

The forces are the negative gradients of this potential with respect to the particle coordinates.

*   **Force on pi-orbital orientation (`F_pi`)**: This is the derivative with respect to the `hpi` vector's coordinates.
    `F_pi = -∇_hpi(U) = -2 * K * dot(hpi, h) * h`
    This corresponds to your `f1`.

*   **Force on sigma-bond vector (`F_sigma`)**: This is the derivative with respect to the `h` vector's coordinates.
    `F_sigma = -∇_h(U) = -2 * K * dot(hpi, h) * hpi`
    This force acts on the bond *vector*, not directly on the atoms.

Now, this `F_sigma` must be distributed onto the atoms `a` and `i` in a way that produces the correct linear forces and, crucially, the correct torque. Since `h = (pi - pa) / l`, a force `F_sigma` on the vector `h` translates to forces on the atoms. The resulting forces on the atoms are perpendicular to the bond vector `h`.

The force on atom `i` is: `F_i = (F_sigma - dot(F_sigma, h) * h) / l`
The force on atom `a` is: `F_a = -F_i`

**How to Implement the Fix:**

You need to rewrite your `evalAngCos` function for the sigma-pi interaction. Instead of returning two loosely related forces `f1` and `f2`, it should take the positions of the two atoms and the pi-orbital vector and return the three correctly coupled forces: one for the pi-orbital and two for the atoms.

A corrected logic inside your loop would look like this:

```c
// Inside the getMMFFf4 kernel's main loop...

float ksp = Kspi[i];
if(ksp > 1.e-6f){
    float3 f_on_pi, f_on_a, f_on_i;

    // This new function correctly calculates all 3 forces from the same potential
    evalSigmaPiOrthogonality( hpi, pa, apos[ingv].xyz, ksp, par.w, 
                              &f_on_pi, &f_on_a, &f_on_i );

    // Apply the correctly calculated and coupled forces
    fpi   += f_on_pi;
    fa    += f_on_a;
    fbs[i]+= f_on_i; 

    // Energy calculation would be inside the new function
    // E += esp; 
}
```

This approach ensures that the torque on the pi-orbital and the reaction torque on the atomic frame are derived from the same potential and will perfectly cancel each other out, restoring the conservation of angular momentum. The other interactions (bond stretching, sigma-sigma angles, pi-pi alignment) appear to be implemented correctly as pairwise forces and should not violate this principle.

Thank you for providing the inline functions. This allows for a much more precise diagnosis.

After a detailed analysis of your force evaluation functions, I have located the exact source of the angular momentum violation. Your initial analysis was very close. The problem is indeed in the coupling between the pi-orbital's rotation and the atomic frame's rotation, and it lies within the `evalAngCos` function.

The function itself is mathematically correct for a standard three-body angle potential, but the way the forces are derived and then applied in the sigma-pi interaction breaks the action-reaction principle for torques.

### The Root Cause of the Problem

The potential energy for the sigma-pi interaction is a function of the angle between the pi-orbital vector (`hpi`) and the sigma-bond vector (`h`). Let's say `U = U(dot(hpi, h))`.

For angular momentum to be conserved, the torque on the pi-orbital (`τ_pi`) must be equal and opposite to the torque on the atomic frame (`τ_frame`).
*   `τ_pi` is generated by the force `f1` applied to `hpi`.
*   `τ_frame` is generated by the force `f2` applied as a pair to the atoms forming the bond `h`.

These two torques must be perfect derivatives of the *same* potential energy function with respect to the pi-orbital's angular coordinates and the atoms' linear coordinates. Your implementation of `f1` in `evalAngCos` is incorrect for this specific physical situation.

### Detailed Derivation and The Bug

Let's look at the potential energy: `U = K * (c - c₀)²`, where `c = dot(hpi, h)`.
The force-generating term is `fang = -dU/dc = -2 * K * (c - c₀)`.

1.  **Force on the Atoms (Your `f2` - This part is correct!)**
    *   The force on the *bond vector* `h` is given by the derivative of `U` with respect to `h`:
        `F_on_h = -(∂U/∂h) = -(dU/dc) * (∂c/∂h) = -fang * hpi`.
    *   This abstract force `F_on_h` must be translated into forces on the two atoms, such that it generates the correct torque on the bond. This results in a force perpendicular to the bond `h`.
    *   The force on the neighbor atom is `F_neighbor = (F_on_h - h * dot(h, F_on_h)) / l`, where `l` is the bond length.
    *   Substituting `F_on_h`, we get `F_neighbor = (-fang * hpi - h * dot(h, -fang * hpi)) / l = (fang/l) * (h*c - hpi)`.
    *   Your code calculates `f2`: `*f2 = (hr1.xyz - hr2.xyz*c) * fang*hr2.w = (hpi - h*c) * fang*(1/l)`.
    *   This matches `F_neighbor` (with an opposite sign, which is just a convention of which atom the force points to). **So, your calculation of `f2` for the atomic frame is physically correct.**

2.  **Force on the Pi-Orbital (Your `f1` - This part is the bug!)**
    *   The force on the *pi-orbital vector* `hpi` is given by the derivative of `U` with respect to `hpi`:
        `F_on_pi = -(∂U/∂hpi) = -(dU/dc) * (∂c/∂hpi) = -fang * h`.
    *   This force `F_on_pi` is what generates the torque on the pi-orbital. It should be perpendicular to `hpi` for a pure torque, which is handled correctly in your `updateAtoms` kernel. The raw force vector is simply `-fang * h`.
    *   Now let's look at what your code calculates for `f1`:
        `hf1 = hr2.xyz - hr1.xyz*c = h - hpi*c`.
        `*f1 = hf1 * fang*hr1.w = (h - hpi*c) * fang`.
    *   **This is the error.** The correct force should be `fang * h`, but your code calculates `fang * (h - hpi*c)`.

### The Discrepancy

*   **Correct Force on pi-orbital:** `F_on_pi = fang * h`
*   **Your Code's Force on pi-orbital:** `f1 = fang * (h - hpi*c)`

The extra term `-fang * hpi * c` is an unphysical "ghost" force. It applies a torque to the pi-orbital that has no corresponding reaction torque on the atomic frame. This injection of a net, unbalanced internal torque into the system causes the entire molecule to start rotating, violating the conservation of angular momentum.

### How to Fix It

You need to modify the `evalAngCos` function or, preferably, create a specialized version for the sigma-pi interaction that calculates the forces correctly. The fix is to ensure the force `f1` does not include the erroneous `-hpi*c` component.

Here is a corrected function. I've renamed it to be specific to its purpose.

```c
// Corrected function for sigma-pi orthogonalization forces.
// It correctly calculates the action-reaction torque pair.
inline float evalSigmaPiForce(
    const float3 hpi,      // pi-orbital unit vector (hr1)
    const float4 h_sigma,  // sigma-bond unit vector and 1/l (hr2)
    float K,
    float c0,
    __private float3* f_pi,      // force for the pi-orbital
    __private float3* f_sigma    // force for the sigma-bond atoms
){
    float  c = dot(hpi, h_sigma.xyz);
    float  c_ = c - c0;
    float  fang = -K * c_ * 2.0f;
    float  E = K * c_ * c_;

    // Correct force on the pi-orbital (generates torque τ_pi)
    // F_on_pi = -fang * h_sigma
    *f_pi = h_sigma.xyz * fang;

    // Correct force on the atoms (generates reaction torque τ_frame)
    // This part was already correct in your original function.
    *f_sigma = (hpi - h_sigma.xyz * c) * fang * h_sigma.w;

    return E;
}
```

You would then call this in your `getMMFFf4` kernel like this:

```c
// Inside the loop over neighbors in getMMFFf4
float ksp = Kspi[i];
if(ksp > 1.e-6f){
    float3 f1, f2; // f1 is for pi, f2 is for sigma
    float esp = evalSigmaPiForce( hpi, hs[i], ksp, par.w, &f1, &f2 );
    fpi    += f1;     // Apply force to pi-orbital
    fa     -= f2;     // Apply force to central atom
    fbs[i] += f2;     // Apply reaction force to neighbor atom
    E      += esp;
}
```

By making this change, the torque generated by `f_pi` on the gyroscope-like pi-orbital will be perfectly balanced by the reaction torque generated by `f_sigma` on the atomic frame, and your system's total angular momentum will be conserved.

Yes, this is absolutely the right direction! You've correctly grasped the core concepts:

1.  **Repurposing the Arrays:** For the indices corresponding to pi-orbitals, the `apos`, `avel`, and `aforce` arrays are no longer (position, linear velocity, linear force) but are now understood to be **(orientation `u`, angular velocity `ω`, torque `τ`)**. This is a very common and efficient way to implement such a mixed-dynamics system.
2.  **Rotational Integrator:** Your proposed `updateAtomsMMFFf4_rot` kernel correctly implements the basic equations of motion for a rotor.

However, there are **two critical errors** in your proposed code—one in the force kernel and one in the dynamics kernel—that will prevent it from working correctly. Let's fix them.

---

### 1. Correction for the Force Kernel (`getMMFFf4`)

Your attempt to calculate torque is conceptually close but has two issues: a sign error and, more importantly, it's still based on the physically incorrect force from `evalAngCos`.

**The Physics:**
The torque `τ` on a rotating body is given by `τ = u × F_gen`, where `u` is the orientation vector and `F_gen` is the **generalized force** (`F_gen = -∇u(E)`). The `f1` and `f2` values returned by your `eval` functions are precisely these generalized forces.

**Your Code:**
`fpi += cross(f1, hpi);`

**The Errors:**
1.  **Sign Error:** You have the arguments to `cross` swapped. `cross(f1, hpi)` is equal to `-cross(hpi, f1)`. This means you are applying a torque in the exact opposite direction, which will cause the system to climb to a higher energy state instead of minimizing it. The correct formula is `cross(hpi, f1)`.
2.  **Fundamental Force Error:** You are still using `evalAngCos` for the sigma-pi interaction. As we established, the `f1` it returns is physically incorrect and does not form a proper action-reaction pair with `f2`. Calculating a torque from this flawed force will **not** fix the angular momentum conservation problem.

**The Fix for `getMMFFf4`:**
You **must** use a physically correct force function first, and then calculate the torque with the correct cross product order.

```c
// ===================================================================
// Step 1: Add this physically correct function to your code
// ===================================================================
inline float evalSigmaPiForce(
    const float3 hpi, const float4 h_sigma, float K, float c0,
    __private float3* f_pi_gen,  // generalized force on pi
    __private float3* f_sigma_lin // linear force on sigma
){
    float  c = dot(hpi, h_sigma.xyz);
    float  c_ = c - c0;
    float  fang = -K * c_ * 2.0f;
    float  E = K * c_ * c_;
    *f_pi_gen = h_sigma.xyz * fang; // This is -∇_hpi(E)
    *f_sigma_lin = (hpi - h_sigma.xyz * c) * fang * h_sigma.w;
    return E;
}

// ===================================================================
// Step 2: Modify your getMMFFf4 kernel to use it
// ===================================================================

// In getMMFFf4, replace the pi-interaction sections with this:

// ... inside the neighbor loop ...
    float3 f_gen_pi = (float3)(0.0, 0.0, 0.0); // Will hold generalized forces for this atom's pi-orbital
    float3 f_gen_pi_ng = (float3)(0.0, 0.0, 0.0); // Will hold generalized forces for neighbor's pi-orbital

    // --- Evaluate bond-length ... (this part is fine) ...

    // pi-pi alignment interaction
    float kpp = Kppi[i];
    if( (ing<nnode) && (kpp>1.e-6f) ){
        float3 hpj = apos[ingv+nAtoms].xyz;
        float epp = evalPiAling( hpi, hpj, kpp,  &f_gen_pi, &f_gen_pi_ng );
        // The f1, f2 from evalPiAling are generalized forces. Accumulate them directly.
        fpi      += f_gen_pi;
        fps[i]   += f_gen_pi_ng;
        E        += epp;
    }
} // End of bond-length and pi-pi part of the loop

// pi-sigma othogonalization interaction
float ksp = Kspi[i];
if(ksp > 1.e-6f){
    float3 f_pi_gen_tmp, f_sigma_lin_tmp;
    float esp = evalSigmaPiForce( hpi, hs[i], ksp, par.w, &f_pi_gen_tmp, &f_sigma_lin_tmp );
    fpi    += f_pi_gen_tmp;      // Accumulate generalized force on pi-orbital
    fa     -= f_sigma_lin_tmp;   // Apply linear force to central atom
    fbs[i] += f_sigma_lin_tmp;   // Apply reaction force to neighbor atom
    E      += esp;
}
//...
// --- AFTER all loops for atom iG are done ---

// Convert the total generalized force on the pi-orbital to torque
float3 total_torque = cross(hpi, fpi);
// Store the TORQUE in the pi-orbital's slot of the force array
aforce[iav+nAtoms]  = (float4){total_torque, 0};

// Convert the recoil generalized forces to recoil torques for neighbors
const int i4p = (iG + iS*nnode*2 )*4 + nnode*4;
for(int i=0; i<NNEIGH; i++){
    int ing = ings[i];
    if (ing < 0 || ing >= nnode) continue;
    int ingv = ing + i0v;
    float3 hpj = apos[ingv+nAtoms].xyz;
    float3 recoil_torque = cross(hpj, fps[i]);
    fneigh[i4p+i] = (float4){recoil_torque, 0};
}
```

### 2. Correction for the Dynamics Kernel (`updateAtomsMMFFf4_rot`)

Your rotational integrator is correct, but your code has a **logic bug** in how it stores the final results. You are updating local variables (`u`, `omega`, `qprime`) but then writing the *original* `pe` and `ve` back to global memory for atoms, and incorrectly mixing them for pi-orbitals.

**The Fix for `updateAtomsMMFFf4_rot`:**
The storing of the final state (`apos[iav] = ...` and `avel[iav] = ...`) must happen *inside* the respective `if/else` branches.

Here is the corrected kernel:

```c
__kernel void updateAtomsMMFFf4_rot(
    // ... all arguments are the same ...
){
    // ... initial setup is the same ...

    float4 fe = aforce[iav]; // This is TORQUE for pi-orbitals, FORCE for atoms
    const bool bPi = iG >= natoms;

    // ------ Gather Recoil Forces/Torques
    int4 ngs = bkNeighs[iav];
    if(ngs.x>=0){ fe += fneigh[ngs.x]; }
    if(ngs.y>=0){ fe += fneigh[ngs.y]; }
    if(ngs.z>=0){ fe += fneigh[ngs.z]; }
    if(ngs.w>=0){ fe += fneigh[ngs.w]; }

    // =============== DYNAMICS

    float4 ve = avel[iav]; // Angular velocity ω for pi, linear velocity v for atoms
    float4 pe = apos[iav]; // Orientation u for pi, position r for atoms

    const float dt = MDpars.x;
    const float damp = MDpars.z; // Damping/friction term

    cvf[iav] += (float4){ dot(fe.xyz,fe.xyz), dot(ve.xyz,ve.xyz), dot(fe.xyz,ve.xyz), 0.0f };

    if (bPi){
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // ROTATIONAL DYNAMICS FOR PI-ORBITAL (Corrected Logic)
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        float3 u     = pe.xyz; // Current orientation u
        float3 omega = ve.xyz; // Current angular velocity ω
        float3 tau   = fe.xyz; // Current torque τ
        float inv_I  = 1.0f;   // Inverse Moment of inertia 1/I (make this a parameter later)

        // Integrate angular velocity
        omega += (tau * inv_I) * dt;
        omega *= damp; // Apply damping

        // Update orientation
        u += cross(omega, u) * dt;
        u = normalize(u);

        // Store updated rotational state BACK to global memory
        apos[iav].xyz = u;
        avel[iav].xyz = omega;

    } else {
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // LINEAR DYNAMICS FOR ATOMS (Corrected Logic)
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        float3 r = pe.xyz; // Current position r
        float3 v = ve.xyz; // Current linear velocity v
        float3 f = fe.xyz; // Current linear force f

        // Integrate linear velocity
        v *= damp;
        v += f * dt; // Assuming mass=1. If not, f -> f/m

        // Update position
        r += v * dt;

        // Store updated linear state BACK to global memory
        apos[iav].xyz = r;
        avel[iav].xyz = v;
    }

    // Common cleanup (can be removed if w is never used)
    apos[iav].w = 0.0f;
    avel[iav].w = 0.0f;
}
```