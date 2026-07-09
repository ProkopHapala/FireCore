/// @file FFfit_lib.cpp @brief C-wrapper for FFfit class (force-field parameter fitting to Hessians)
/// @ingroup Classical_Molecular_Mechanics

#include "FFfit.h"
#include <cstring>

// Handle-based multi-instance: each Python FFfit() gets its own C++ instance
extern "C" {

// === Lifecycle ===

void* fffit_create() {
    return new FFfit();
}

void fffit_destroy(void* handle) {
    delete (FFfit*)handle;
}

// === Geometry / topology ===

void fffit_set_geometry(void* h, int natoms, double* positions) {
    ((FFfit*)h)->set_geometry(natoms, (Vec3d*)positions);
}

void fffit_add_bond(void* h, int i, int j, double l0) {
    ((FFfit*)h)->add_bond(i, j, l0);
}

void fffit_add_angle(void* h, int i, int j, int k, double theta0) {
    ((FFfit*)h)->add_angle(i, j, k, theta0);
}

void fffit_set_symbols(void* h, int nsym, char** syms) {
    std::vector<std::string> symbols(nsym);
    for (int i = 0; i < nsym; i++) symbols[i] = syms[i];
    ((FFfit*)h)->set_symbols(symbols);
}

void fffit_auto_assign_types(void* h) {
    ((FFfit*)h)->auto_assign_types();
}

void fffit_set_n_free(void* h, int n_free) {
    ((FFfit*)h)->param_map.n_free_params = n_free;
}

void fffit_set_bond_param(void* h, int ibond, int param_idx) {
    ((FFfit*)h)->set_bond_param(ibond, param_idx);
}

void fffit_set_angle_param(void* h, int iangle, int param_idx) {
    ((FFfit*)h)->set_angle_param(iangle, param_idx);
}

int fffit_get_n_free(void* h) {
    return ((FFfit*)h)->param_map.n_free_params;
}

int fffit_get_natoms(void* h) {
    return ((FFfit*)h)->natoms;
}

int fffit_get_nbonds(void* h) {
    return (int)((FFfit*)h)->bonds.size();
}

int fffit_get_nangles(void* h) {
    return (int)((FFfit*)h)->angles.size();
}

// === Fitting ===

/// Fit via linear least-squares. H_ref: (3N x 3N) flattened row-major.
/// Returns fitted params in k_out (must be pre-allocated with n_free doubles).
void fffit_fit_lsq(void* h, double* H_ref, double* H_0, double* weight, double* k_out) {
    auto k = ((FFfit*)h)->fit_hessian(H_ref, H_0, weight);
    memcpy(k_out, k.data(), k.size() * sizeof(double));
}

void fffit_fit_gd(void* h, double* H_ref, double* H_0, double* weight,
                  double lr, double momentum, int max_iter, double tol,
                  int verbose, double* k_out) {
    auto k = ((FFfit*)h)->fit_gradient_descent(H_ref, H_0, weight, lr, momentum, max_iter, tol, verbose);
    memcpy(k_out, k.data(), k.size() * sizeof(double));
}

void fffit_compute_model_hessian(void* h, double* H_out) {
    auto H = ((FFfit*)h)->compute_model_hessian();
    memcpy(H_out, H.data(), H.size() * sizeof(double));
}

void fffit_get_params(void* h, double* k_out) {
    auto& k = ((FFfit*)h)->k_free;
    memcpy(k_out, k.data(), k.size() * sizeof(double));
}

void fffit_set_params(void* h, double* k_in, int np) {
    ((FFfit*)h)->k_free.assign(k_in, k_in + np);
}

void fffit_get_bond_param_idx(void* h, int* idx_out) {
    auto& idx = ((FFfit*)h)->param_map.bond_param_idx;
    memcpy(idx_out, idx.data(), idx.size() * sizeof(int));
}

void fffit_get_angle_param_idx(void* h, int* idx_out) {
    auto& idx = ((FFfit*)h)->param_map.angle_param_idx;
    memcpy(idx_out, idx.data(), idx.size() * sizeof(int));
}

void fffit_print_params(void* h) {
    ((FFfit*)h)->print_params();
}

// === Multi-system accumulation ===

/// Accumulate one system's contribution to normal equations G (n_free*n_free), y (n_free).
/// G and y must be pre-allocated and zero-initialized by caller.
/// H_ref: (3N x 3N) flattened, weight: (3N*3N) flattened or nullptr.
void fffit_accumulate_normal_equations(void* h, double* G, double* y, double* H_ref, double* H_0, double* weight) {
    FFfit* ff = (FFfit*)h;
    int np = ff->param_map.n_free_params;
    int n3 = ff->natoms * 3;
    auto A = ff->build_sensitivity_matrices();
    for (int p = 0; p < np; p++) {
        for (int q = p; q < np; q++) {
            double dot = 0.0;
            for (int i = 0; i < n3 * n3; i++) {
                double w = weight ? weight[i] * weight[i] : 1.0;
                dot += w * A[p][i] * A[q][i];
            }
            G[p * np + q] += dot;
            if (q != p) G[q * np + p] += dot;
        }
        for (int i = 0; i < n3 * n3; i++) {
            double w = weight ? weight[i] * weight[i] : 1.0;
            double h0 = H_0 ? H_0[i] : 0.0;
            y[p] += w * A[p][i] * (H_ref[i] - h0);
        }
    }
}

/// Compute gradient and loss for current system with given k (n_free).
/// Accumulates into grad_out (n_free) and loss_out (scalar).
/// Caller must pre-initialize grad_out and loss_out to zero.
void fffit_compute_gradient_loss(void* h, double* H_ref, double* H_0, double* weight,
                                  double* k, double* grad_out, double* loss_out) {
    FFfit* ff = (FFfit*)h;
    double loss = 0.0;
    auto grad = ff->compute_gradient(H_ref, H_0, weight,
                                     std::vector<double>(k, k + ff->param_map.n_free_params), &loss);
    int np = ff->param_map.n_free_params;
    for (int f = 0; f < np; f++) grad_out[f] += grad[f];
    *loss_out += loss;
}

/// Solve normal equations G k = y (G is n_free*n_free, y is n_free).
/// k_out must be pre-allocated (n_free).
void fffit_solve_normal_equations(double* G, double* y, double* k_out, int np) {
    // Gaussian elimination with partial pivoting
    for (int i = 0; i < np; i++) {
        int p = i;
        double maxval = fabs(G[i * np + i]);
        for (int r = i + 1; r < np; r++) {
            double val = fabs(G[r * np + i]);
            if (val > maxval) { maxval = val; p = r; }
        }
        if (maxval < 1e-15) continue;
        if (p != i) {
            for (int c = 0; c < np; c++) std::swap(G[i * np + c], G[p * np + c]);
            std::swap(y[i], y[p]);
        }
        for (int r = i + 1; r < np; r++) {
            double factor = G[r * np + i] / G[i * np + i];
            for (int c = i; c < np; c++) G[r * np + c] -= factor * G[i * np + c];
            y[r] -= factor * y[i];
        }
    }
    for (int i = np - 1; i >= 0; i--) {
        k_out[i] = y[i];
        for (int j = i + 1; j < np; j++) k_out[i] -= G[i * np + j] * k_out[j];
        k_out[i] /= G[i * np + i];
    }
}

} // extern "C"
