import numpy as np

def compute_pareto_front(scores, z_spans):
    """
    Find the indices of Pareto optimal points.
    We want to MINIMIZE both scores (clash) and z_spans.
    """
    pts = np.vstack((scores, z_spans)).T
    is_pareto = np.ones(pts.shape[0], dtype=bool)
    for i, p in enumerate(pts):
        if is_pareto[i]:
            # p is Pareto if no other point is strictly better in all dims,
            # or better in one and equal in rest.
            # a dominates b if: a_x <= b_x and a_y <= b_y and (a_x < b_x or a_y < b_y)
            # pts <= p  is a bool array of shape (N, 2)
            dominating = np.all(pts <= p, axis=1) & np.any(pts < p, axis=1)
            is_pareto[i] = not np.any(dominating)
    return is_pareto

# test
s = np.array([1, 2, 3, 2, 1.5])
z = np.array([3, 2, 1, 3, 1.5])
p = compute_pareto_front(s, z)
print(p)
