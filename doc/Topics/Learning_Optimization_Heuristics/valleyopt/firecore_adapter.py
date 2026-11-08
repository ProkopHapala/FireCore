"""Thin illustrative bridge to FireCore-style writable NumPy DOF views.

Usage for eFF after setup and eFF.getBuffs():
    from pyBall import eFF as eff
    adapter = ArrayEnergy(eff.pDOFs, eff.eval, free_mask=~fixed_mask)
    solver = minimize(adapter, RidgeTR, np.zeros(adapter.ndim), ...)
    adapter.commit(solver.best_x)

Only the adapter knows which physical variables the engine uses.
"""
import numpy as np

class ArrayEnergy:
    def __init__(self, position_view, eval_function, free_mask=None,
                 scale=1., positive_indices=None, min_positive=1.e-4):
        arr = np.asarray(position_view)
        if arr.dtype != np.float64 or not arr.flags.c_contiguous or not arr.flags.writeable:
            raise TypeError('This adapter expects a writable contiguous float64 NumPy view')
        self.view = arr
        self.shape = arr.shape
        self.base = arr.ravel().copy()
        self.eval_function = eval_function
        if free_mask is None: free_mask = np.ones_like(self.base,dtype=bool)
        mask = np.asarray(free_mask,dtype=bool).ravel()
        if mask.shape!=self.base.shape: raise ValueError('free_mask shape mismatch')
        self.free = np.flatnonzero(mask)
        self.ndim = len(self.free)
        self.scale = np.broadcast_to(np.asarray(scale,dtype=float), (self.ndim,)).copy()
        if np.any(self.scale<=0): raise ValueError('scale must be positive')
        self.positive_indices = np.array([] if positive_indices is None else positive_indices,dtype=int)
        self.min_positive=float(min_positive)

    def coordinates(self, z):
        z=np.asarray(z,dtype=float).ravel()
        if z.shape!=(self.ndim,): raise ValueError('Optimization coordinate shape mismatch')
        x=self.base.copy()
        x[self.free]=self.base[self.free]+self.scale*z
        return x

    def __call__(self,z):
        x=self.coordinates(z)
        if not np.all(np.isfinite(x)) or (len(self.positive_indices) and
                 np.any(x[self.positive_indices]<=self.min_positive)):
            return float('inf')
        self.view.ravel()[:] = x
        e=float(self.eval_function())
        return e if np.isfinite(e) else float('inf')

    def commit(self,z):
        """Set FireCore state to best solution; force/energy buffers updated."""
        self.view.ravel()[:] = self.coordinates(z)
        return float(self.eval_function())
