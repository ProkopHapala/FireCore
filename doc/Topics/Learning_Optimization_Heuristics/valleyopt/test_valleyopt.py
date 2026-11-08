import unittest
import numpy as np
from valleyopt import METHODS, RidgeTR, WallProbeWalk, SecantSketch, MomentumWalk, CurvedWalk, minimize
from firecore_adapter import ArrayEnergy

class TestBase(unittest.TestCase):
    def test_smoke_and_memory(self):
        def energy(x): return float(.5*np.dot(x-1.,x-1.))
        for name, method in METHODS.items():
            with self.subTest(method=name):
                s=minimize(energy,method,np.array([2.,0.,-1.]),max_evals=110,step=.2,seed=11)
                self.assertLessEqual(s.e,energy(np.array([2.,0.,-1.])))
                self.assertLessEqual(s.best_e,s.e+1.e-12)
                self.assertEqual(len(s.trace.best_energy),s.trace.nfev)
                self.assertTrue(np.all(np.diff(s.trace.best_energy)<=1.e-12))
                self.assertTrue(np.all(np.isfinite(s.x)))
                self.assertTrue(np.all(np.isfinite(s.best_x)))

    def test_ridge_probe_keeps_anchor(self):
        def energy(x): return (x[0]-2)**2+(x[1]+1)**2
        s=RidgeTR([0.,0.],energy(np.zeros(2)),step=.1,seed=1)
        # discover direction, accept or reject; next calls initialize the model
        trial=s.ask(); s.tell(energy(trial))
        x=s.x.copy(); e=s.e
        trial=s.ask(); self.assertEqual(s._pending[0],'probe')
        s.tell(energy(trial))
        np.testing.assert_allclose(s.x,x)
        self.assertEqual(s.e,e)

    def test_wall_probe_keeps_anchor(self):
        # deterministic manual force initialization to ensure a rejected forward step
        def energy(x): return float(np.dot(x,x))
        s=WallProbeWalk([0.,0.],0.,step=.2,seed=1)
        s.t=np.array([1.,0.])
        trial=s.ask(); s.tell(energy(trial))
        self.assertEqual(s.phase,'plus')
        x=s.x.copy(); e=s.e
        for _ in range(2):
            trial=s.ask(); s.tell(energy(trial))
            np.testing.assert_allclose(s.x,x)
            self.assertEqual(s.e,e)

    def test_exact_quadratic_model(self):
        H = np.array([[2.0,0.4],[0.4,7.0]])
        b = np.array([-.5,1.2])
        def f(x): return float(.5*x@H@x + b@x)
        s=RidgeTR(np.array([.3,-.2]), f(np.array([.3,-.2])), step=.1, seed=17)
        y=s.ask(); s.tell(f(y))  # discover local basis and recenter
        anchor=s.x.copy()
        for _ in range(5):
            y=s.ask()
            self.assertEqual(s._pending[0], 'probe')
            s.tell(f(y))
        np.testing.assert_allclose(s.x,anchor)
        np.testing.assert_allclose(s.g,s.U.T@(H@anchor+b),rtol=1.e-8,atol=1.e-8)
        np.testing.assert_allclose(s.H,s.U.T@H@s.U,rtol=1.e-8,atol=1.e-8)

    def test_seed_isolation_in_benchmark(self):
        from bench import rotated_valley
        _, x0=rotated_valley(n=1000,seed=10006)
        s=RidgeTR(x0,1.2,step=.1,seed=3)
        direction=s.ask()-x0
        # Avoid artificially aligned random probes (a prior bug in bench.py).
        self.assertLess(abs(float(np.dot(direction,x0))/(np.linalg.norm(direction)*np.linalg.norm(x0))),.2)

    def test_array_adapter(self):
        physical=np.array([3.,2.,.5],dtype=np.float64)
        def energy(): return float(np.dot(physical,physical))
        obj=ArrayEnergy(physical,energy,[True,False,True],scale=[.1,1.],positive_indices=[2])
        self.assertEqual(obj.ndim,2)
        obj(np.array([-2.,1.]))
        np.testing.assert_allclose(physical,[2.8,2.,1.5])
        self.assertEqual(obj(np.array([0.,-2.])),float('inf'))

if __name__=='__main__': unittest.main()
