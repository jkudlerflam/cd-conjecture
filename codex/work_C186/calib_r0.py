"""C186 calibration (R3): r = 0 must reproduce the collaborator's binding boundary v_c = g u^2/[2(2-g)] (= u^2/4 at g=2/3)."""
import numpy as np, time
from mixed_fewbody import binding, Model
print("hermiticity check:", end=" ")
m = Model(1.0, np.pi/2, 0.3, mu=-2.1)
for Nmax in (2, 3, 4):
    m.lowest(0.7, Nmax, 14, check=True)
print("ok")
t0 = time.time()
for u in (1.0, 2.0, 0.5):
    vc = u**2/4
    for v in (0.0, 0.5*vc, 0.9*vc, 1.1*vc, 1.5*vc):
        for R in (200, 400, 800):
            binding(u, np.pi/2, v, Nmax=2, R=R, r=0.0)
print("time", time.time()-t0)
