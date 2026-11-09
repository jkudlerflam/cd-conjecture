"""C182: balanced six-particle cluster (aaabbb) at phi = pi/2, K = 3(k_a + k_b) = pi, box R (Dirichlet upper bound).
Usage: python3 ladder6.py wa R"""
import numpy as np, sys, json, time
import nbody
from pairband import valleys
wa = float(sys.argv[1]); R = int(sys.argv[2])
w = 1j*wa; ea = valleys(w)[0][1]
t0 = time.time()
e = nbody.lowest(6, np.pi, w, R)[0][0]
print(json.dumps(dict(wa=wa, R=R, N=6, K=np.pi, E=float(e), per_particle=float(e/6-ea), secs=time.time()-t0)))
