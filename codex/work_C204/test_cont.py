import numpy as np, time
from c204lib import *
for eta in [0.05, 0.01, 1e-3]:
    b = PI/(PI-eta)
    for hf in [5,10,20]:
        t=time.time()
        c = Continuum(eta, 14.0, hfac=hf)
        v = np.array([0.0, 1.0, 3.0, 6.0])
        r0 = b/(2*PI*np.cosh(b*v))
        s = c.s(v)
        N = np.arctan(np.exp(b*v))/PI
        print(eta, hf, c.N, '%.2e'%np.max(np.abs(s/r0-1)), '%.2e'%np.max(np.abs(c.N0(v)-N)), 'mass-1/2 %.2e'%(c.mass()-0.5), '%.1fs'%(time.time()-t))
# finite Q: convergence in h
for eta in [1e-3]:
    for Q in [0.3, 2.0]:
        vals=[]
        for hf in [5,10,20]:
            c=Continuum(eta,Q,hfac=hf)
            vals.append((c.mass(), c.N0(np.array([Q*0.7]))[0], c.s(np.array([Q]))[0]))
        print(Q, np.array(vals))
