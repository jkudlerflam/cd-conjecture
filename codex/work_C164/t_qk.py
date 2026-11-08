import numpy as np, time, qband, est, sys
kappa, mu = 0.05, 2.3
L = int(sys.argv[1])
for e in (0.3, 0.5):
    for h in (0.0, 0.03, 0.08, 0.15):
        v = np.array([-1, -1 + mu*e*e/2, -1 + kappa*e*e/2, 0, 0, -e, -h*e*e, 0, 0])
        J, f = est.J_f(v); hh = est.bond_from_JB(J, f)
        t = time.time(); o = qband.analyze(hh, L)
        print(f"e={e} h={h} m0={o['m0']} gap={o['gap']:.5f} E1(k)={np.round(o['E1'],4).tolist()} mins(k/pi,E)={[(round(k/np.pi,3), round(E,5)) for k,E in o['mins'][:3]]} {time.time()-t:.1f}s", flush=True)
