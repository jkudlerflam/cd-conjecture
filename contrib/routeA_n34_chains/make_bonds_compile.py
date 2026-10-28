"""Compute compile-cell bonds K2, K3 for a grid of (tau, phiT) at fixed arm asymmetry xL.
Usage: python3 make_bonds_compile.py m xL "tau1,tau2,..." "phi1,phi2,..."
Writes data/bond_compile_n{m}_tau{tau}_phi{phi}_xL{xL}.npz with K2, K3 (traceless, lambda = 1,
not normalized).  Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, os, time
import numpy as np
import bonds34 as B
m = int(sys.argv[1]); xL = float(sys.argv[2])
taus = [float(x) for x in sys.argv[3].split(',')]
phis = [float(x) for x in sys.argv[4].split(',')]
ch, P = B.compile_prepare(m)
for tau in taus:
    for phi in phis:
        t = time.time()
        B.COMPILE_BREAK.update(xL=xL, phiT=phi)
        K2, K3, info = B.compile_bonds(ch, P, tau)
        fn = os.path.join(B.HERE, 'data', 'bond_compile_n%d_tau%g_phi%g_xL%g.npz' % (m, tau, phi, xL))
        np.savez_compressed(fn, K2=K2, K3=K3, tau=tau, phiT=phi, xL=xL, wD=info['wD'])
        print(fn, '|K2| %.3e |K3| %.3e  %.0fs' % (np.linalg.norm(K2, 2), np.linalg.norm(K3, 2), time.time() - t), flush=True)
