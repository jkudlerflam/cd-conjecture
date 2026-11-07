"""Full crossed spectra of the checkerboard transfer, resolved by charge Q and crossed momentum n.
Usage: python3 run_spectra.py DELTA MMAX [betas...]
Writes data/spec_D{Delta}_m{m}.npz with, for each beta, arrays Q, n, u (eigenvalues of U_sym).
"""
import sys, os, time
import numpy as np
from cflib import Column, Transfer, bond_xxz_rot, sector_spectrum

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)


def spectrum(Delta, m, beta, col):
    T = Transfer(bond_xxz_rot(Delta), 2, beta, m, col)
    Qs = col.charge()
    out_Q, out_n, out_u = [], [], []
    for q in np.unique(Qs):
        idx = np.where(Qs == q)[0]
        for n in range(m):
            u = sector_spectrum(T, idx, n)
            out_Q.append(np.full(len(u), q)); out_n.append(np.full(len(u), n)); out_u.append(u)
    return np.concatenate(out_Q), np.concatenate(out_n), np.concatenate(out_u)


if __name__ == '__main__':
    Delta = float(sys.argv[1]); mmax = int(sys.argv[2])
    betas = [float(b) for b in sys.argv[3:]] or [0.5, 1, 2, 3, 4, 6, 8, 12, 16]
    mmin = int(os.environ.get('MMIN', 2))
    for m in range(mmin, mmax + 1):
        col = Column(2, m)
        res = {}
        t0 = time.time()
        for b in betas:
            Q, n, u = spectrum(Delta, m, b, col)
            assert len(u) == col.N
            res['Q_%g' % b] = Q.astype(np.int8); res['n_%g' % b] = n.astype(np.int16); res['u_%g' % b] = u
        res['betas'] = np.array(betas)
        np.savez_compressed(os.path.join(HERE, 'data', 'spec_D%g_m%d.npz' % (Delta, m)), **res)
        print('Delta', Delta, 'm', m, 'done in %.1fs' % (time.time() - t0), flush=True)
