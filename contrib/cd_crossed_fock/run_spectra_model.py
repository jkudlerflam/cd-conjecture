"""Full crossed spectra for a registry model (models.py), resolved by charge and crossed momentum.
Usage: python3 run_spectra_model.py NAME MMAX [betas...]; writes data/spec_NAME_m{m}.npz (MMIN env, default 2)."""
import sys, os, time
import numpy as np
from cflib import Column, Transfer, sector_spectrum
from models import get

HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    name = sys.argv[1]; mmax = int(sys.argv[2])
    betas = [float(b) for b in sys.argv[3:]] or [0.5, 1, 2, 3, 4, 6, 8, 12, 16]
    h, mode, M, A, label = get(name)
    for m in range(int(os.environ.get('MMIN', 2)), mmax + 1):
        col = Column(2, m); Qs = col.charge(mode)
        res = {}; t0 = time.time()
        for b in betas:
            T = Transfer(h, 2, b, m, col)
            oq, on, ou = [], [], []
            for q in np.unique(Qs):
                idx = np.where(Qs == q)[0]
                for n in range(m):
                    u = sector_spectrum(T, idx, n)
                    oq.append(np.full(len(u), q)); on.append(np.full(len(u), n)); ou.append(u)
            u = np.concatenate(ou); assert len(u) == col.N
            res['Q_%g' % b] = np.concatenate(oq).astype(np.int8); res['n_%g' % b] = np.concatenate(on).astype(np.int16); res['u_%g' % b] = u
        res['betas'] = np.array(betas)
        np.savez_compressed(os.path.join(HERE, 'data', 'spec_%s_m%d.npz' % (name, m)), **res)
        print(name, 'm', m, 'done in %.1fs' % (time.time() - t0), flush=True)
