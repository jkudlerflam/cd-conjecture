"""Symmetry-resolved entanglement of saved iMPS: S = sum_q p_q S_q + H({p_q}) on each bond (parity q = 0,1).
Writes cache out3/<name>.sq.npz with p (bond,q), Sq (bond,q), S, Sconf = sum_q p_q S_q.  Cheap; run on a node.
usage: python sconf.py [substr ...]"""
import sys, glob, os
import numpy as np
import h5py
from tenpy.tools import hdf5_io

for fn in sorted(glob.glob("out3/*.h5")):
    if fn.endswith(".part.h5") or (sys.argv[1:] and not all(a in fn for a in sys.argv[1:])):
        continue
    out = fn[:-3] + ".sq.npz"
    if os.path.exists(out) and os.path.getmtime(out) > os.path.getmtime(fn):
        continue
    with h5py.File(fn, "r") as f:
        psi = hdf5_io.load_from_hdf5(f)["psi"]
    p = np.zeros((2, 2)); Sq = np.zeros((2, 2)); S = np.zeros(2)
    for b in range(2):
        for q, ent in psi.entanglement_spectrum(by_charge=True)[b]:
            w = np.exp(-np.asarray(ent))      # Schmidt weights s^2 (TeNPy returns -log(s^2))
            qq = int(np.asarray(q).ravel()[0]) % 2
            p[b, qq] += w.sum()
            Sq[b, qq] += -np.sum(w * np.log(np.maximum(w, 1e-300)))   # accumulate -sum w log w, normalize below
        for qq in range(2):
            if p[b, qq] > 0:
                Sq[b, qq] = Sq[b, qq] / p[b, qq] + np.log(p[b, qq])    # entropy of normalized block
        S[b] = -np.sum([pp * np.log(pp) for pp in p[b] if pp > 0]) + np.sum(p[b] * Sq[b])
    Sconf = np.sum(p * Sq, axis=1)
    np.savez(out, p=p, Sq=Sq, S=S, Sconf=Sconf)
    print(fn, "p0=%.4f,%.4f S=%.5f Sconf=%.5f" % (p[0, 0], p[1, 0], S.mean(), Sconf.mean()), flush=True)
