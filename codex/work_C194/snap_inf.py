"""C194: correlation lengths from an iDMRG checkpoint taken while a chi stage is still iterating (snapshot, not
converged by run_inf's rule). Usage: python3 -E -s snap_inf.py <ckpt.pkl> <cfg.json>"""
import sys, json, pickle, shutil, tempfile
import numpy as np
from run_inf import make_model, spec_rec
c = json.load(open(sys.argv[2]))
tmp = tempfile.mktemp(suffix='.pkl'); shutil.copy(sys.argv[1], tmp)
dm = pickle.load(open(tmp, 'rb')); M = make_model(c); dm.M = M; dm.W.M = M
n2 = 2 * M.q
ev0 = dm.transfer_eigs(0, k=10); T0 = spec_rec(ev0, ev0, n2)
T1 = spec_rec(dm.transfer_eigs(1, string=True, k=6), ev0, n2); T2 = spec_rec(dm.transfer_eigs(2, string=False, k=6), ev0, n2)
print(f"t={dm.t} chi_mid={dm.hist[-1]['chi_mid'] if dm.hist else '?'} xi0={T0[1][2]:.1f} xi1={T1[0][2]:.1f} xi2={T2[0][2]:.1f} ratio={T1[0][2]/T2[0][2]:.3f}")
