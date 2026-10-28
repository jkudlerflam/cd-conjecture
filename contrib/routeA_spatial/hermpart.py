"""Hermitian part of the phase-rotated hopping matrix, from one-particle
matrix elements of the exact second-order bond (no scalar subtracted).
Positive definiteness means every band crosses a small common field.
Usage: python3 hermpart.py.  Writes hermpart.json."""
import json
import numpy as np
from ycell import params, build_cell, code, second_order, one_particle_T
from logical_chain import IR_PARAMS

out = []
for n, Dx in ((2, 4.0), (3, 4.0), (4, 4.0), (4, 10.0), (4, 20.0)):
    P = params(td=0.0, **dict(IR_PARAMS, Dx=Dx, DM=Dx))
    c = build_cell(n, P); E0, G, gap = code(c); H2, k = second_order(c, G, E0)
    T = one_particle_T(H2, n); sc = abs(T).max()
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    ev = np.linalg.eigvalsh((R + R.conj().T) / 2) / sc
    rec = dict(n=n, Dx=Dx, cell_gap=gap, poles=k, T_scale=sc,
               asym=float(np.linalg.norm(R - R.T) / np.linalg.norm(R)), herm_eigs_over_T=list(ev))
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open('hermpart.json', 'w'), indent=1)
