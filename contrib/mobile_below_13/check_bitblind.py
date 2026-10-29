"""Bit blindness of the vacancy dynamics: K = infinity (legal projection), zero
port, zero field, t != 0, w != 0.  Two p=3 cells (two active bits each) and two
vacancies on a 10-site ring.  When the two cells carry different bit strings the
legal closures must be isospectral for every choice; when they are equal the
ring has an extra symmetry, fewer distinct words, and its spectrum must be a
subset of the generic one."""
import os, json, itertools
import numpy as np
import bond12 as b
from closure import closure

HERE = os.path.dirname(os.path.abspath(__file__))
R0, E, V = ('R', 0), ('e', 0), b.V
out = []
for t, w in [(1., 0.), (.3, -.05)]:
    p = b.params(n=3, K=0., musel=.93, field=0., w=w)
    spec = {}
    for bits in itertools.product((0, 1), repeat=4):
        a0, a1, c0, c1 = bits
        seed = [('L', a0), R0, ('R', a1), E, V, ('L', c0), R0, ('R', c1), E, V]
        H, _, _, st = closure(seed, [b.bond(p, lam=0., t=t)] * 10, legal_only=True)
        spec[bits] = np.linalg.eigvalsh(H.toarray())
    gen = [k for k in spec if k[:2] != k[2:]]
    ref = spec[gen[0]]
    d_gen = max(float(np.max(abs(spec[k] - ref))) if len(spec[k]) == len(ref) else np.inf for k in gen)
    d_sub = max(float(max(np.min(abs(ref - x)) for x in spec[k])) for k in spec if k[:2] == k[2:])
    rec = dict(t=t, w=w, generic_dimension=len(ref), symmetric_dimension=len(spec[(0, 0, 0, 0)]),
               max_difference_distinct_cells=d_gen, max_distance_equal_cells_to_generic=d_sub)
    assert d_gen < 1e-11 and d_sub < 1e-11
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open(os.path.join(HERE, 'bitblind_checks.json'), 'w'), indent=1)
