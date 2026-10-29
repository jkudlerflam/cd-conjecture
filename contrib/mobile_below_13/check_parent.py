"""Cross-check bond12 against C88's bond13 (read-only import) and basic algebra."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'codex', 'work_C88'))
import bond13 as c88
sys.path.pop(0)
import bond12 as b

out = {}
rng = np.random.default_rng(1)
errs = []
for trial in range(8):
    n = int(rng.integers(2, 15))
    kw = dict(K=float(rng.uniform(5, 50)), musel=float(rng.uniform(0, 2)), field=float(rng.uniform(0, .1)),
              flux=float(rng.uniform(0, 1)), theta=float(rng.uniform(-1, 1)), g=float(rng.uniform(0, .3)),
              td=float(rng.uniform(0, .1)), D=float(rng.uniform(1, 10)))
    lam = float(rng.uniform(0, .5))
    J = float(rng.uniform(.5, 2))
    p12 = b.params(n, **kw); p13 = c88.params(n, **kw)
    h12 = b.bond(p12, lam=lam, J=J, t=.3)
    h13 = c88.bond(p13, lam=lam, J=J, t=.3)
    # material-material block: identical (C88 embeds on-site as contact on that block)
    occ12 = [12 * i + j for i in range(11) for j in range(11)]
    occ13 = [13 * i + j for i in range(11) for j in range(11)]
    errs.append(float(np.max(abs(h12[np.ix_(occ12, occ12)] - h13[np.ix_(occ13, occ13)]))))
    assert np.max(abs(h12 - h12.conj().T)) == 0
out['material_block_vs_C88_max_error'] = max(errs)
# conserved additive quantities: N (occupied), M (markers), Q (charge)
h = b.bond(b.params(5, K=37., musel=.93, field=.021, flux=.37, theta=-.21, g=.11, td=.05), lam=.17, t=.4)
for name, q in [('N', [int(a[0] != 'V') for a in b.LET]), ('M', [int(b.marker(a)) for a in b.LET]),
                ('Q', [b.charge(a) for a in b.LET]), ('excited_minus_markers_ok', None)]:
    if q is None:
        continue
    qb = np.add.outer(q, q).ravel()
    out['commutator_' + name] = float(np.max(abs((qb[:, None] - qb[None, :]) * h)))
    assert out['commutator_' + name] == 0
out['nonzero_entries'] = int(np.count_nonzero(h))
out['offdiagonal_nonzero_entries'] = int(np.count_nonzero(h - np.diag(np.diag(h))))
# legal vacancy pairs
out['letters_allowed_left_of_V'] = [str(a) for a in b.MAT if b.allowed(a, b.V)]
out['letters_allowed_right_of_V'] = [str(a) for a in b.MAT if b.allowed(b.V, a)]
# squeezing lemma: a V b legal (both pairs) implies a b legal, for all material a,b
viol = [(a, c) for a in b.MAT for c in b.MAT if b.allowed(a, b.V) and b.allowed(b.V, c) and not b.allowed(a, c)]
out['squeeze_violations'] = len(viol)
assert not viol
# material pairs (C88 pinless relation) equal C88's
assert all(b.allowed(a, c) == c88.allowed(a, c) for a in b.MAT for c in b.MAT)
out['material_grammar_equals_C88'] = True
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, 'parent_checks.json'), 'w'), indent=1)
