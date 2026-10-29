"""(1) Bond rebuild and zero-sector equality.
 (a) my13.bond against Astra bond13.bond entry by entry at generic points, coloured
     and uncoloured; Hermiticity; commutation with N, M, Q.
 (b) Bond level: h11(C88) + K Pi_(L1,marker) equals the C84 material bond
     (Astra bond11 in work_C84, and the previous referee's independent my11) and the
     occupied block of C85's bond14 divided by J.
 (c) Ring level, exhaustive over all 11^L material words for L = 2..5: the C88 ring
     Hamiltonian and the C85/C84 ring Hamiltonian (with pin) both preserve the F = 0
     subspace and agree on it entry by entry, at finite K and nonzero port.
Usage: python3 check_bond13.py <dir containing work_C84, work_C85, work_C88 copies>"""
import sys, os, json, itertools, importlib.util
import numpy as np
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13

ROOT = sys.argv[1]


def load(name, path):
    sys.path.insert(0, os.path.dirname(path))
    for k in ('parent_shift', 'bond18', 'bond11'):
        sys.modules.pop(k, None)
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); sys.path.pop(0); return m


b13 = load('b13', os.path.join(ROOT, 'work_C88', 'bond13.py'))
b11 = load('b11', os.path.join(ROOT, 'work_C84', 'bond11.py'))
b14 = load('b14', os.path.join(ROOT, 'work_C85', 'bond14.py'))
sys.path.insert(0, os.path.join(HERE, '..', 'referee_C84_C85'))
import my11

out = {'bond_comparison': [], 'parent_comparison': [], 'ring_zero_sector': []}
pts = [dict(p=2, K=37., musel=.93, field=.021, flux=.37, theta=-.21, g=.13, td=.21, D=3.1, lam=.17, J=.7, t=1.3, mu=.4),
       dict(p=3, K=1e3, musel=0., field=0., flux=np.pi / 6, theta=-np.pi / 4, g=.07, td=.003, D=4., lam=.1, J=1., t=1., mu=0.),
       dict(p=14, K=2. ** 20, musel=1.889, field=1e-3, flux=np.pi / 28, theta=-np.pi / 4, g=.01, td=.002, D=1000., lam=.3, J=.2, t=2.5, mu=-.3),
       dict(p=5, K=11., musel=-.4, field=.5, flux=1.1, theta=.8, g=.9, td=.6, D=2., lam=.9, J=2., t=.3, mu=1.1)]
for q in pts:
    n = q['p']
    pa = b13.params(n, K=q['K'], musel=q['musel'], field=q['field'], flux=q['flux'], theta=q['theta'], g=q['g'],
                    td=q['td'], D=q['D'])
    pm = my13.params(n, K=q['K'], musel=q['musel'], field=q['field'], phiA=q['flux'], phiB=0.,
                     psiA=n * q['flux'] + np.pi / 2, psiB=np.pi / 2, theta=q['theta'], g=q['g'], td=q['td'], D=q['D'],
                     Delta=pa['Delta'])
    rec = dict(point={k: v for k, v in q.items()})
    for col in (True, False):
        ha = b13.bond(pa, lam=q['lam'], J=q['J'], t=q['t'], mu=q['mu'], colored=col)
        hm, ls = my13.bond(pm, q['lam'], q['J'], q['t'], q['mu'], coloured=col)
        rec['max_entry_diff_coloured' if col else 'max_entry_diff_uncoloured'] = float(abs(ha - hm).max())
        if col:
            rec['hermiticity'] = float(abs(hm - hm.conj().T).max())
            rec['nonzero_offdiag'] = int(np.count_nonzero(hm - np.diag(np.diag(hm))))
            for name, qq in [('N', [int(not my13.is_vac(a)) for a in ls]),
                             ('M', [int(a in my13.MARKERS) for a in ls]), ('Q', [my13.charge(a) for a in ls])]:
                qb = np.add.outer(qq, qq).ravel()
                rec['commutator_' + name] = float(abs((qb[:, None] - qb[None, :]) * hm).max())
    rec['allowed_table_agrees'] = all(my13.allowed(a, b) == b13.allowed(b13.LET[i], b13.LET[j])
                                      for i, a in enumerate(my13.LET) for j, b in enumerate(my13.LET))
    # parents
    pin = np.diag([float(a == 'L1' and b in my13.MARKERS) for a in my13.MAT for b in my13.MAT])
    h88 = my13.h11(pm, q['lam'])
    p84 = b11.bond(pa, lam=q['lam'])
    h14 = b14.bond(pa, lam=q['lam'], J=q['J'], t=q['t'], mu=0., colored=True)
    occ = [14 * i + j for i in range(11) for j in range(11)]
    h85 = h14[np.ix_(occ, occ)] / q['J']
    pr = my11.default_params(n, D=q['D'], Delta=pa['Delta'], g=q['g'], td=q['td'], phiA=q['flux'], phiB=0.,
                             psiA=n * q['flux'] + np.pi / 2, psiB=np.pi / 2, K=q['K'], musel=q['musel'],
                             field=q['field'])
    ok_theta = abs(q['theta'] + np.pi / 4) < 1e-15
    rec['C84_astra_minus_pin'] = float(abs(h88 + q['K'] * pin - p84).max())
    rec['C85_occupied_block'] = float(abs(h88 + q['K'] * pin - h85).max())
    if ok_theta:
        rec['C84_referee_my11'] = float(abs(h88 + q['K'] * pin - my11.bond(pr, q['lam'])).max())
    out['bond_comparison'].append(rec)
    print(json.dumps(rec), flush=True)


# (c) ring level, exhaustive
def ring_matrix(h, L, d):
    """Sparse ring Hamiltonian sum_i h_{i,i+1} on (C^d)^L, periodic, from a d^2 bond."""
    N = d ** L
    idx = np.arange(N)
    digs = np.array([(idx // d ** (L - 1 - k)) % d for k in range(L)])  # digs[k] = letter at site k
    hs = sp.coo_matrix(h)
    rows, cols, vals = [], [], []
    for i in range(L):
        j = (i + 1) % L
        pair = digs[i] * d + digs[j]
        # for every nonzero h[r, c], words with pair == c map to words with pair r
        for r, c, v in zip(hs.row, hs.col, hs.data):
            ws = idx[pair == c]
            if len(ws) == 0:
                continue
            ar, br = divmod(r, d)
            new = ws + (ar - digs[i][ws]) * d ** (L - 1 - i) + (br - digs[j][ws]) * d ** (L - 1 - j)
            rows.append(new); cols.append(ws); vals.append(np.full(len(ws), v))
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)).tocsr(), digs


q = pts[0]
pa = b13.params(3, K=q['K'], musel=q['musel'], field=q['field'], flux=q['flux'], theta=q['theta'], g=q['g'], td=q['td'], D=q['D'])
pm = my13.params(3, K=q['K'], musel=q['musel'], field=q['field'], phiA=q['flux'], phiB=0., psiA=3 * q['flux'] + np.pi / 2,
                 psiB=np.pi / 2, theta=q['theta'], g=q['g'], td=q['td'], D=q['D'], Delta=pa['Delta'])
h88 = my13.h11(pm, q['lam'])
h14 = b14.bond(pa, lam=q['lam'], J=1., t=1., mu=0., colored=True)
occ = [14 * i + j for i in range(11) for j in range(11)]
h85 = h14[np.ix_(occ, occ)]
for L in (2, 3, 4, 5):
    H88, digs = ring_matrix(h88, L, 11)
    H85, _ = ring_matrix(h85, L, 11)
    i1 = my13.IXM['L1']
    mk = [my13.IXM[m] for m in my13.MARKERS]
    Fv = sum(((digs[i] == i1) & np.isin(digs[(i + 1) % L], mk)).astype(int) for i in range(L))
    Hc = H88.tocoo()
    cons88 = bool(np.all(Fv[Hc.row] == Fv[Hc.col]))
    Hc5 = H85.tocoo()
    cons85 = bool(np.all(Fv[Hc5.row] == Fv[Hc5.col]))
    z0 = np.flatnonzero(Fv == 0)
    D = (H88[z0][:, z0] - H85[z0][:, z0])
    diff = float(abs(D).max()) if D.nnz else 0.
    full = H85 - H88 - sp.diags(q['K'] * Fv.astype(float))
    rec = dict(L=L, words=11 ** L, F0_words=len(z0), C88_conserves_F=cons88, C85_conserves_F=cons85,
               max_diff_on_F0=diff, identity_H85_eq_H88_plus_KF=float(abs(full).max()) if full.nnz else 0.)
    out['ring_zero_sector'].append(rec)
    print(json.dumps(rec), flush=True)
json.dump(out, open(os.path.join(HERE, 'check_bond13.json'), 'w'), indent=1)
