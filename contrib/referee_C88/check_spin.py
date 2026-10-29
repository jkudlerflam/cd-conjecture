"""(5) S_pin: pin defects on small mobile rings, and the variational replacement bound.
For a ring of two material cells of p data sites each plus v vacancies, build the full
legal (K = infinity) sector at fixed cell lengths for each pin number F = 0, 1, 2 (union
of closures from every legal squeezed word with that F, so no charge or bit restriction),
and the ground energies E_F. For the F = 1 ground state psi_pin, measure the contact
probabilities <c_j> on the p bonds of the pinned cell (data_j, data_j+1 and last data,
marker) and evaluate the referee's replacement bound

    E_0 <= E_1 - J lambda_max(P(<c>)) + (port correction),

where P(<c>) is the path matrix with off-diagonal weights <c_j>: replacing the frozen
pinned word by the clock state phi on the all-zero-bit vacuum code L0^j R0^(p-j) E keeps
every vacancy amplitude and changes the energy by -J sum_j 2 phi_j phi_j+1 <c_j>.
Also: frozen rings (v = 0) and a finite-K check."""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13
from my13 import MAT, MARKERS, allowed, is_vac

DATA = [l for l in MAT if l not in MARKERS]


def squeezed_words(p, ncell=2):
    cells = [w + (m,) for w in itertools.product(DATA, repeat=p) for m in MARKERS
             if all(allowed(a, b) for a, b in zip(w + (m,), (w + (m,))[1:]))]
    out = {}
    for cs in itertools.product(cells, repeat=ncell):
        s = sum(cs, ())
        if all(allowed(s[i], s[(i + 1) % len(s)]) for i in range(len(s))):
            out.setdefault(my13.F(s), []).append(s)
    return out


def multi_closure(seeds, h, ls, project=True, coloured=True):
    """Union of ring closures; returns states, H (penalty stripped), pen."""
    dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    T = my13.transitions(h); diag = np.real(np.diag(h))
    states = []; idx = {}
    for s in seeds:
        if s not in idx:
            idx[s] = len(states); states.append(s)
    L = len(states[0])
    r, c, v = [], [], []
    cur = 0
    while cur < len(states):
        w = states[cur]; dg = 0.
        for i in range(L):
            a, b = w[i], w[(i + 1) % L]
            col = ix[a] * dd + ix[b]; dg += diag[col]
            for row, val in T[col]:
                zz = list(w); zz[i] = ls[row // dd]; zz[(i + 1) % L] = ls[row % dd]; zz = tuple(zz)
                if project and my13.nbad(zz, coloured):
                    continue
                if zz not in idx:
                    idx[zz] = len(states); states.append(zz)
                r.append(idx[zz]); c.append(cur); v.append(val)
        if dg:
            r.append(cur); c.append(cur); v.append(dg)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsc()
    pen = np.array([my13.nbad(w, coloured) for w in states], float)
    return states, H, pen


def ground(H, k=1):
    if H.shape[0] < 3000:
        e, U = np.linalg.eigh(H.toarray()); return e[:k], U[:, :k]
    e, U = eigsh(H, k=k, which='SA', tol=1e-11, ncv=40)
    o = np.argsort(e); return e[o], U[:, o]


def pinned_contacts(w, p):
    """Contacts on the p bonds of the (single) pinned cell of a legal F = 1 word."""
    L = len(w)
    mpos = [i for i in range(L) if not is_vac(w[i])]
    n = len(mpos)
    letters = [w[i] for i in mpos]
    k = [j for j in range(n) if letters[j] == 'L1' and letters[(j + 1) % n] in MARKERS]
    assert len(k) == 1
    last = k[0]
    cell = [(last - p + 1 + q) % n for q in range(p + 1)]  # p data letters then the marker
    return np.array([int((mpos[cell[q + 1]] - mpos[cell[q]]) % L == 1) for q in range(p)])


def path_lmax(cw):
    n = len(cw) + 1
    P = np.zeros((n, n))
    for j, x in enumerate(cw):
        P[j, j + 1] = P[j + 1, j] = x
    return float(np.linalg.eigvalsh(P)[-1])


def components(states, ls, hs):
    """Matrices of several bonds on a fixed (invariant) state list."""
    dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    idx = {w: k for k, w in enumerate(states)}
    L = len(states[0])
    Ts = [my13.transitions(h) for h in hs]; Ds = [np.real(np.diag(h)) for h in hs]
    acc = [([], [], []) for _ in hs]
    for cur, w in enumerate(states):
        cols = [ix[w[i]] * dd + ix[w[(i + 1) % L]] for i in range(L)]
        for m, (T, D) in enumerate(zip(Ts, Ds)):
            r, c, v = acc[m]
            dg = sum(D[col] for col in cols)
            if dg:
                r.append(cur); c.append(cur); v.append(dg)
            for i, col in enumerate(cols):
                for row, val in T[col]:
                    zz = list(w); zz[i] = ls[row // dd]; zz[(i + 1) % L] = ls[row % dd]; zz = tuple(zz)
                    k = idx.get(zz)
                    if k is None:
                        continue  # only projected-out illegal words are missing
                    r.append(k); c.append(cur); v.append(val)
    N = len(states)
    return [sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsc() for r, c, v in acc]


CACHE = {}


def sector(p, v, ncell, project):
    key = (p, v, ncell, project)
    if key not in CACHE:
        pm = my13.params(p, K=0., musel=.5)
        hg, ls = my13.bond(pm, .1, 1., 1., 0.)
        hJ, _ = my13.bond(pm, 0., 1., 0., 0.)
        ht, _ = my13.bond(pm, 0., 0., 1., 0.)
        hl = my13.bond(pm, 1., 1., 0., 0.)[0] - hJ
        sq = squeezed_words(p, ncell)
        data = {}
        for Fv in sorted(sq):
            seeds = [my13.colour(s + ('V',) * v) if v else s for s in sq[Fv]]
            st, _, pen = multi_closure(seeds, hg, ls, project=project)
            HJ, Ht, Hl = components(st, ls, [hJ, ht, hl])
            cm = np.array([pinned_contacts(w, p) for w in st]) if (Fv == 1 and project) else None
            data[Fv] = dict(states=st, HJ=HJ, Ht=Ht, Hl=Hl, pen=pen, cmat=cm)
        CACHE.clear(); CACHE[key] = data
    return CACHE[key]


def run(p, v, J, t, lam, K=None, ncell=2):
    data = sector(p, v, ncell, K is None)
    res = dict(p=p, cells=ncell, v=v, J=J, t=t, lam=lam, K='infinity' if K is None else K)
    E = {}
    for Fv, dct in data.items():
        H = J * dct['HJ'] + t * dct['Ht'] + J * lam * dct['Hl']
        if K is not None:
            H = H + sp.diags(J * K * dct['pen'])
        e, U = ground(H)
        E[Fv] = float(e[0]); res[f'dim_F{Fv}'] = len(dct['states'])
        if dct['cmat'] is not None:
            psi = np.abs(U[:, 0]) ** 2
            cav = psi @ dct['cmat']
            res['pinned_cell_contacts'] = cav.tolist()
            res['bound_uniform_phi'] = float(J * 2 / (p + 1) * cav.sum())
            res['bound_optimal_phi'] = float(J * path_lmax(cav))
    res['E'] = E
    res['pin_cost_1'] = E[1] - E[0]
    res['pin_cost_2_minus_2x1'] = (E[2] - E[0]) - 2 * (E[1] - E[0]) if 2 in E else None
    res['frozen_clock_energy_2cos'] = 2 * np.cos(np.pi / (p + 2))
    res['cost_over_J_2cos'] = res['pin_cost_1'] / (J * res['frozen_clock_energy_2cos'])
    if 'bound_optimal_phi' in res:
        res['bound_satisfied'] = bool(res['pin_cost_1'] >= res['bound_optimal_phi'] - 3 * np.sqrt(2) * lam - 1e-10)
        res['bound_satisfied_lam0_form'] = bool(res['pin_cost_1'] >= res['bound_optimal_phi'] - 1e-10)
    print(json.dumps(res), flush=True)
    return res


if __name__ == '__main__':
    out = []
    for (p, vs) in [(2, (0, 1, 2, 3)), (3, (0, 1))]:
        for v in vs:
            ts = (1.,) if v == 0 else ((0.1, 0.3, 1., 3., 10.) if (p, v) in ((2, 1), (2, 2)) else (0.3, 1., 3.))
            for t in ts:
                for lam in ((0., 0.1) if t == 1. or (p, v) in ((2, 1), (2, 2)) else (0.,)):
                    out.append(run(p, v, 1., t, lam))
    # finite K, small
    for K in (30., 100., 1000.):
        out.append(run(2, 1, 1., 1., 0.1, K))
    json.dump(out, open(os.path.join(HERE, 'check_spin.json'), 'w'), indent=1)
