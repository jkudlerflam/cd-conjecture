"""Does the cell pattern survive vacancies?  Exact diagonalisation of h10 on a ring: one cell (walker plus n data
sites, bits 0, port and conversion off) with N_v vacancies, against the headless ring R_0^(n+1) with the same N_v.
eps_cell = E0(one cell) - E0(headless) is the energy of putting a head into the material and must be negative.
Both on-site conventions are run: 'contact' (C88's, the definition) and 'genuine' (the first version, which fails)."""
import numpy as np, itertools, json, sys, scipy.sparse as sp, scipy.sparse.linalg as spla
import bond9, bond10
from bond10 import IDX, d
def nu_sel(b):
    a = lambda x: 2*(1-np.cos(np.pi/(x+1))); dl = lambda x: (x+1)*a(x) - x*a(x+1)
    return 2 - (dl(b) + dl(b-1))/2
def ground(material_words, Nv, h):
    """material_words: list of cyclic material words spanning the sector; returns E0 on the ring with Nv vacancies"""
    M = len(material_words[0]); L = M + Nv; V = IDX['V']
    states = []
    for mw in material_words:
        for vac in itertools.combinations(range(L), Nv):
            w = [V]*L; it = iter(mw)
            for i in range(L):
                if i not in vac: w[i] = next(it)
            states.append(tuple(w))
    # close under cyclic rotations of the material word by including all rotations of each seed
    index = {s: i for i, s in enumerate(states)}
    nz = {}
    for r, c in zip(*np.nonzero(np.abs(h) > 1e-14)): nz.setdefault(int(c), []).append((int(r), h[r, c]))
    rows, cols, vals = [], [], []; leak = 0
    for i, s in enumerate(states):
        for site in range(L):
            j = (site+1) % L
            for r, v in nz.get(s[site]*d + s[j], ()):
                t2 = list(s); t2[site] = r//d; t2[j] = r % d
                k = index.get(tuple(t2))
                if k is None: leak += 1
                else: rows.append(k); cols.append(i); vals.append(v)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(states),)*2)
    assert leak == 0, leak
    if H.shape[0] <= 600: return float(np.linalg.eigvalsh(H.toarray())[0]), H.shape[0]
    return float(spla.eigsh(H, k=1, which='SA', tol=1e-11)[0][0].real), H.shape[0]
def rotations(word):
    return [tuple(word[k:] + word[:k]) for k in range(len(word))]
if __name__ == '__main__':
    out = []
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    b = n + 1; nu = nu_sel(b); frozen = nu - 2*np.cos(np.pi/(b+1))
    W, R0, L0 = IDX['W'], IDX['R0'], IDX['L0']
    cell_words = []
    for j in range(n+1):                                   # walker position j: L0^j W R0^(n-j), all cyclic rotations
        cell_words += rotations([L0]*j + [W] + [R0]*(n-j))
    cell_words = sorted(set(cell_words))
    dead_words = [tuple([R0]*b)]
    for conv in ('contact', 'genuine'):
        for t in ((25.0, 5.0, 1.0, 0.2) if n <= 5 else (25.0, 1.0)):
            p = dict(lam=0.0, g=0.0, Delta=0.5, D=3.0, nu=nu, f=0.0, td=0.0, J=1.0, t=t, mu=0.0, K=2.0)
            h = bond10.bond(p, conv)
            for Nv in range(0, 6 if n <= 5 else 5):
                EA, dA = ground(cell_words, Nv, h); EB, dB = ground(dead_words, Nv, h)
                rec = dict(n=n, onsite=conv, t_over_J=t, Nv=Nv, vacancy_density=Nv/(b+Nv), dim=dA, eps_cell=EA-EB, frozen_value=float(frozen),
                           ratio_to_frozen=(EA-EB)/frozen)
                out.append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open(f'data/mobile_selector_n{n}.json', 'w'), indent=1)
