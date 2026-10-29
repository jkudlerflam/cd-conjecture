"""(6) Vacancy lift in the zero-excitation sector at zero port (data bits are
spectators and are dropped). Bond: -t vacancy hopping on every bond, plus J times
the material bond on adjacent occupied pairs only, as in COUNTEREXAMPLE_D24.md.
15-letter reduction: letters L, R, M (= M_0), moves RR<->LR, RM<->LM, penalty K on
(R,L). 18-letter parent: letters L, R, d (vacuum token), e (empty marker),
m (vacuum marker), moves dR<->Ld, de<->Lm, arrow-grammar penalty.
Observable: ground-state expectation of the number of illegal MATERIAL pairs
(neighbours after deleting vacancies), starting from a legal seed."""
import json, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def model(kind):
    if kind == 15:
        moves = {('R', 'R'): [('L', 'R')], ('L', 'R'): [('R', 'R')], ('R', 'M'): [('L', 'M')], ('L', 'M'): [('R', 'M')]}
        bad = lambda a, b: a == 'R' and b == 'L'
        cell = lambda n: ['L'] + ['R'] * (n - 1) + ['M']          # wall after first site
    else:
        AR = {'L': 'LL', 'R': 'RR', 'd': 'LR', 'e': 'RL', 'm': 'LL'}
        moves = {('d', 'R'): [('L', 'd')], ('L', 'd'): [('d', 'R')], ('d', 'e'): [('L', 'm')], ('L', 'm'): [('d', 'e')]}
        bad = lambda a, b: AR[a][1] != AR[b][0]
        cell = lambda n: ['L', 'd'] + ['R'] * (n - 2) + ['e']
    return moves, bad, cell

def material_bad(w, bad):
    m = [x for x in w if x != 'V']
    return sum(bad(m[i], m[(i + 1) % len(m)]) for i in range(len(m)))

def run(kind, n, cells, nv, J, K, t=1.):
    moves, bad, cell = model(kind)
    mat = sum((cell(n) for _ in range(cells)), [])
    L = len(mat) + nv
    # spread vacancies roughly evenly
    seed = []; vi = set(int(round(i * L / nv)) % L for i in range(nv)) if nv else set()
    while len(vi) < nv: vi.add(max(set(range(L)) - vi))
    it = iter(mat)
    for i in range(L): seed.append('V' if i in vi else next(it))
    states = [tuple(seed)]; ix = {states[0]: 0}; r_, c_, v_ = [], [], []; cur = 0
    while cur < len(states):
        w = states[cur]
        def emit(z, val):
            z = tuple(z)
            if z not in ix: ix[z] = len(states); states.append(z)
            r_.append(ix[z]); c_.append(cur); v_.append(val)
        diag = 0.
        for i in range(L):
            j = (i + 1) % L; a, b = w[i], w[j]
            if (a == 'V') != (b == 'V'):
                z = list(w); z[i], z[j] = b, a; emit(z, -t)
            elif a != 'V':
                for a2, b2 in moves.get((a, b), []):
                    z = list(w); z[i], z[j] = a2, b2; emit(z, -J)
                if bad(a, b): diag += J * K
        if diag: r_.append(cur); c_.append(cur); v_.append(diag)
        cur += 1
        if len(states) > 300000: raise RuntimeError('cap')
    N = len(states)
    H = sp.coo_matrix((v_, (r_, c_)), shape=(N, N)).tocsc()
    if N < 1500:
        ev, U = np.linalg.eigh(H.toarray()); E, psi = ev[0], U[:, 0]
    else:
        ev, U = eigsh(H, k=1, which='SA', tol=1e-12); E, psi = ev[0], U[:, 0]
    nb = np.array([material_bad(w, bad) for w in states])
    return dict(kind=kind, n=n, cells=cells, vacancies=nv, J_over_t=J, K=K, JK_over_t=J * K, dim=N,
                reachable_illegal_states=int(np.sum(nb > 0)), E0=float(E),
                ground_material_bad_pairs=float(np.sum(abs(psi) ** 2 * nb)),
                ground_illegal_weight=float(np.sum(abs(psi[nb > 0]) ** 2)))

out = []
for kind, (n, cells, nv), J, K in itertools.product([15, 18], [(3, 2, 1), (2, 2, 6)], [1., 1e-2], [10., 100., 1e3, 1e4, 1e5]):
    r = run(kind, n, cells, nv, J, K); out.append(r); print(r, flush=True)
json.dump(out, open('check_vacancy.json', 'w'), indent=1)
