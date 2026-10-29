"""(5) C85 leading grammar correction and K suppression, from my14.py.
 (a) On every legal coloured ring word of lengths 4, 5, 6 (nonempty occupancy),
     build P B Q N_bad^-1 Q B P directly from the literal bond (penalty removed) and
     compare with diag(t^2 C + J^2 A + 2 J^2 lambda^2 Z) (C85 eq. (8)); also the mixed
     hopping-material block and the number of mismatches of each intermediate.
 (b) ED energy test on a ring of two p = 2 cells plus one vacancy, generic J, t,
     lambda: K (E_K - E_inf) against <psi| -(t^2 C + J^2 A + 2 J^2 lambda^2 Z)/J |psi>.
 (c) Illegal-material weight of the ground state on a ring of two p = 3 cells with
     nonzero bits plus two vacancies, at nonzero port, coloured and uncoloured,
     including the weight on squeezed words that violate only the boundary pin."""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my11, my14
from my14 import O, I, ext_allowed, is_vac

out = {}
cnum = {a: len(I[a]) - 1 for a in my11.LETTERS}


def counters(w):
    L = len(w)
    C = sum(cnum[w[i]] for i in range(L) if not is_vac(w[i]) and is_vac(w[(i + 1) % L]))
    A = sum(1 for i in range(L) if O[w[i - 1]] == 1 and w[i] == 'R0' and w[(i + 1) % L] in ('R0', 'R1', 'E'))
    Z = sum(1 for i in range(L) if O[w[i - 1]] == 1 and w[i] == 'E' and w[(i + 1) % L] == 'R0')
    return C, A, Z


def nbad(w):
    return sum(not ext_allowed(w[i], w[(i + 1) % len(w)]) for i in range(len(w)))


def legal_rings(L):
    res = []

    def rec(pref):
        if len(pref) == L:
            if ext_allowed(pref[-1], pref[0]) and any(not is_vac(x) for x in pref):
                res.append(tuple(pref))
            return
        for x in my14.LET:
            if not pref or ext_allowed(pref[-1], x):
                rec(pref + [x])
    rec([])
    return res


def transitions(h, ls, coloured=True):
    dd = len(ls)
    T = {}
    for col in range(dd * dd):
        for row in np.flatnonzero(abs(h[:, col]) > 0):
            if row != col:
                T.setdefault(col, []).append((row, h[row, col]))
    return T


pm = my11.default_params(3, D=3., g=.13, td=.21, phiA=.37, phiB=0., psiA=1.3, psiB=np.pi / 2, K=0., musel=.3,
                         field=.11)
J, t, lam = .7, 1.3, .37
h, ls = my14.bond(pm, lam, J, t, 0.)
ix = {l: i for i, l in enumerate(ls)}; dd = len(ls)
T = transitions(h, ls)
is_hop = lambda a, b: is_vac(a) or is_vac(b)
gc = {}
for L in (4, 5, 6):
    words = legal_rings(L)
    widx = {w: i for i, w in enumerate(words)}
    inter = {}  # z -> list of (word index, amplitude, kind)
    maxnb = 0
    for k, w in enumerate(words):
        for i in range(L):
            a, b = w[i], w[(i + 1) % L]
            for row, val in T.get(ix[a] * dd + ix[b], []):
                z = list(w); z[i] = ls[row // dd]; z[(i + 1) % L] = ls[row % dd]; z = tuple(z)
                nb = nbad(z)
                if nb:
                    maxnb = max(maxnb, nb)
                    inter.setdefault(z, []).append((k, val, 'T' if is_hop(a, b) else 'W'))
    Nw = len(words)
    Mfull = sp.lil_matrix((Nw, Nw), dtype=complex); Mmix = sp.lil_matrix((Nw, Nw), dtype=complex)
    for z, lst in inter.items():
        nb = nbad(z)
        for (k1, v1, s1) in lst:
            for (k2, v2, s2) in lst:
                Mfull[k1, k2] += np.conj(v1) * v2 / nb
                if s1 == 'T' and s2 == 'W':
                    Mmix[k1, k2] += np.conj(v1) * v2 / nb
    Mfull = Mfull.tocsr()
    pred = np.array([t ** 2 * c + J ** 2 * a + 2 * J ** 2 * lam ** 2 * zz for c, a, zz in map(counters, words)])
    D = Mfull - sp.diags(pred)
    gc[L] = dict(legal_words=Nw, illegal_intermediates=len(inter), max_mismatch_of_intermediate=maxnb,
                 max_abs_discrepancy=float(abs(D).max()), mixed_block_max=float(abs(Mmix.tocsr()).max()) if Mmix.nnz else 0.)
    print(L, gc[L], flush=True)
out['leading_correction_operator'] = gc


# generic ring closure from a bond matrix
def ring_closure(seed, h, ls, Kpen, coloured=True, project=False, cap=500000):
    dd = len(ls); ixx = {l: i for i, l in enumerate(ls)}
    T = transitions(h, ls)
    diag = np.real(np.diag(h))
    L = len(seed)
    allowedf = (lambda a, b: ext_allowed(a, b)) if coloured else (
        lambda a, b: True if (a == 'V' or b == 'V') else my11.allowed(a, b))
    nbf = lambda w: sum(not allowedf(w[i], w[(i + 1) % L]) for i in range(L))
    states = [tuple(seed)]; idx = {states[0]: 0}
    r, c, v = [], [], []
    cur = 0
    while cur < len(states):
        w = states[cur]
        dg = 0.
        for i in range(L):
            a, b = w[i], w[(i + 1) % L]
            col = ixx[a] * dd + ixx[b]
            dg += diag[col]
            for row, val in T.get(col, []):
                z = list(w); z[i] = ls[row // dd]; z[(i + 1) % L] = ls[row % dd]; z = tuple(z)
                if project and nbf(z):
                    continue
                if z not in idx:
                    idx[z] = len(states); states.append(z)
                    if len(states) > cap:
                        raise RuntimeError('cap')
                r.append(idx[z]); c.append(cur); v.append(val)
        r.append(cur); c.append(cur); v.append(dg)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsc()
    pen = np.array([nbf(w) for w in states], float)
    return states, H, pen


def lowest(H, k=4):
    if H.shape[0] < 1500:
        e, U = np.linalg.eigh(H.toarray())
        return e[:k], U[:, :k]
    sigma = -float(abs(H).sum(axis=1).max()) - 1
    e, U = eigsh(H, k=k, sigma=sigma, which='LM', tol=1e-13)
    o = np.argsort(e)
    return e[o], U[:, o]


# (b) energy test
seed = ('R0', 'R1', 'E', 'R0', 'R0', 'E', 'V0')
pmb = my11.default_params(2, D=3., g=.13, td=.21, phiA=.37, phiB=0., psiA=1.3, psiB=np.pi / 2, K=0., musel=.3,
                          field=.11)
J, t, lam = .7, 1.3, .3
h0, ls = my14.bond(pmb, lam, J, t, 0.)
stP, HP, _ = ring_closure(seed, h0, ls, 0., project=True)
eP, UP = lowest(HP, 3)
psi = UP[:, 0]
corr = np.array([-(t ** 2 * c + J ** 2 * a + 2 * J ** 2 * lam ** 2 * zz) / J for c, a, zz in map(counters, stP)])
slope_pred = float(np.vdot(psi, corr * psi).real)
stF, HF, pen = ring_closure(seed, h0, ls, 0.)
en = []
for K in (3e2, 1e3, 3e3, 1e4):
    e, _ = lowest(HF + sp.diags(J * K * pen), 2)
    en.append(dict(K=K, K_times_shift=float(K * (e[0] - eP[0]))))
out['energy_test'] = dict(seed=seed, J=J, t=t, lam=lam, projected_dim=len(stP), full_dim=len(stF),
                          projected_gap=float(eP[1] - eP[0]), predicted_slope=slope_pred, runs=en)
print(out['energy_test'], flush=True)

# (c) leakage with nonzero bits and nonzero port
def colour(seed):
    seed = list(seed)
    for i, l in enumerate(seed):
        if is_vac(l):
            seed[i] = f'V{O[seed[i - 1]]}'
    seed = tuple(seed)
    assert nbad(seed) == 0, seed
    return seed


c1 = my11.encode_word(0, 'v', (1, 1)); c2 = my11.encode_word(1, 'v', (0, 1))
d1 = my11.encode_word(0, 'v', (1,)); d2 = my11.encode_word(2, 'v', (1,))
geoms = dict(two_p3_cells_one_vacancy=(3, 0., colour(c1[:2] + ('V',) + c1[2:] + c2)),
             two_p2_cells_two_vacancies=(2, .1, colour(d1[:1] + ('V',) + d1[1:] + d2 + ('V',))),
             two_p2_cells_one_vacancy=(2, .1, colour(d1 + d2[:2] + ('V',) + d2[2:])))
leak = []


def squeezed_bad(w):
    m = [l for l in w if not (l == 'V' or is_vac(l))]
    bad = [(m[i], m[(i + 1) % len(m)]) for i in range(len(m)) if not my11.allowed(m[i], m[(i + 1) % len(m)])]
    return len(bad), all(a == 'L1' and b in my11.MARKERS for a, b in bad) and len(bad) > 0


for gname, (pp, lamc, seed) in geoms.items():
    pmc = my11.default_params(pp, D=4., g=.07, td=.003, phiA=np.pi / (2 * pp), phiB=0., psiA=np.pi,
                              psiB=np.pi / 2, K=0.)
    for coloured in (True, False):
        for J, Ks in [(1., (1e2, 1e3, 1e4)), (.1, (1e3, 1e4, 1e5))]:
            if not coloured and J != 1.:
                continue
            sd = seed if coloured else tuple('V' if is_vac(l) else l for l in seed)
            for K in Ks:
                pk = dict(pmc); pk['K'] = K
                hk, ls = my14.bond(pk, lamc, J, 1., 0., coloured)
                st, H, pen = ring_closure(sd, hk, ls, 0., coloured=coloured)
                e, U = lowest(H, 1)
                v = U[:, 0]
                sb = [squeezed_bad(w) for w in st]
                ill = np.array([x[0] > 0 for x in sb]); pin = np.array([x[1] for x in sb])
                res = float(np.linalg.norm(H @ v - e[0] * v))
                leak.append(dict(geometry=gname, lam=lamc, seed=sd, coloured=coloured, J=J, K=K, JK=J * K, dim=len(st),
                                 illegal_material_weight=float(np.sum(abs(v[ill]) ** 2)),
                                 pin_only_weight=float(np.sum(abs(v[pin]) ** 2)),
                                 pin_only_states=int(pin.sum()), residual=res))
                print(leak[-1], flush=True)
out['leakage'] = leak
json.dump(out, open(os.path.join(HERE, 'check_grammar14.json'), 'w'), indent=1)
