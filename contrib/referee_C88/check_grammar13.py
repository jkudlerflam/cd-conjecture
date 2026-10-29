"""(3) Leading grammar correction (PROOF eq. (10)) and K^-2 suppression, from my13.py.
 (a) Every legal coloured ring word of lengths 4, 5, 6 with N > 0: build
     P B Q N_bad^-1 Q B P from the literal bond (penalty stripped) and compare with
     diag(t^2 C + J^2 A + 2 J^2 lambda^2 Z); mixed hopping-material block; mismatch
     number of each intermediate. Generic J, t, lambda, g, t_d, field, selector.
 (b) Energy test with nonzero port: ring of two p = 2 cells and one vacancy,
     K (E_K - E_inf) against <psi_inf| -(t^2 C + J^2 A + 2 J^2 lambda^2 Z)/J |psi_inf>.
 (c) Ground-state illegal-material weight in own geometries with nonzero bits,
     one and two vacancies and nonzero port, coloured and uncoloured."""
import sys, os, json
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13
from my13 import LET, allowed, is_vac, counters, nbad

out = {}
pm = my13.params(3, D=3., g=.13, td=.21, phiA=.37, phiB=0., psiA=1.3, psiB=np.pi / 2, K=0., musel=.3, field=.11)
J, t, lam = .7, 1.3, .37
h, ls = my13.bond(pm, lam, J, t, 0.)
dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
T = my13.transitions(h)


def legal_rings(L):
    res = []

    def rec(pref):
        if len(pref) == L:
            if allowed(pref[-1], pref[0]) and any(not is_vac(x) for x in pref):
                res.append(tuple(pref))
            return
        for x in LET:
            if not pref or allowed(pref[-1], x):
                rec(pref + [x])
    rec([])
    return res


gc = {}
for L in (4, 5, 6):
    words = legal_rings(L)
    inter = {}
    maxnb = 0
    for k, w in enumerate(words):
        for i in range(L):
            a, b = w[i], w[(i + 1) % L]
            for row, val in T[ix[a] * dd + ix[b]]:
                zz = list(w); zz[i] = ls[row // dd]; zz[(i + 1) % L] = ls[row % dd]; zz = tuple(zz)
                nb = nbad(zz)
                if nb:
                    maxnb = max(maxnb, nb)
                    inter.setdefault(zz, []).append((k, val, 'T' if (is_vac(a) or is_vac(b)) else 'W'))
    Nw = len(words)
    r, c, v, rm, cm, vm = [], [], [], [], [], []
    for zz, lst in inter.items():
        nb = nbad(zz)
        for (k1, v1, s1) in lst:
            for (k2, v2, s2) in lst:
                r.append(k1); c.append(k2); v.append(np.conj(v1) * v2 / nb)
                if s1 != s2:
                    rm.append(k1); cm.append(k2); vm.append(np.conj(v1) * v2 / nb)
    M = sp.coo_matrix((v, (r, c)), shape=(Nw, Nw)).tocsr()
    Mx = sp.coo_matrix((vm, (rm, cm)), shape=(Nw, Nw)).tocsr() if vm else sp.csr_matrix((Nw, Nw))
    Mx.sum_duplicates()
    pred = np.array([t ** 2 * cc + J ** 2 * aa + 2 * J ** 2 * lam ** 2 * zc for cc, aa, zc in map(counters, words)])
    D = M - sp.diags(pred)
    gc[L] = dict(legal_words=Nw, illegal_intermediates=len(inter), max_mismatch_of_intermediate=maxnb,
                 max_abs_discrepancy=float(abs(D).max()), mixed_block_max=float(abs(Mx).max()) if Mx.nnz else 0.,
                 words_with_C_positive=int(sum(counters(w)[0] > 0 for w in words)),
                 words_with_pin=int(sum(my13.F(w) > 0 for w in words)))
    print(L, gc[L], flush=True)
out['leading_correction_operator'] = gc


def lowest(H, k=2):
    if H.shape[0] < 2500:
        e, U = np.linalg.eigh(H.toarray())
        return e[:k], U[:, :k]
    sigma = -float(abs(H).sum(axis=1).max()) - 1
    e, U = eigsh(H, k=k, sigma=sigma, which='LM', tol=1e-13)
    o = np.argsort(e)
    return e[o], U[:, o]


# (b) energy test
seed = my13.colour(('R0', 'R1', 'E', 'R0', 'R0', 'E', 'V'))
pmb = my13.params(2, D=3., g=.13, td=.21, phiA=.37, phiB=0., psiA=1.3, psiB=np.pi / 2, K=0., musel=.3, field=.11)
J, t, lam = .7, 1.3, .3
h0, ls = my13.bond(pmb, lam, J, t, 0.)
stP, HP, _, _ = my13.ring_closure(seed, h0, ls, project=True)
eP, UP = lowest(HP, 3)
psi = UP[:, 0]
corr = np.array([-(t ** 2 * c + J ** 2 * a + 2 * J ** 2 * lam ** 2 * zz) / J for c, a, zz in map(counters, stP)])
slope = float(np.vdot(psi, corr * psi).real)
stF, HF, pen, mb = my13.ring_closure(seed, h0, ls)
en = []
for K in (3e2, 1e3, 3e3, 1e4, 3e4):
    e, _ = lowest(HF + sp.diags(J * K * pen), 2)
    en.append(dict(K=K, K_times_shift=float(K * (e[0] - eP[0]))))
out['energy_test'] = dict(seed=seed, J=J, t=t, lam=lam, projected_dim=len(stP), full_dim=len(stF),
                          projected_gap=float(eP[1] - eP[0]), predicted_K_times_shift=slope, runs=en)
print(out['energy_test'], flush=True)

# (c) leakage
c1 = my13.encode(0, 'v', (1, 1)); c2 = my13.encode(1, 'v', (0, 1))
d1 = my13.encode(0, 'v', (1,)); d2 = my13.encode(2, 'v', (1,))
geoms = dict(two_p3_cells_bits_11_01_one_vacancy=(3, 0., c1[:2] + ('V',) + c1[2:] + c2),
             two_p2_cells_two_vacancies_port=(2, .1, d1[:1] + ('V',) + d1[1:] + d2 + ('V',)),
             two_p2_cells_one_vacancy_port=(2, .1, d1 + d2[:2] + ('V',) + d2[2:]))
leak = []
for gname, (pp, lamc, seed0) in geoms.items():
    pmc = my13.params(pp, K=0.)
    for coloured in (True, False):
        sd = my13.colour(seed0) if coloured else seed0
        for Jc in ((1., .1) if coloured else (1.,)):
            hk, ls = my13.bond(pmc, lamc, Jc, 1., 0., coloured)
            st, H, pen, mb = my13.ring_closure(sd, hk, ls, coloured=coloured)
            for JK in (1e2, 1e3, 1e4):
                K = JK / Jc
                A = H + sp.diags(Jc * K * pen) if coloured else H + sp.diags(Jc * K * pen)
                e, U = lowest(A, 1)
                vv = U[:, 0]
                rec = dict(geometry=gname, lam=lamc, coloured=coloured, J=Jc, K=K, JK=JK, dim=len(st),
                           illegal_material_weight=float(np.sum(abs(vv[mb > 0]) ** 2)),
                           pin_weight=float(np.sum(abs(vv[np.array([my13.F(w) for w in st]) > 0]) ** 2)),
                           residual=float(np.linalg.norm(A @ vv - e[0] * vv)))
                leak.append(rec); print(rec, flush=True)
out['leakage'] = leak
json.dump(out, open(os.path.join(HERE, 'check_grammar13.json'), 'w'), indent=1, default=str)
