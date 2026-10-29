"""(3) C84 first order, second-order Gaussian structure, and small ED, from my11.py.
For p = 2, 3, 4 data sites (n = 1, 2, 3 active bits):
 (a) projected (K = infinity) single-cell matrices of C84 equal those of the C80
     rebuild (my12.py) with its last conversion deleted, restricted to last bit 0;
 (b) code = constrained vacuum, first-order block G^dag V G = 0;
 (c) exact Taylor coefficients of H_2(g) = -W^dag (H0 - 2E_v)^-1 W, W = V G:
     g^0 scalar, g^1 zero, g^2 Gaussian (span of 1 and c_i^dag c_j, Jordan-Wigner
     order left bits then right bits);
 (d) projected two-cell ED of H0 + lambda V: code levels vs 2E_v + lambda^2 eig(H_2);
 (e) finite-K ED (p = 2, 3) with the quasidegenerate formula (C84 eq. (15))."""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, splu
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'referee_C80_C82'))
import my11, my12

out = {}


def params(p, g=.07, K=np.inf):
    return my11.default_params(p, D=4., Delta=4 - 2 * np.cos(np.pi / (p + 2)), g=g, td=.003,
                               phiA=np.pi / (2 * p), phiB=0., psiA=np.pi, psiB=np.pi / 2, K=K)


def cell_matrix(p, pm, mod=my11, words=None, rules=None):
    words = words if words is not None else mod.legal_words(p)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(p)]
    left = 'E' if mod is my11 else 'E0'
    st, H, V = mod.closure(words, pm, bonds, left=left, project=True, rules=rules)
    assert len(st) == len(words)
    ix = {w: i for i, w in enumerate(st)}
    perm = [ix[w] for w in words]
    return H[np.ix_(perm, perm)].toarray()


def jw_ops(m):
    """annihilators c_0..c_{m-1}, occupation basis index = sum bit_i 2^(m-1-i)."""
    I2 = np.eye(2); Zs = np.diag([1., -1.]); a = np.array([[0., 1.], [0., 0.]])
    ops = []
    for i in range(m):
        mats = [Zs] * i + [a] + [I2] * (m - i - 1)
        M = np.array([[1.]])
        for x in mats:
            M = np.kron(M, x)
        ops.append(M)
    return ops


def gaussian_residual(H, m):
    c = jw_ops(m)
    basis = [np.eye(2 ** m)] + [c[i].T @ c[j] for i in range(m) for j in range(m)]
    A = np.array([b.ravel() for b in basis]).T
    coef, *_ = np.linalg.lstsq(A, H.ravel(), rcond=None)
    r = np.linalg.norm(A @ coef - H.ravel())
    tl = np.linalg.norm(H - np.trace(H) / len(H) * np.eye(len(H)))
    return dict(relative_to_traceless=float(r / tl), relative_to_full=float(r / np.linalg.norm(H)),
                traceless_over_full=float(tl / np.linalg.norm(H)))


for p in (2, 3, 4):
    n = p - 1
    rec = {}
    pm = params(p)
    W11 = my11.legal_words(p)
    # (a) parent comparison
    pm12 = my12.default_params(p, D=4., Delta=pm['Delta'], g=pm['g'], td=pm['td'], phiA=pm['phiA'], phiB=0.,
                               psiA=np.pi, psiB=np.pi / 2, K=np.inf)
    rules12 = [r for r in my12.move_rules(pm12) if not (r[0] == 'DA' and r[1] == 'E1')]
    W12 = [tuple('E0' if l == 'E' else l for l in w) for w in W11]
    assert all(w in set(my12.legal_words(p)) for w in W12)
    H11 = cell_matrix(p, pm)
    H12 = cell_matrix(p, pm12, my12, W12, rules12)
    rec['cell_matrix_vs_C80_minus_last_conversion'] = float(np.max(abs(H11 - H12)))
    # cell decomposition
    exc = np.array([sum(my11.excited(l) for l in w) for w in W11])
    S0 = np.flatnonzero(exc == 0); S1 = np.flatnonzero(exc == 1)
    assert np.max(abs(H11[np.ix_(S0, S1)])) == 0
    Ev = -2 * np.cos(np.pi / (p + 2))
    v = np.sqrt(2 / (p + 2)) * np.sin(np.pi * (np.arange(p + 1) + 1) / (p + 2))
    ix = {w: i for i, w in enumerate(W11)}
    bitsl = list(itertools.product((0, 1), repeat=n))
    G = np.zeros((len(W11), len(bitsl)))
    for k, b in enumerate(bitsl):
        for j in range(p + 1):
            G[ix[my11.encode_word(j, 'v', b)], k] = v[j]
    rec['code_eigen_residual'] = float(np.max(abs(H11 @ G - Ev * G)))
    e0 = np.linalg.eigvalsh(H11[np.ix_(S0, S0)])
    rec['vacuum_gap'] = float(sorted(set(np.round(e0, 10)))[1] - Ev)
    # two-cell port: V on (left marker, right first letter)
    prt = my11.rule_table(my11.port_rules(pm))
    nW = len(W11)
    rows, cols, vals = [], [], []
    for a, wa in enumerate(W11):
        for b, wb in enumerate(W11):
            for a2, b2, val in prt.get((wa[-1], wb[0]), []):
                za = wa[:-1] + (a2,); zb = (b2,) + wb[1:]
                if za in ix and zb in ix:  # projected: keep legal outputs only
                    rows.append(ix[za] * nW + ix[zb]); cols.append(a * nW + b); vals.append(val)
    V2 = sp.coo_matrix((vals, (rows, cols)), shape=(nW * nW, nW * nW)).tocsc()
    GG = np.kron(G, G)
    Wv = V2 @ GG
    rec['first_order_GVG'] = float(np.max(abs(GG.T @ Wv)))
    # (1,1) sector
    idx11 = (S1[:, None] * nW + S1[None, :]).ravel()
    rec['V_image_in_11_sector'] = float(np.linalg.norm(Wv[idx11]) / np.linalg.norm(Wv))
    W1 = Wv[idx11]
    h1 = H11[np.ix_(S1, S1)]
    pmg0 = params(p, g=0.); C1 = (cell_matrix(p, pm) - cell_matrix(p, pmg0))[np.ix_(S1, S1)] / pm['g']
    h10 = h1 - pm['g'] * C1
    e, U = np.linalg.eigh(h10)
    den = (e[:, None] + e[None, :] - 2 * Ev).ravel()
    rec['min_11_denominator'] = float(den.min())

    def R0(x):  # (h10 x I + I x h10 - 2Ev)^-1 x, x shape (m^2, k)
        m = len(e)
        X = x.reshape(m, m, -1)
        Y = np.einsum('ia,jb,ijk->abk', U.conj(), U.conj(), X, optimize=True)
        Y = Y / den.reshape(m, m)[:, :, None]
        return np.einsum('ia,jb,abk->ijk', U, U, Y, optimize=True).reshape(m * m, -1)

    def Cx(x):
        m = len(e)
        X = x.reshape(m, m, -1)
        return (np.einsum('ij,jbk->ibk', C1, X) + np.einsum('bj,ijk->ibk', C1, X)).reshape(m * m, -1)
    x = R0(W1)
    H20 = -W1.conj().T @ x
    H21 = x.conj().T @ Cx(x)
    y = Cx(x)
    H22 = -y.conj().T @ R0(y)
    rec['H2_g0_scalar_dev'] = float(np.max(abs(H20 - np.trace(H20) / len(H20) * np.eye(len(H20)))))
    rec['H2_g1_max'] = float(np.max(abs(H21)))
    rec['H22_hermiticity'] = float(np.max(abs(H22 - H22.conj().T)))
    # code basis order: (left bits, right bits) with left major -> occupation index on 2n modes
    rec['H22_gaussian_residual'] = gaussian_residual(H22, 2 * n)
    rec['H22_norm'] = float(np.linalg.norm(H22, 2))
    # exact H2 at this g
    eg, Ug = np.linalg.eigh(h1)
    deng = (eg[:, None] + eg[None, :] - 2 * Ev).ravel()
    m = len(eg)
    Y = np.einsum('ia,jb,ijk->abk', Ug.conj(), Ug.conj(), W1.reshape(m, m, -1), optimize=True) / deng.reshape(m, m)[:, :, None]
    H2g = -W1.conj().T @ np.einsum('ia,jb,abk->ijk', Ug, Ug, Y, optimize=True).reshape(m * m, -1)
    rec['H2_exact_vs_taylor_rel'] = float(np.linalg.norm(H2g - H20 - pm['g'] ** 2 * H22) / np.linalg.norm(pm['g'] ** 2 * H22))
    # (d) projected two-cell ED in the (0,0)+(1,1) legal space
    idx00 = (S0[:, None] * nW + S0[None, :]).ravel()
    sub = np.concatenate([idx00, idx11])
    Hs_ = sp.csr_matrix(H11); Id = sp.identity(nW, format='csr')
    H0 = (sp.kron(Hs_, Id) + sp.kron(Id, Hs_)).tocsr()
    eds = []
    h2e = np.linalg.eigvalsh(H2g)
    for lam in (.2, .1, .05):
        Hs = (H0[sub][:, sub] + lam * V2[sub][:, sub]).tocsr()
        qc = np.array([sum(my11.charge(l) for l in w) for w in W11])
        qs = (qc[:, None] + qc[None, :]).ravel()[sub]
        allev = []
        for Qv in sorted(set(qs)):
            ii = np.flatnonzero(qs == Qv)
            allev.append(np.linalg.eigvalsh(Hs[ii][:, ii].toarray()))
        ee = np.sort(np.concatenate(allev))[:4 ** n]
        err = float(np.max(abs(ee - (2 * Ev + lam ** 2 * h2e))))
        eds.append(dict(lam=lam, dim=len(sub), err_over_lam4=err / lam ** 4))
    rec['projected_ED'] = eds
    out[f'p{p}_n{n}'] = rec
    print(p, rec, flush=True)

# (e) finite-K ED with quasidegenerate formula, p = 2, 3
fk = []
for p in (2, 3):
    n = p - 1
    for K in (100., 1000.):
        pm = params(p, K=K)
        Ev = -2 * np.cos(np.pi / (p + 2))
        v = np.sqrt(2 / (p + 2)) * np.sin(np.pi * (np.arange(p + 1) + 1) / (p + 2))
        seeds = []
        for a in itertools.product((0, 1), repeat=n):
            for b in itertools.product((0, 1), repeat=n):
                seeds.append(my11.encode_word(0, 'v', a) + my11.encode_word(0, 'v', b))
        bonds, ports = my11.cell_bonds(p, 2)
        st, H0, V = my11.closure(seeds, pm, bonds, left='E', port_bonds=ports)
        ex = np.array([(sum(my11.excited(l) for l in w[:p + 1]), sum(my11.excited(l) for l in w[p + 1:])) for w in st])
        i00 = np.flatnonzero((ex[:, 0] == 0) & (ex[:, 1] == 0)); i11 = np.flatnonzero((ex[:, 0] == 1) & (ex[:, 1] == 1))
        h00 = H0[np.ix_(i00, i00)].toarray()
        e00, u00 = np.linalg.eigh(h00)
        Ea, Ga = e00[:4 ** n], u00[:, :4 ** n]
        Wa = V[np.ix_(i11, i00)] @ Ga
        h11 = H0[np.ix_(i11, i11)].tocsc()
        sol = {}
        for k, E in enumerate(Ea):
            lu = splu((h11 - E * sp.identity(len(i11), format='csc')).tocsc())
            sol[k] = lu.solve(Wa)
        H2 = np.zeros((4 ** n, 4 ** n), complex)
        for a_ in range(4 ** n):
            for b_ in range(4 ** n):
                H2[a_, b_] = -0.5 * (np.vdot(Wa[:, a_], sol[a_][:, b_]) + np.vdot(Wa[:, a_], sol[b_][:, b_]))
        rows = []
        for lam in (.2, .1, .05):
            A = (H0 + lam * V).tocsc()
            ee = np.sort(eigsh(A, k=4 ** n, sigma=2 * Ev - .3, which='LM')[0])
            pr = np.linalg.eigvalsh(np.diag(Ea) + lam ** 2 * H2)
            rows.append(dict(lam=lam, err_over_lam4=float(np.max(abs(ee - pr))) / lam ** 4))
        rec = dict(p=p, n=n, K=K, closure_dim=len(st), code_band_width=float(Ea[-1] - Ea[0]),
                   code_gap=float(e00[4 ** n] - Ea[-1]), runs=rows)
        print(rec, flush=True)
        fk.append(rec)
out['finite_K_ED'] = fk
json.dump(out, open(os.path.join(HERE, 'check_second11.json'), 'w'), indent=1)
