"""(2) Legal cells, code, first-order port, and the second-order Gaussian
matrices, from the referee rebuild only. Optional argument: path to a scratch
copy of codex/work_C33, used solely to evaluate C33's common-loop matrices
(fast_loop.fast) for the 1/2 comparison."""
import sys, os, json, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu, eigsh
import my15
from my15 import LETTERS, allowed

out = {}
EXC = lambda l: l[0] in 'DX' or l in ('MA', 'MB')

# (a) legal single cells with a frozen marker on the left
leg = []
for n in range(1, 6):
    data = [l for l in LETTERS if l[0] != 'M']
    cnt = 0; forms = set()
    for w in itertools.product(data, repeat=n):
        ok = allowed('M0', w[0]) and all(allowed(a, b) for a, b in zip(w, w[1:]))
        if not ok: continue
        for m in my15.MARKERS:
            if allowed(w[-1], m):
                cnt += 1
                forms.add(''.join(l[0] if l[0] in 'LR' else l[:2] for l in w) + m)
    exc_counts = sorted({sum(c in 'DX' for c in f[:-2]) + (f[-2:] in ('MA', 'MB')) for f in forms})
    leg.append(dict(n=n, legal_dim=cnt, formula=(5 * n + 3) * 2 ** n, orientation_forms=len(forms),
                    excitation_numbers=exc_counts))
    assert cnt == (5 * n + 3) * 2 ** n
out['legal_cells'] = leg; print(leg)

def code_words(n, bits):
    return [tuple([f'L{bits[i]}' for i in range(j)] + [f'R{bits[i]}' for i in range(j, n)] + ['M0']) for j in range(n + 1)]

# (b) vacuum path at K = infinity and at finite K (bits are spectators)
vac = []
for n in range(1, 7):
    p = my15.default_params(n)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(n)]
    st, H, _ = my15.closure(code_words(n, (0,) * n), p, bonds, project=True, left='M0')
    ev = np.linalg.eigvalsh(H.toarray())
    Ev = -2 * np.cos(np.pi / (n + 2))
    r = dict(n=n, zero_exc_legal_dim=len(st), E0=float(ev[0]), Ev=Ev, gap=float(ev[1] - ev[0]),
             gap_formula=2 * (np.cos(np.pi / (n + 2)) - np.cos(2 * np.pi / (n + 2))))
    # finite K, full orientation space (illegal strings included): spectra for all bit patterns
    pK = my15.default_params(n, K=37.)
    specs = []
    for bits in itertools.product((0, 1), repeat=n):
        stK, HK, _ = my15.closure(code_words(n, bits), pK, bonds, project=False, left='M0')
        specs.append(np.linalg.eigvalsh(HK.toarray()))
    r.update(finiteK_dim=len(stK), finiteK_dim_expected=2 ** n,
             finiteK_bit_spectator_spread=float(max(np.max(abs(s - specs[0])) for s in specs)),
             finiteK_E0=float(specs[0][0]), finiteK_gap=float(specs[0][1] - specs[0][0]),
             bound13_lower=Ev - max(n - 1, 0) / (37. - 4))
    vac.append(r); print(r)
out['vacuum'] = vac

# (c) two cells: first order and the g^2 Gaussian matrices
def two_cell(n, p, project):
    v = np.sqrt(2 / (n + 2)) * np.sin(np.pi * np.arange(1, n + 2) / (n + 2))
    L = 2 * (n + 1)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(L - 1)]
    seeds = []; G = {}
    for b1 in itertools.product((0, 1), repeat=n):
        for b2 in itertools.product((0, 1), repeat=n):
            idx = sum(b << i for i, b in enumerate(b1 + b2))
            for j1, w1 in enumerate(code_words(n, b1)):
                for j2, w2 in enumerate(code_words(n, b2)):
                    seeds.append(w1 + w2); G[w1 + w2] = (idx, v[j1] * v[j2])
    st, H, V = my15.closure(seeds, p, bonds, project=project, left='M0', with_port_at=(n, n + 1))
    return st, H, V, G

def fermions(m):
    out_ = []
    for j in range(m):
        c = np.zeros((2 ** m, 2 ** m))
        for b in range(2 ** m):
            if (b >> j) & 1: c[b ^ (1 << j), b] = (-1) ** bin(b & ((1 << j) - 1)).count('1')
        out_.append(c)
    return out_

def gaussian_residual(Hc, n):
    m = 2 * n; inds = [1 << i for i in range(m)]
    h = Hc[np.ix_(inds, inds)] - Hc[0, 0] * np.eye(m); cs = fermions(m)
    ref = Hc[0, 0] * np.eye(len(Hc), dtype=complex)
    for i in range(m):
        for j in range(m): ref = ref + h[i, j] * cs[i].T @ cs[j]
    return float(np.linalg.norm(Hc - ref, 2) / np.linalg.norm(Hc, 2))

def formula3(n, p):
    Ev = -2 * np.cos(np.pi / (n + 2)); a = p['Delta'] - Ev
    S = 2 * n + 1  # order: M, X_{n-1..0}, D_0..D_{n-1}
    Hg = np.zeros((S, S), complex)
    for i in range(S): Hg[i, i] = a + (p['D'] if i <= n else 0)
    for i in range(S - 1): Hg[i, i + 1] = Hg[i + 1, i] = -1
    Dj = lambda j: n + 1 + j
    lp = -p['td'] * np.exp(1j * p['psiB']); Hg[0, Dj(n - 1)] += lp; Hg[Dj(n - 1), 0] += np.conj(lp)
    I = np.eye(S); A = np.kron(Hg, I) + np.kron(I, Hg); Rr = np.linalg.inv(A)
    e = np.zeros(S * S); e[0 * S + Dj(0)] = 1; x = Rr @ e
    U = [np.kron(np.outer(I[Dj(i)], I[Dj(i)]), I) @ x for i in range(n)]
    Vv = [np.kron(I, np.outer(I[Dj(j)], I[Dj(j)])) @ x for j in range(n)]
    C = np.array([[U[i].conj() @ Rr @ Vv[j] for j in range(n)] for i in range(n)])
    Am = np.array([[U[i].conj() @ Rr @ U[j] + Vv[i].conj() @ Rr @ Vv[j] for j in range(n)] for i in range(n)])
    al4 = (np.sqrt(2 / (n + 2)) * np.sin(np.pi / (n + 2))) ** 4
    return -2 * al4 * C.real, -2j * al4 * Am.imag

c33 = None
if len(sys.argv) > 1:
    sys.path.insert(0, os.path.abspath(sys.argv[1])); sys.dont_write_bytecode = True
    cwd = os.getcwd(); os.chdir(os.path.abspath(sys.argv[1]))
    from fast_loop import fast as c33
    os.chdir(cwd)

sec = []
for n, phases in itertools.product([1, 2, 3], ['C75', 'real_common_loop']):
    Ev = -2 * np.cos(np.pi / (n + 2))
    ps = dict(psiA=np.pi, psiB=np.pi / 2) if phases == 'C75' else dict(psiA=np.pi / 2, psiB=0.)
    p = my15.default_params(n, D=6., Delta=4 + Ev + 0.3, td=0.07, g=0.5, **ps)
    st, H, V, G = two_cell(n, p, True)
    p0 = dict(p, g=0.)
    _, H0, _ = my15.closure(st, p0, [(-1, 0)] + [(i, i + 1) for i in range(2 * n + 1)], project=True, left='M0',
                            with_port_at=(n, n + 1))
    C = (H - H0) / p['g']
    N = len(st); Gm = np.zeros((N, 4 ** n), complex)
    for i, w in enumerate(st):
        if w in G: Gm[i, G[w][0]] += G[w][1]
    exc = np.array([[sum(EXC(l) for l in w[:n + 1]), sum(EXC(l) for l in w[n + 1:])] for w in st])
    code_res = float(np.max(abs(H0 @ Gm - 2 * Ev * Gm)))
    first_order = float(np.max(abs(Gm.conj().T @ (V @ Gm))))
    sel = np.flatnonzero((exc[:, 0] == 1) & (exc[:, 1] == 1))
    leak = float(abs(V @ Gm)[np.setdiff1d(np.arange(N), sel)].max())
    S = (V @ Gm)[sel]
    A0 = (H0[sel][:, sel] - 2 * Ev * sp.eye(len(sel))).tocsc(); lu = splu(A0)
    X = lu.solve(S); Y = C[sel][:, sel] @ X
    H22 = -Y.conj().T @ lu.solve(Y); H20 = -S.conj().T @ X
    inds = [1 << i for i in range(2 * n)]
    h = H22[np.ix_(inds, inds)] - H22[0, 0] * np.eye(2 * n)
    T = h[:n, n:]; mu = h[:n, :n] + h[n:, n:]
    Tf, muf = formula3(n, p)
    r = dict(n=n, phases=phases, two_cell_dim=N, source_dim=len(sel), code_residual=code_res, first_order_max=first_order,
             port_leak_outside_11=leak,
             H20_spread=float(np.max(abs(H20 - H20[0, 0] * np.eye(len(H20))))),
             H22_hermiticity=float(np.max(abs(H22 - H22.conj().T))),
             gaussian_residual=gaussian_residual(H22, n),
             T_vs_formula3=float(np.max(abs(T - Tf)) / np.max(abs(T))),
             mu_vs_formula3=float(np.max(abs(mu - muf)) / max(np.max(abs(mu)), 1e-300)),
             mu_over_T=float(np.max(abs(mu)) / np.max(abs(T))),
             T=[[str(np.round(x, 14)) for x in row] for row in T])
    if c33 is not None and phases == 'real_common_loop':
        # C33 fast(): neutral carrier loop phase loopflux-pi/2, charged carriers loop phase 0.
        # loopflux=pi/2 makes the loop common (real, entry -t_d), matching psi_B=0 here.
        T33, M33 = c33(n, D=p['D'], a=p['Delta'] - Ev, td=p['td'], loopflux=np.pi / 2)
        r.update(T_over_C33_T=[float(x) for x in (T.real / T33).ravel()],
                 T_minus_half_C33=float(np.max(abs(T - T33 / 2)) / np.max(abs(T))),
                 mu_max=float(np.max(abs(mu))), C33_mu_max=float(np.max(abs(M33))))
    sec.append(r); print({k: v for k, v in r.items() if k != 'T'})
out['second_order'] = sec

# (d) finite K: exact H2(K) at g fixed, compare with K = infinity; and ED in lambda at n=1,2
fk = []; ed = []
for n in [1, 2]:
    Ev = -2 * np.cos(np.pi / (n + 2))
    base = my15.default_params(n, D=4., Delta=4 + Ev, td=1e-2, g=0.04)
    st, H, V, G = two_cell(n, base, True)
    Ninf = len(st); Gi = np.zeros((Ninf, 4 ** n), complex)
    for i, w in enumerate(st):
        if w in G: Gi[i, G[w][0]] += G[w][1]
    exc = np.array([[sum(EXC(l) for l in w[:n + 1]), sum(EXC(l) for l in w[n + 1:])] for w in st])
    sel = np.flatnonzero((exc[:, 0] == 1) & (exc[:, 1] == 1)); S = (V @ Gi)[sel]
    H2inf = -S.conj().T @ splu((H[sel][:, sel] - 2 * Ev * sp.eye(len(sel))).tocsc()).solve(S)
    for K in [100., 1000., 10000.]:
        p = dict(base, K=K)
        st, H, V, G = two_cell(n, p, False)
        N = len(st)
        exc = np.array([[sum(EXC(l) for l in w[:n + 1]), sum(EXC(l) for l in w[n + 1:])] for w in st])
        bits = np.array([sum(my15.bit(l) << k for k, l in enumerate([x for x in w if my15.bit(x) is not None])) for w in st])
        z0 = np.flatnonzero((exc[:, 0] == 0) & (exc[:, 1] == 0))
        # zero-excitation blocks per bit pattern, orientation-sorted
        orient = lambda w: tuple(l[0] for l in w)
        blocks = [sorted(z0[bits[z0] == b], key=lambda i: orient(st[i])) for b in range(4 ** n)]
        Hb = [H[b][:, b].toarray() for b in blocks]
        spread = max(float(np.max(abs(x - Hb[0]))) for x in Hb)
        ev, U = np.linalg.eigh(Hb[0]); EK = ev[0]
        Gk = np.zeros((N, 4 ** n), complex)
        for b, bl in enumerate(blocks): Gk[bl, b] = U[:, 0]
        # fix sign so that the legal all-R-left... overlap with K=inf code is positive
        first = float(np.max(abs(Gk.conj().T @ (V @ Gk))))
        selK = np.flatnonzero((exc[:, 0] == 1) & (exc[:, 1] == 1))
        leak = float(abs(V @ Gk)[np.setdiff1d(np.arange(N), selK)].max())
        SK = (V @ Gk)[selK]
        H2K = -SK.conj().T @ splu((H[selK][:, selK] - EK * sp.eye(len(selK))).tocsc()).solve(SK)
        # code vectors share one sign convention (positive scalar ground vector), so compare directly
        diff = float(np.linalg.norm(H2K - H2inf, 2))
        r = dict(n=n, K=K, dim=N, bit_spectator_spread=spread, EK_minus_2Ev=float(EK - 2 * Ev),
                 first_order_max=first, port_leak_outside_11=leak, H2_diff=diff, K_times_diff=K * diff)
        fk.append(r); print(r)
        if K == 100.:
            for lam in [0.2, 0.1, 0.05]:
                Hl = (H + lam * V).tocsc()
                if N < 2000:
                    vals = np.linalg.eigvalsh(Hl.toarray())[:4 ** n]
                else:
                    vals = np.sort(eigsh(Hl, k=4 ** n, sigma=EK - 0.25, which='LM', ncv=80, tol=1e-12,
                                         return_eigenvectors=False))
                pred = EK + lam ** 2 * np.linalg.eigvalsh(H2K)
                e = float(np.max(abs(np.sort(vals) - pred)))
                ed.append(dict(n=n, K=K, lam=lam, dim=N, max_err=e, err_over_lam4=e / lam ** 4)); print(ed[-1])
out['finite_K'] = fk; out['ed'] = ed
json.dump(out, open('check_cells.json', 'w'), indent=1)
