#!/usr/bin/env python3
"""Toy infrared test for the encoded free-fermion target (contrib/ir_control/README.md).

Logical chain: n complex fermion flavours per site, written as a spin chain on C^(2^n)
by Jordan-Wigner (fermion parity is exact, every bond term is parity even).  The bond
operator h (d^2 x d^2) is passed to the two-site iDMRG of contrib/counterexample_search.

Models
  free   : H_* = sum_a -t_a (e^{-i th_a} c+_{j,a} c_{j+1,a} + h.c.) - mu_a n_{j,a}
  +rand  : H_* + g * (random parity-even Hermitian two-site operator, fixed seed, unit norm)
  sym    : the symmetric Majorana point th_a = k_a = pi/2, t_a = 1 (all Fermi points at 0, pi);
           its bond is (1/2) T_j with T_j = sum_{a=1}^{2n} i gamma_{j,a} gamma_{j+1,a}
  +gn    : + g T_j^2, the SO(2n)-symmetric quartic of C8 RESULT.md Sec. 11 (mass channel
           coefficient a - 2b = g with b = 0, so g < 0 is the attractive sign)
  leak   : each site also carries an aux qubit (d = 2^n * 2) with H_cell = Delta |1><1|_aux,
           and the bond is lam * [h_* on the code] + lam * (random even leakage terms that
           flip aux qubits), so the code chain is reached only through second-order
           Schrieffer-Wolff.  The logical target appears at first order in lam.

Central charge: finite-entanglement scaling, S = (c/6) log xi + const over the chi sweep
(Pollmann et al. 2009).  Calibrated on free chains below.

usage: python3 toy_ir.py --case NAME [--g 0.3] [--seed 1] [--chis 16 24 32 48 64 96 128]
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '2')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'counterexample_search'))
sys.path.insert(0, HERE)  # on the cluster idmrg.py is copied next to this file
import idmrg  # noqa: E402
from idmrg import run  # noqa: E402
from scipy.sparse.linalg import eigsh as _eigsh, ArpackNoConvergence  # noqa: E402


def _robust_eigsh(op, k=2, which='SA', v0=None, tol=0.0, ncv=None, **kw):
    """ARPACK occasionally stalls on near-degenerate gapless spectra; retry, then fall back."""
    try:
        return _eigsh(op, k=k, which=which, v0=v0, tol=tol, ncv=ncv, **kw)
    except ArpackNoConvergence:
        try:
            return _eigsh(op, k=k, which=which, v0=v0, tol=max(tol, 1e-8), ncv=min(op.shape[0] - 1, 60), maxiter=20000)
        except ArpackNoConvergence as e:
            if len(e.eigenvalues):
                return e.eigenvalues, e.eigenvectors
            return _eigsh(op, k=1, which=which, v0=v0, tol=1e-6, ncv=min(op.shape[0] - 1, 80), maxiter=50000)


idmrg.eigsh = _robust_eigsh


def local_ops(n):
    """Annihilators b_a on C^(2^n) with intra-site Jordan-Wigner, and the parity P."""
    I2 = np.eye(2); Z = np.diag([1., -1.]); a = np.array([[0., 1.], [0., 0.]])
    ops = []
    for k in range(n):
        f = [Z] * k + [a] + [I2] * (n - k - 1)
        m = f[0]
        for x in f[1:]:
            m = np.kron(m, x)
        ops.append(m.astype(complex))
    P = np.eye(1)
    for _ in range(n):
        P = np.kron(P, Z)
    return ops, P.astype(complex)


def free_bond(n, t, th, mu, aux=False):
    """Two-site bond of H_* (on-site chemical potential split evenly over the two sites)."""
    b, P = local_ops(n)
    d = 2 ** n
    I = np.eye(d)
    h = np.zeros((d * d, d * d), complex)
    for k in range(n):
        bd = b[k].conj().T
        # c+_{j,k} c_{j+1,k} = (b+_k P)_j (x) (b_k)_{j+1}
        hop = np.kron(bd @ P, b[k])
        h += -t[k] * (np.exp(-1j * th[k]) * hop + np.exp(1j * th[k]) * hop.conj().T)
        nk = bd @ b[k]
        h += -mu[k] * 0.5 * (np.kron(nk, I) + np.kron(I, nk))
    return h


def random_even_bond(n, seed, extra_dim=1):
    """Random Hermitian two-site operator, even under P (x) P, unit operator norm."""
    _, P = local_ops(n)
    P = np.kron(P, np.eye(extra_dim))
    d = P.shape[0]
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(d * d, d * d)) + 1j * rng.normal(size=(d * d, d * d))
    X = (X + X.conj().T) / 2
    PP = np.kron(P, P)
    X = (X + PP @ X @ PP) / 2              # keep the parity-even part
    X -= np.trace(X) / (d * d) * np.eye(d * d)
    return X / np.linalg.norm(X, 2)


def fermi_data(t, th, mu):
    k = np.arccos(-np.asarray(mu) / (2 * np.asarray(t)))
    pts = np.concatenate([np.asarray(th) + k, np.asarray(th) - k])
    v = 2 * np.asarray(t) * np.sin(k)
    X = np.concatenate([pts, -pts]) % (2 * np.pi)
    Xs = np.sort(X)
    gaps = np.diff(np.concatenate([Xs, [Xs[0] + 2 * np.pi]]))
    return dict(k=k.tolist(), p=(pts % (2 * np.pi)).tolist(), v=v.tolist(), delta2=float(gaps.min()))


def c_fit(res, drop=1):
    S = np.array([0.5 * (r['S_AB'] + r['S_BA']) for r in res])
    xi = np.array([r['xi'] for r in res])
    m = np.isfinite(xi) & (xi > 0)
    S, xi = S[m][drop:], xi[m][drop:]
    if len(S) < 2:
        return float('nan')
    A = np.stack([np.log(xi), np.ones_like(xi)], 1)
    return float(6 * np.linalg.lstsq(A, S, rcond=None)[0][0])


# Two-flavour incommensurate target: Fermi points spread out, velocities about 1 to 1.6
# chosen by a random search maximizing delta_2 (0.597) and the four-index margin (0.434)
INC = dict(t=[1.0, 1.0], th=[-2.617, 0.298], mu=[-2 * np.cos(0.962), -2 * np.cos(0.742)])
SYM = dict(t=[1.0, 1.0], th=[np.pi / 2, np.pi / 2], mu=[0.0, 0.0])


def build(case, g, seed, Delta=1.0, lam=0.25):
    n = 2
    if case.startswith('xx'):            # n = 1 calibration, c = 1
        h = free_bond(1, [1.0], [0.3], [-2 * np.cos(1.1)])
        return h, 2, fermi_data([1.0], [0.3], [-2 * np.cos(1.1)])
    if case in ('inc_gn', 'sym_gn'):
        par = INC if case == 'inc_gn' else SYM
        T = 2 * free_bond(n, **SYM)                      # sum_a i gamma gamma on the bond
        h = free_bond(n, **par) + g * (T @ T)
        return (h + h.conj().T) / 2, 2 ** n, fermi_data(**par)
    if case in ('inc', 'inc_rand'):
        h = free_bond(n, **INC)
        fd = fermi_data(**INC)
    elif case in ('sym', 'sym_rand'):
        h = free_bond(n, **SYM)
        fd = fermi_data(**SYM)
    elif case in ('inc_leak', 'sym_leak'):
        par = INC if case == 'inc_leak' else SYM
        fd = fermi_data(**par)
        d0 = 2 ** n
        hl = free_bond(n, **par)                        # target on the code
        Pc = np.diag([1., 0.])                          # aux in |0> = code
        # embed: logical (x) aux per site, ordering (log1, aux1, log2, aux2)
        def emb(op2):                                    # op on log1 (x) log2 -> full, aux untouched
            op = op2.reshape(d0, d0, d0, d0)
            I2 = np.eye(2)
            full = np.einsum('abcd,ef,gh->aebgcfdh', op, I2, I2)
            return full.reshape(4 * d0 * d0, 4 * d0 * d0)
        def proj_code(op_full):
            Pfull = np.kron(np.kron(np.eye(d0), Pc), np.kron(np.eye(d0), Pc))
            return Pfull @ op_full @ Pfull
        L = random_even_bond(n, seed, extra_dim=2)       # generic even operator on (log (x) aux)^2
        Pfull = np.kron(np.kron(np.eye(d0), Pc), np.kron(np.eye(d0), Pc))
        L = L - Pfull @ L @ Pfull                        # no first-order action on the code
        L = L / np.linalg.norm(L, 2)
        Hc = Delta * np.kron(np.eye(d0), np.diag([0., 1.]))
        dd = 2 * d0
        Icell = np.eye(dd)
        # the target acts on the logical factor whatever the aux state; L has no code-code block,
        # so the first-order effective chain is exactly lam * H_*, and L enters at order lam^2/Delta
        h = lam * emb(hl) + lam * g * L + 0.5 * (np.kron(Hc, Icell) + np.kron(Icell, Hc))
        return (h + h.conj().T) / 2, dd, fd
    else:
        raise ValueError(case)
    if case.endswith('_rand'):
        h = h + g * random_even_bond(n, seed)
    return (h + h.conj().T) / 2, 2 ** n, fd


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True)
    ap.add_argument('--g', type=float, default=0.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--lam', type=float, default=0.25)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64, 96, 128])
    ap.add_argument('--steps', type=int, default=200)
    ap.add_argument('--out', default=os.path.join(HERE, 'toy_ir_results.jsonl'))
    a = ap.parse_args()
    h, d, fd = build(a.case, a.g, a.seed, lam=a.lam)
    t0 = time.time()
    res = run(h, d, a.chis, steps_per_chi=a.steps, log=True)
    rec = dict(case=a.case, g=a.g, seed=a.seed, lam=a.lam, d=d, fermi=fd, c_fit=c_fit(res),
               c_fit_last4=c_fit(res[-4:], drop=0),
               chi=[r['chi'] for r in res], S=[0.5 * (r['S_AB'] + r['S_BA']) for r in res],
               xi=[r['xi'] for r in res], e=[r['e_site'] for r in res], dimer=[r['dimer'] for r in res],
               secs=time.time() - t0)
    print(json.dumps(rec))
    with open(a.out, 'a') as f:
        f.write(json.dumps(rec) + '\n')
