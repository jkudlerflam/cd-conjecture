"""Is the stated n = 2, 3 working point inside the weak-coupling window?

From the exact second-order bond H2 (rt_cell.second_order) we fit
  H2 = e0 + sum h_ij f_i^dag f_j + sum V_ijkl f_i^dag f_j^dag f_k f_l + rest
(Hilbert-Schmidt projection, so V is normal ordered at half filling), then run
self-consistent Hartree-Fock for the infinite chain at fixed logical chemical
potential eps, and count Fermi points and centres of the dressed bands.  This
is a mean-field diagnostic, not an infrared theorem: if HF already rearranges
the narrow band, the stated point is not perturbatively controlled.
Usage: python3 rt_ir.py n td g [g ...]
"""
import sys, json, itertools
import numpy as np
from rt_bond import P
from rt_cell import second_order, jw, bands, herm_part_analysis


def fit_two_body(H2, n):
    m = 2 * n
    c = jw(m)
    ops, labels = [np.eye(2 ** m)], [('id',)]
    for i in range(m):
        for j in range(m):
            ops.append(c[i].T @ c[j]); labels.append(('2', i, j))
    pairs = list(itertools.combinations(range(m), 2))
    for (i, j) in pairs:
        for (k, l) in pairs:
            ops.append(c[i].T @ c[j].T @ c[k] @ c[l]); labels.append(('4', i, j, k, l))
    A = np.array([o.reshape(-1) for o in ops]).T
    coef, *_ = np.linalg.lstsq(A, H2.reshape(-1), rcond=None)
    rest = H2 - (A @ coef).reshape(H2.shape)
    # HS projection on the quadratic span alone (half-filling normal order)
    A2 = A[:, :1 + m * m]
    c2, *_ = np.linalg.lstsq(A2, H2.reshape(-1), rcond=None)
    h_hs = c2[1:].reshape(m, m)
    V = np.zeros((m, m, m, m), complex)
    for cf, lb in zip(coef, labels):
        if lb[0] == '4':
            V[lb[1], lb[2], lb[3], lb[4]] = cf
    h_lin = coef[1:1 + m * m].reshape(m, m)     # quadratic part in the joint fit (vacuum normal order)
    V = 0.5 * (V + np.transpose(V, (2, 3, 0, 1)).conj())   # remove fit noise (H2 is Hermitian)
    h_lin = 0.5 * (h_lin + h_lin.conj().T)
    h_hs = 0.5 * (h_hs + h_hs.conj().T)
    return h_hs, h_lin, V, float(np.linalg.norm(rest, 2))


def hf_onebody(V, rho):
    """Mean field of sum V_ijkl f_i^dag f_j^dag f_k f_l for <f_a^dag f_b> = rho[a,b]."""
    h = (np.einsum('ijkl,jk->il', V, rho) - np.einsum('ijkl,ik->jl', V, rho)
         - np.einsum('ijkl,jl->ik', V, rho) + np.einsum('ijkl,il->jk', V, rho))
    herr = np.abs(h - h.conj().T).max() / max(1e-300, np.abs(V).max())
    assert herr < 1e-8, ('non-Hermitian mean field', herr)
    return 0.5 * (h + h.conj().T)


def run_hf(n, h_lin, V, eps, nk=4001, iters=200, mix=0.5):
    """Self-consistent HF with the quadratic part h_lin (vacuum normal order, so the
    full quartic V is added in mean field) at chemical potential -eps."""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    T0 = h_lin[:n, n:]
    mu0 = h_lin[:n, :n] + h_lin[n:, n:] + eps * np.eye(n)
    T, mu = T0.copy(), mu0.copy()
    for it in range(iters):
        Hk = mu[None] + T[None] * np.exp(1j * ks)[:, None, None] + T.conj().T[None] * np.exp(-1j * ks)[:, None, None]
        E, U = np.linalg.eigh(Hk)
        occ = (E < 0).astype(float)
        Pk = np.einsum('kab,kb,kcb->kac', U, occ, U.conj())          # P(k) = sum_occ u u^dag
        rho0 = np.mean(np.transpose(Pk, (0, 2, 1)), axis=0)           # <f^dag_{c,x} f_{c,y}>
        rho1 = np.mean(np.exp(1j * ks)[:, None, None] * np.transpose(Pk, (0, 2, 1)), axis=0)   # <f^dag_{c,x} f_{c+1,y}>
        rhom1 = rho1.conj().T                                          # <f^dag_{c+1,x} f_{c,y}>
        rho = np.block([[rho0, rho1], [rhom1, rho0]])
        hm = hf_onebody(V, rho)
        newmu = mu0 + hm[:n, :n] + hm[n:, n:]
        newT = T0 + hm[:n, n:]
        d = max(np.abs(newmu - mu).max(), np.abs(newT - T).max())
        mu = mix * newmu + (1 - mix) * mu
        T = mix * newT + (1 - mix) * T
        if d < 1e-10 * np.abs(T0).max():
            break
    return T, mu, it, d


def analyse(n, td, g):
    p = P(td=td, g=g)
    r = second_order(n, p)
    H2 = r['H2']
    h_hs, h_lin, V, rest = fit_two_body(H2, n)
    T_hs, mu_hs = h_hs[:n, n:], h_hs[:n, :n] + h_hs[n:, n:]
    sc = np.abs(T_hs).max()
    ev, cr, asym = herm_part_analysis(T_hs)
    eps = 0.25 * 2 * np.abs(ev).min() * sc          # the README's choice of logical field
    out = dict(n=n, td=td, g=g, T_scale=sc, rest_6body_over_T=rest / sc,
               V_max_over_T=float(np.abs(V).max() / sc), Rs_eigs=list(ev), eps_over_T=eps / sc)
    b0 = bands(T_hs, mu_hs + eps * np.eye(n), nk=40001)
    out['bare'] = dict(roots=[len(b['roots']) for b in b0], centres=[b['centre'] for b in b0],
                       fillings=[b['filling'] for b in b0])
    T, mu, it, d = run_hf(n, h_lin, V, eps)
    b1 = bands(T, mu, nk=40001)
    out['hf'] = dict(iters=it, conv=float(d / sc), roots=[len(b['roots']) for b in b1],
                     centres=[b['centre'] for b in b1], fillings=[b['filling'] for b in b1],
                     dmu_over_T=float(np.abs(mu - mu_hs - eps * np.eye(n)).max() / sc),
                     dT_over_T=float(np.abs(T - T_hs).max() / sc))
    return out


if __name__ == '__main__':
    n, td = int(sys.argv[1]), float(sys.argv[2])
    res = []
    for g in [float(x) for x in sys.argv[3:]]:
        o = analyse(n, td, g)
        print(json.dumps(o, default=float), flush=True)
        res.append(o)
    json.dump(res, open('ir_hf_n%d.json' % n, 'w'), indent=1, default=float)
