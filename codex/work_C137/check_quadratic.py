"""Numerical checks of Theorem Q (PROOF.md section 4) on random quadratic qubit bonds.

Checks, for the h-independent sine window c_j = v_j v_{j+1}/s with no single-site terms
(the bond written in the parity-symmetric gauge, equal Z coefficients on both sites):
  (a) W' = a + sum_j v_j v_{j+1} h^cross_j + z sum_k v_k^2 Z_k >= e0(h)   [Lemma Q1]
  (b) lambda_min(sum_j c_j h_j) >= e0(h) - (1 - s) ||h_cross||          [Theorem Q]
  (c) ||h_cross|| <= M(h)
by the Majorana method for n up to 400 and by exact diagonalization for n <= 10."""
import warnings; warnings.filterwarnings('ignore')
from common import *
from free import *

rng = np.random.default_rng(137)
rows = []
worst_a = worst_b = worst_c = -np.inf
for trial in range(40):
    J = rng.normal(size=4); z = rng.normal() * rng.choice([0.1, 1, 3]); a0 = rng.normal()
    if trial < 5:   # include critical and near-critical TFIM, XX, XY in field
        J = [[-1, 0, 0, 0], [0.5, 0, 0, 0.5], [1, 0.3, -0.3, 0.6], [-1, 0, 0, 0], [0.7, 0.2, 0.1, 0.7]][trial]
        z = [-1, 0, 0.4, -1.03, 0.9][trial]
    h = a0 * np.eye(4) + J[0] * kron(PX, PX) + J[1] * kron(PX, PY) + J[2] * kron(PY, PX) + J[3] * kron(PY, PY) \
        + z * (kron(PZ, I2) + kron(I2, PZ)) / 2
    hc = J[0] * kron(PX, PX) + J[1] * kron(PX, PY) + J[2] * kron(PY, PX) + J[3] * kron(PY, PY)
    nc = float(np.max(np.abs(la.eigvalsh(hc))))
    e0 = e0_quadratic(h)
    gm = gauge_M(h)
    worst_c = max(worst_c, nc - gm['M_up'])
    for n in [2, 3, 4, 5, 8, 16, 33, 64, 128, 400]:
        v = np.sqrt(2 / (n + 1)) * np.sin(np.pi * np.arange(1, n + 1) / (n + 1)); s = np.cos(np.pi / (n + 1))
        c = v[:-1] * v[1:] / s
        lam = window_energy(h, c)
        # W' : bonds v_j v_{j+1} on the cross part, on-site v_k^2 z Z_k, scalar a0
        p = pauli_coeffs(h)
        bonds = [{k: v[j] * v[j + 1] * p[k] for k in ('XX', 'XY', 'YX', 'YY')} for j in range(n - 1)]
        A = build_A(n, bonds, v ** 2 * (p['ZI'] + p['IZ']))
        Wp = p['II'] - 0.25 * np.sum(np.abs(la.eigvalsh(1j * A)))
        da = e0 - Wp; db = (e0 - (1 - s) * nc) - lam
        worst_a = max(worst_a, da); worst_b = max(worst_b, db)
        rec = dict(trial=trial, n=n, e0=e0, M=gm['M_up'], hcross=nc, lam=lam, Wprime=Wp,
                   K_window=n * n * (e0 - lam) / gm['M_up'], bound=n * n * (1 - s) * nc / gm['M_up'])
        if n <= 10:
            fam = WindowFamily(h, n, udirs=(), cfix=c)
            rec['lam_ED'] = float(lowest(fam.matrix(np.zeros(0)), 1)[0][0])
            assert abs(rec['lam_ED'] - lam) < 1e-9, (rec['lam_ED'], lam)
        rows.append(rec)
write_json(HERE / 'results' / 'check_quadratic.json', rows)
print('max e0 - Wprime (should be <= 0):', worst_a)
print('max violation of Theorem Q (should be <= 0):', worst_b)
print('max ||h_cross|| - M (should be <= 0):', worst_c)
print('max n^2 (e0 - lam)/M:', max(r['K_window'] for r in rows), ' pi^2/2 =', np.pi ** 2 / 2)
