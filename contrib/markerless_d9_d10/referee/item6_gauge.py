"""Item 6(c),(e): does the certificate's scalar graph + C84-form formula describe this bond?

Literal letter model (second_order.py, two padded cells) with the MATCHED express phase
phi_A = pi/(2n), psi_A = pi, psi_B = pi/2, against the certificate's dictionary:
  one scalar graph for both species, path  X_n ... X_0 D_0 ... D_n  (vertex 0 = X_n),
  diagonal a+D on X, a on D, path entries -1, closing entry H[X_n, D_n] = -i t_d,
  x = R |X_n, D_0>,  U_i = row D_i of x, V_j = column D_j of x,
  C = U^dag R V,  A = U^dag R U + V^dag R V,   T_c = -2 Re C,  mu_c = -2 i Im A      (units alpha^4)
which is what certificate/full_loop_integer.py encloses (F = 4 Re(...), M antisymmetrised).
Also: size of the lam^2 g^4 coefficient against the inherited bound alpha^4/4 (C84 PROOF (10)) at a = 4.
"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import second_order as so

def c84_form(n, a, D, td):
    m = n + 1; S = 2 * n + 2
    H = np.zeros((S, S), complex)
    for v in range(S):
        H[v, v] = a + D if v < m else a
    for v in range(S - 1):
        H[v, v + 1] = H[v + 1, v] = -1
    H[0, S - 1] += -1j * td; H[S - 1, 0] += 1j * td
    I = np.eye(S)
    R = np.linalg.inv(np.kron(H, I) + np.kron(I, H))
    src = np.zeros(S * S, complex); src[0 * S + m] = 1
    x = (R @ src).reshape(S, S)
    U, V = [], []
    for j in range(n):
        u = np.zeros((S, S), complex); u[m + j, :] = x[m + j, :]; U.append(u.reshape(-1))
        v = np.zeros((S, S), complex); v[:, m + j] = x[:, m + j]; V.append(v.reshape(-1))
    U = np.array(U).T; V = np.array(V).T
    C = U.conj().T @ R @ V
    A = U.conj().T @ R @ U + V.conj().T @ R @ V
    return -2 * C.real, -2j * A.imag

res = []
for (n, a, D, td) in [(2, 2.6, 3.0, 0.07), (2, 2.2, 0.4, 0.3), (3, 2.6, 3.0, 0.07), (3, 2.2, 0.4, 0.3), (2, 4.0, 3.0, 0.05), (3, 4.0, 3.0, 0.05)]:
    al4 = (2 / (n + 2)) ** 2 * np.sin(np.pi / (n + 2)) ** 4
    out, h, e, s = so.analyse(n, a, D, td, phiA=np.pi / (2 * n))
    T_bf, mu_bf = h[:n, n:] / al4, (h[:n, :n] + h[n:, n:]) / al4
    Tc, muc = c84_form(n, a, D, td)
    sc = abs(T_bf).max()
    # best global phase chi with T_bf = e^{i chi} T_c
    chi = np.angle(np.vdot(Tc.reshape(-1), T_bf.reshape(-1)))
    d = dict(n=n, a=a, D=D, td=td, phiA="pi/(2n)", gauss_resid_over_T=out["gauss_resid"] / (sc * al4),
             chi_over_pi=float(chi / np.pi),
             T_minus_phase_times_Tc_over_T=float(abs(T_bf - np.exp(1j * chi) * Tc).max() / sc),
             imag_part_of_gauged_T_over_T=float(abs((T_bf * np.exp(-1j * chi)).imag).max() / sc),
             real_part_of_mu_over_T=float(abs(mu_bf.real).max() / sc),
             mu_minus_muc_over_T=float(abs(mu_bf - muc).max() / sc),
             mu_plus_muc_over_T=float(abs(mu_bf + muc).max() / sc),
             mu_over_T=float(abs(mu_bf).max() / sc),
             lam2g4_coefficient_over_alpha4=float(out["M4_max"] / al4), inherited_bound=0.25)
    # spectra of the two Bloch matrices must agree after k -> k + chi (and possibly k -> -k)
    ks = np.linspace(-np.pi, np.pi, 41)
    def spec(T, mu, k): return np.linalg.eigvalsh(mu + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k))
    d["spectrum_mismatch_same_orientation_over_T"] = float(max(abs(spec(T_bf, mu_bf, k) - spec(Tc, muc, k + chi)).max() for k in ks) / sc)
    d["spectrum_mismatch_mirrored_over_T"] = float(max(abs(spec(T_bf, mu_bf, k) - spec(Tc, muc, -(k + chi))).max() for k in ks) / sc)
    res.append(d); print(json.dumps(d), flush=True)
json.dump(res, open(os.path.join(HERE, "item6_gauge.json"), "w"), indent=1)
