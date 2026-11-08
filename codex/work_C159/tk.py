"""C159 library for the strata T and K (Claude Opus agent, 2026-11-07).

Convention (as C155 est.py / contrib/d2_grid): v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz),
H = sum_j [s_j^T J s_{j+1} + f . s_j], J = diag(Jx, Jy, Jz) + A(D), A_ab = eps_abc D_c (Pauli matrices).

T: f = 0 (coordinates Jx, Jy, Jz, Dx, Dy, Dz).
K: D = (0, 0, Dz), f = (fx, fy, 0); the antiunitary is Theta R_pi(z) (reflection through the xy plane).
"""
import numpy as np
from scipy.linalg import polar

from est import J_f, bond_from_JB, PAULI


def norm(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def rand_T(rng):
    v = np.zeros(9); v[:6] = rng.normal(size=6)
    return norm(v)


def rand_K(rng):
    v = np.zeros(9); v[[0, 1, 2, 5, 6, 7]] = rng.normal(size=6)
    return norm(v)


def rot_axis_angle(O):
    th = np.arccos(np.clip((np.trace(O) - 1) / 2, -1, 1))
    if th < 1e-10:
        return np.array([0, 0, 1.]), 0.0
    if np.pi - th < 1e-7:
        w, V = np.linalg.eigh((O + O.T) / 2)
        return V[:, np.argmax(w)], th
    m = np.array([O[2, 1] - O[1, 2], O[0, 2] - O[2, 0], O[1, 0] - O[0, 1]]) / (2 * np.sin(th))
    return m / np.linalg.norm(m), th


def spiral_parent(J, f=np.zeros(3)):
    """Polar decomposition J = S O (O in SO(3)), O = rotation by theta about m.
    Parent (C154 section 2): XXZ about m with J_perp = s_perp, J_par = s_m, field m.f.
    Returns dict with Delta = s_m / s_perp, theta, m, harmonic fractions."""
    U, P = polar(J, side='left')          # J = P U, P symmetric psd, U orthogonal
    if np.linalg.det(U) < 0:
        U, P = -U, -P
    S, O = P, U
    m, th = rot_axis_angle(O)
    s_m = m @ S @ m
    s_perp = (np.trace(S) - s_m) / 2
    S0 = s_m * np.outer(m, m) + s_perp * (np.eye(3) - np.outer(m, m))
    f0 = (m @ f) * m
    # Delta in the antiferromagnetic convention: a staggered pi rotation about m flips the sign of J_perp,
    # so the parent (s_perp, s_m) is equivalent to (|s_perp|, s_m); Delta = +1 AF Heisenberg, -1 ferromagnet.
    return dict(S=S, O=O, m=m, theta=th, s_m=s_m, s_perp=s_perp, Delta=s_m / abs(s_perp) if abs(s_perp) > 1e-12 else np.inf,
                Delta_raw=s_m / s_perp if abs(s_perp) > 1e-12 else np.inf,
                harm=np.linalg.norm(S - S0) / np.linalg.norm(S), fharm=np.linalg.norm(f - f0) / max(np.linalg.norm(f), 1e-14),
                hm=m @ f)


def luttinger_K(Delta):
    """Zero-field XXZ (J_perp > 0 convention, Delta = J_par / J_perp), |Delta| < 1."""
    if not (-1 < Delta < 1):
        return np.nan
    g = np.arccos(Delta)
    return np.pi / (2 * (np.pi - g))


def symmetry_check(J, f, tol=1e-9):
    """Uniform on-site symmetries among the 8 principal-axis sign actions of the symmetric part
    (exhaustive when S has distinct eigenvalues, C68), plus -1 (Kramers Theta).
    Returns list of (Q, det) with Q^T J Q = J, Q^T f' = f where for antiunitary Q (det -1) the condition
    on f is Q f = f with the Pauli action -R (C68 convention: improper Q)."""
    Ssym = (J + J.T) / 2
    w, V = np.linalg.eigh(Ssym)
    out = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                Q = V @ np.diag([s1, s2, s3]) @ V.T
                if np.allclose(Q, np.eye(3)):
                    continue
                if np.allclose(Q.T @ J @ Q, J, atol=tol) and np.allclose(Q.T @ f, f, atol=tol):
                    out.append((np.round(Q, 6), int(round(np.linalg.det(Q)))))
    return out


def corr_matrix(psi, L, r):
    """<s^a_0 s^b_r> for the ring ground state psi (computational basis, site 0 most significant)."""
    C = np.zeros((3, 3), complex)
    t = psi.reshape([2] * L)
    for a in range(3):
        ta = np.tensordot(PAULI[a], t, axes=([1], [0]))
        for b in range(3):
            tb = np.moveaxis(np.tensordot(PAULI[b], np.moveaxis(ta, r, 0), axes=([1], [0])), 0, r)
            C[a, b] = np.vdot(t.reshape(-1), tb.reshape(-1))
    return C


def magnetization(psi, L):
    t = psi.reshape([2] * L)
    return np.array([np.vdot(t.reshape(-1), np.tensordot(PAULI[a], t, axes=([1], [0])).reshape(-1)).real for a in range(3)])
