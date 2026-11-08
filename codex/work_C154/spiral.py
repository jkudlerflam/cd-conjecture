"""C154: canonical spiral frame of a qubit bond and its U(1) parent (Claude Opus agent, 2026-11-07).
J = S O (O in SO(3), S symmetric) by polar decomposition; O = rotation by theta about m.
In the frame sigma_j = O^{-j} tau_j the bond is O^j S O^-j and the field O^j B.  The parent is the
rotation-about-m average: S0 = s_m m m^T + s_perp (1 - m m^T), s_m = m.S.m, s_perp = (tr S - s_m)/2, B0 = (m.B) m.
J(lam) = (S0 + lam (S - S0)) O, B(lam) = B0 + lam (B - B0) joins the U(1)-class parent (lam = 0) to the bond (lam = 1).
Functions: polar(J), parent(J, B), path(J, B, lam), u1_check(J, B), harmonic content."""
import numpy as np

def rot_axis_angle(O):
    th = np.arccos(np.clip((np.trace(O) - 1) / 2, -1, 1))
    if th < 1e-12:
        return np.array([0, 0, 1.0]), 0.0
    if abs(th - np.pi) < 1e-9:
        w, V = np.linalg.eigh((O + O.T) / 2); m = V[:, np.argmax(w)]
        return m / np.linalg.norm(m), np.pi
    m = np.array([O[2, 1] - O[1, 2], O[0, 2] - O[2, 0], O[1, 0] - O[0, 1]]) / (2 * np.sin(th))
    return m / np.linalg.norm(m), th

def polar(J):
    U, s, Vt = np.linalg.svd(J)
    O = U @ Vt; S = U @ np.diag(s) @ U.T
    if np.linalg.det(O) < 0:
        O, S = -O, -S
    return S, O

def parent(J, B):
    S, O = polar(J)
    m, th = rot_axis_angle(O)
    if th < 1e-12:   # symmetric J: no spiral, choose m as the S-eigenaxis closest to B (no momentum protection)
        w, V = np.linalg.eigh(S); m = V[:, np.argmax(np.abs(V.T @ B))] if np.linalg.norm(B) > 0 else V[:, 0]
    sm = m @ S @ m; sp = (np.trace(S) - sm) / 2
    P = np.outer(m, m)
    S0 = sm * P + sp * (np.eye(3) - P)
    B0 = (m @ B) * m
    return dict(S=S, O=O, m=m, theta=th, S0=S0, B0=B0, s_m=sm, s_perp=sp)

def path(J, B, lam):
    p = parent(J, B)
    return (p['S0'] + lam * (p['S'] - p['S0'])) @ p['O'], p['B0'] + lam * (np.asarray(B) - p['B0'])

def u1_check(J, B, m):
    """max deviation of R^T J R - J and R^T B - B for rotations about m (U(1) about m)."""
    def R(a):
        K = np.array([[0, -m[2], m[1]], [m[2], 0, -m[0]], [-m[1], m[0], 0]])
        return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K
    return max(max(np.abs(R(a).T @ J @ R(a) - J).max(), np.abs(R(a).T @ B - B).max()) for a in (0.3, 1.1, 2.0))

if __name__ == "__main__":
    import json
    d = json.load(open('../work_C68/checks.json'))['point_1282']
    J = np.array(d['J']); B = np.array(d['B'])
    p = parent(J, B)
    print("point 1282: theta/pi =", p['theta'] / np.pi, " m =", np.round(p['m'], 4), " s_m =", p['s_m'], " s_perp =", p['s_perp'])
    print("  S eigenvalues", np.linalg.eigvalsh(p['S']), " |S - S0| / |S| =", np.linalg.norm(p['S'] - p['S0']) / np.linalg.norm(p['S']),
          " |B - B0| / |B| =", np.linalg.norm(B - p['B0']) / np.linalg.norm(B))
    J0, B0 = path(J, B, 0.0)
    print("  parent U(1) deviation:", u1_check(J0, B0, p['m']), " reconstruct lam=1:", np.abs(np.array(path(J, B, 1.0)[0]) - J).max())
    # parent in its own frame: XXZ with J_perp = |s_perp| (after removing the DM by the spiral), Delta = s_m/|s_perp|, field m.B
    print("  parent XXZ: J_perp_eff =", abs(p['s_perp']), " Delta =", p['s_m'] / abs(p['s_perp']), " h =", p['m'] @ B,
          " h/J_perp =", (p['m'] @ B) / abs(p['s_perp']))
    # test on the cyclic stratum: XZ + 0.1 YX
    Jc = np.zeros((3, 3)); Jc[0, 2] = 1; Jc[1, 0] = 0.1
    pc = parent(Jc, np.zeros(3))
    print("XZ+0.1YX: theta/pi =", pc['theta'] / np.pi, " m =", np.round(pc['m'], 4), " S =", np.round(pc['S'], 4).tolist())
