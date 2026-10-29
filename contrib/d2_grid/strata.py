"""Residual strata of the d = 2 family and random great 2-spheres through them.

Convention (as in codex/work_C6_d2 and contrib/d2_classification/ed.py):
H = sum_j [s_j^T J s_{j+1} + f . s_j], J = diag(Jx, Jy, Jz) + A(D), A_ab = eps_abc D_c,
so D_z multiplies X_j Y_{j+1} - Y_j X_{j+1}.  A point is the 9-vector
v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz) in the principal frame of the symmetric part.
Positive rescaling is irrelevant, so each stratum is a sphere in its coordinate subspace.

Strata (C6 classification, principal frame, generic points):
  NS  no uniform internal symmetry: all nine coordinates            (dimension 8)
  Z2  one unitary Z2 about z: (Jx, Jy, Jz, Dz, fz)                  (dimension 4)
  TR  time reversal only: f = 0, (Jx, Jy, Jz, Dx, Dy, Dz)           (dimension 5)
  RA  plane-reflection antiunitary: D || z, f in the xy plane:
      (Jx, Jy, Jz, Dz, fx, fy)                                      (dimension 5)
The staggered-field problem a(XX - YY) + u ZZ + b Z (a = 1) is a separate 2D grid.

Written for contrib/d2_grid (Claude Opus 5.5 subagent, 2026-10-28).
"""
import numpy as np

COORDS = {
    'NS': [0, 1, 2, 3, 4, 5, 6, 7, 8],
    'Z2': [0, 1, 2, 5, 8],
    'TR': [0, 1, 2, 3, 4, 5],
    'RA': [0, 1, 2, 5, 6, 7],
}


def J_f(v):
    v = np.asarray(v, float)
    J = np.diag(v[:3]).astype(float)
    Dx, Dy, Dz = v[3:6]
    # A_ab = eps_abc D_c
    J[0, 1] += Dz; J[1, 0] -= Dz
    J[1, 2] += Dx; J[2, 1] -= Dx
    J[2, 0] += Dy; J[0, 2] -= Dy
    return J, v[6:9].copy()


def sphere_basis(stratum, seed):
    """Three orthonormal vectors in the coordinate subspace of the stratum."""
    rng = np.random.default_rng(seed)
    idx = COORDS[stratum]
    M = rng.normal(size=(len(idx), 3))
    Q, _ = np.linalg.qr(M)
    U = np.zeros((3, 9))
    U[:, idx] = Q.T
    return U


def fibonacci(N):
    """N nearly uniform points on the unit 2-sphere."""
    k = np.arange(N) + 0.5
    z = 1 - 2 * k / N
    phi = np.pi * (1 + 5 ** 0.5) * k
    r = np.sqrt(1 - z ** 2)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], 1)


def sphere_points(stratum, seed, N):
    U = sphere_basis(stratum, seed)
    X = fibonacci(N)
    return X, X @ U


def staggered(u, b, a=1.0):
    """a(XX - YY) + u ZZ + b Z in the J, f convention."""
    return np.array([a, -a, u, 0, 0, 0, 0, 0, b], float)


def symmetry_report(v, tol=1e-9):
    """Uniform internal symmetries of a principal-frame point (C6 eqs. Gu, Ga), checked by brute force
    over the principal-axis pi rotations and the improper involutions that preserve diag S."""
    J, f = J_f(v)
    S = (J + J.T) / 2
    D = np.array([J[1, 2] - J[2, 1], J[2, 0] - J[0, 2], J[0, 1] - J[1, 0]]) / 2
    out = []
    for name, R in [('Rx', np.diag([1, -1, -1])), ('Ry', np.diag([-1, 1, -1])), ('Rz', np.diag([-1, -1, 1]))]:
        if np.allclose(R @ S @ R.T, S, atol=tol) and np.allclose(R @ D, D, atol=tol) and np.allclose(R @ f, f, atol=tol):
            out.append(name)
    for name, O in [('T', -np.eye(3)), ('Kx', np.diag([-1, 1, 1])), ('Ky', np.diag([1, -1, 1])), ('Kz', np.diag([1, 1, -1]))]:
        if np.allclose(O @ S @ O.T, S, atol=tol) and np.allclose(O @ D, -D, atol=tol) and np.allclose(O @ f, f, atol=tol):
            out.append(name)
    return out
