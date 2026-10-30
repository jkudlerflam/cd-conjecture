"""Bonds (same conventions and random seeds as codex/work_C117/bootstrap.py) and the
gauge-optimized half-spread M."""
import numpy as np


def bonds():
    x = np.array([[0, 1], [1, 0.]])
    z = np.diag([1., -1.]); y = np.array([[0, -1j], [1j, 0]])
    xy = (np.kron(x, x) + np.kron(y, y)).real / 4
    out = [('Heisenberg', 2, xy + np.kron(z, z) / 4), ('XXZ_half', 2, xy + np.kron(z, z) / 8)]
    p = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            p[3 * i + j, 3 * j + i] = 1
    out.append(('SU3_swap', 3, p))
    rng = np.random.default_rng(117)
    for d in (2, 3):
        a = rng.integers(-8, 9, size=(d * d, d * d))
        out.append(('random_real_d' + str(d), d, (a + a.T) / 16))
    return out


# Known infrared data (A = sum c v), where available in closed form.
KNOWN_A = {
    'Heisenberg': np.pi / 2,                       # c=1, v=pi/2
    'XXZ_half': (np.pi / 2) * (np.sqrt(3) / 2) / (np.pi / 3),   # v = (pi/2) sin g/g, g=pi/3
    'SU3_swap': 2 * np.pi * 2 / 3,                 # c=2, v=2pi/3
}


def gauge_M(h, d):
    """min over real traceless o of (lmax-lmin)/2 of h + o x 1 - 1 x o."""
    import cvxpy as cp
    o = cp.Variable((d, d), symmetric=True)
    I = np.eye(d)
    K = h + cp.kron(o, I) - cp.kron(I, o)
    t1 = cp.Variable(); t2 = cp.Variable()
    n = d * d
    prob = cp.Problem(cp.Minimize((t1 - t2) / 2),
                      [t1 * np.eye(n) - K >> 0, K - t2 * np.eye(n) >> 0, cp.trace(o) == 0])
    prob.solve(solver='CLARABEL')
    return float(prob.value)
