"""Table of A, gauge-optimized M, the proved lower bound on M, and A/M for reference chains."""
import numpy as np
from ff_bound import casimir_A, trace_norm, bond_mj, nc_to_mj, half_spread
from gauge_M import gauge_M

rows = []

def add_nc(name, A0, B):
    A0 = np.array(A0, dtype=complex); B = np.array(B, dtype=complex)
    A0m, Bm = nc_to_mj(A0, B); n = A0m.shape[0]
    H = bond_mj(A0m, Bm, n); d = 2 ** (n // 2)
    M, _ = gauge_M(H, d)
    A = casimir_A("nc", A0, B)
    rows.append((name, d, A, trace_norm(B), half_spread(H), M, A / M))

def add_mj(name, A0, B):
    A0 = np.array(A0, float); B = np.array(B, float); n = A0.shape[0]
    H = bond_mj(A0, B, n); d = 2 ** (n // 2)
    M, _ = gauge_M(H, d)
    A = casimir_A("mj", A0, B)
    rows.append((name, d, A, trace_norm(B) / 2, half_spread(H), M, A / M))

add_nc("XX, half filling (t=1)", [[0]], [[-1]])
add_nc("XX, k_F = pi/3 (mu = -1)", [[-1.0]], [[-1]])
add_mj("critical Ising (J=g=1)", [[0, 1], [-1, 0]], [[0, 0], [1, 0]])
add_mj("critical Kitaev t=Delta=1, mu=2", [[0, 2], [-2, 0]], [[0, 0], [2, 0]])
add_mj("single Majorana chain x2 (B=1)", [[0, 0], [0, 0]], [[1, 0], [0, 1]])
add_nc("two XX chains (m=2)", np.zeros((2, 2)), -np.eye(2))
add_nc("two XX chains, on-site mixing 0.5", [[0, 0.5], [0.5, 0]], -np.eye(2))
add_nc("two XX, hoppings (1, 0.3), mixing 0.5", [[0, 0.5], [0.5, 0]], -np.diag([1, 0.3]))
add_nc("SSH critical (m=2, rank-one B)", [[0, 1], [1, 0]], [[0, 0], [1, 0]])
add_nc("XX, complex hopping e^{i pi/5}", [[0]], [[np.exp(1j * np.pi / 5)]])
rng = np.random.default_rng(5)
X = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)); A0 = (X + X.conj().T) / 2
B = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
add_nc("random m=2 chain (seed 5)", A0, B)

print(f"{'chain':42s} {'d':>2s} {'A':>8s} {'LB on M':>8s} {'bare M':>8s} {'M':>8s} {'A/M':>7s}")
for r in rows:
    print(f"{r[0]:42s} {r[1]:2d} {r[2]:8.4f} {r[3]:8.4f} {r[4]:8.4f} {r[5]:8.4f} {r[6]:7.4f}")
print("\nLB on M = ||B||_1 (number conserving) or ||B||_1/2 (Majorana form); bare M = half-spread with the on-site term split evenly.")
print("Proved: A <= pi * LB <= pi M.  Observed: A/M <= 2, equality for XX-type chains.")
