"""Screw-symmetric two-site chains are uniform after the gauge R_j = P^j: check Kitaev-XX-Gamma and the angle-theta compass chain."""
import numpy as np
from ed import build
def comp(Jo, Je, Bo, Be, Ju, Bu, L=8):
    eo = np.linalg.eigvalsh(build(lambda j: Jo if j % 2 == 0 else Je, lambda j: Bo if j % 2 == 0 else Be, L).toarray())
    eu = np.linalg.eigvalsh(build(lambda j: Ju, lambda j: Bu, L).toarray())
    return np.max(np.abs(eo - eu))
# Kitaev-XX-Gamma (arXiv:2509.21901): x-bond K xx + J(xx+yy) + G(xy+yx), y-bond K yy + same; plus a field along (1,1,0)
K, J, G, hh = 1.0, 0.37, -0.61, 0.23
base = J * np.diag([1, 1, 0.]) + G * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.]])
Jx = base + K * np.diag([1, 0, 0.]); Jy = base + K * np.diag([0, 1, 0.])
P = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1.]])     # pi rotation about (1,1,0)
B = hh * np.array([1, 1, 0.])
print("Kitaev-XX-Gamma (+[110] field): max diff", comp(Jx, Jy, B, B, Jx @ P, B))
# compass chain with angle theta: bonds n_+ n_+^T and n_- n_-^T, n_pm = (sin th, 0, +-cos th)... in xz plane
th = 0.4
n1 = np.array([np.sin(th), 0, np.cos(th)]); n2 = np.array([-np.sin(th), 0, np.cos(th)])
R = np.array([[np.cos(2*th), 0, -np.sin(2*th)], [0, 1, 0], [np.sin(2*th), 0, np.cos(2*th)]])  # rotation about y
Q = R if np.allclose(R @ n1, n2) else R.T
assert np.allclose(Q @ n1, n2)
Q2 = Q @ Q
# need P with P n1 = n2 and P n2 = n1 and P in SO(3): pi rotation about z-axis? use reflection-free: pi rotation about (0,0,1)? maps (s,0,c)->(-s,0,c). yes
Pc = np.diag([-1, -1, 1.])
assert np.allclose(Pc @ n1, n2) and np.allclose(Pc @ n2, n1)
J1 = np.outer(n1, n1); J2 = np.outer(n2, n2)
print("angle compass chain: max diff", comp(J1, J2, np.zeros(3), np.zeros(3), J1 @ Pc, np.zeros(3)))
