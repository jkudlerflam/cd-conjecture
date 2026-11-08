"""C164: check onset_dmrg.correlations against exact diagonalization on an open chain of 10 sites."""
import numpy as np, est, dmrg as D
from dmrg_c import chain_bond
import onset_dmrg as O
rng = np.random.default_rng(5); v = rng.normal(size=9); v /= np.linalg.norm(v)
h, edge = chain_bond(v); L = 10
r = D.dmrg(h, L, chis=(16, 32, 64, 64), edge_field=edge)
C, mag = O.correlations(r['M'], 2, 7)
# exact open chain
import scipy.sparse as sp, scipy.sparse.linalg as spla
def op(o, j):
    return sp.kron(sp.kron(sp.identity(2 ** j), sp.csr_matrix(o)), sp.identity(2 ** (L - j - 1)))
H = sum(sp.kron(sp.kron(sp.identity(2 ** j), sp.csr_matrix(h)), sp.identity(2 ** (L - j - 2))) for j in range(L - 1))
H = H + op(edge, 0) + op(edge, L - 1)
w, V = spla.eigsh(H.tocsr(), k=1, which='SA'); psi = V[:, 0]
ex = lambda A: np.vdot(psi, A @ psi).real
m = [[ex(op(est.PAULI[a], j)) for a in range(3)] for j in range(L)]
Ce = [sum(ex(op(est.PAULI[a], 2) @ op(est.PAULI[a], 2 + rr)) for a in range(3)) - np.dot(m[2], m[2 + rr]) for rr in range(8)]
print('E dmrg - ED', r['E'] - w[0]); print('max |C - C_ED|', np.max(np.abs(C - np.array(Ce)))); print('max |mag - ED|', np.max(np.abs(mag - np.array(m[2:10]))))
