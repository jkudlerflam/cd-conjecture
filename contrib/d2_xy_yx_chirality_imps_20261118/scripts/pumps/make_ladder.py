"""Decoupled two-leg calibration ladders for pumps.py (model=file), d = 4, site = (leg a) x (leg b).
h2 = ha(a_j, a_{j+1}) x 1_b + 1_a x hb(b_j, b_{j+1}), index order [(s_j, s_{j+1}), (t_j, t_{j+1})], s = 2 a + b."""
import numpy as np
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1.]).astype(complex)
def xxz(Jz): return np.kron(X, X) + np.kron(Y, Y) + Jz * np.kron(Z, Z)
def ladder(ha, hb):
    A = ha.reshape(2, 2, 2, 2); B = hb.reshape(2, 2, 2, 2); I = np.eye(2)
    H = (np.einsum('pqrs,bB,cC->pbqcrBsC', A, I, I) + np.einsum('pqrs,aA,dD->apdqArDs', B, I, I))
    return H.reshape(16, 16)
for name, (ja, jb) in {'xx_xx': (0.0, 0.0), 'xx_xxz05': (0.0, 0.5)}.items():
    h2 = ladder(xxz(ja), xxz(jb)); assert np.allclose(h2, h2.conj().T); np.save(name + '.npy', h2)
    print(name, h2.shape)
