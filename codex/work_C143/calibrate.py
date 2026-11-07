"""Calibrate the search objective g = log S_L2 - log S_L1 on known regular and singular bonds."""
import time, sys
import numpy as np
from search import objective, softcount

def X(J, B=(0, 0, 0)):
    return np.concatenate([np.asarray(J, float).ravel(), np.asarray(B, float)])

def proj_x(psi):
    # Pauli coefficients of |psi><psi| : h = sum c_ab s^a s^b + fields; return (J,B) with gauge part dropped
    from ed import PAULI, I2
    psi = np.asarray(psi, complex); psi /= np.linalg.norm(psi)
    P = np.outer(psi, psi.conj())
    J = np.array([[np.trace(P @ np.kron(PAULI[a], PAULI[b])).real / 4 for b in range(3)] for a in range(3)])
    l = np.array([np.trace(P @ np.kron(PAULI[a], I2)).real / 4 for a in range(3)])
    r = np.array([np.trace(P @ np.kron(I2, PAULI[a])).real / 4 for a in range(3)])
    return X(J, l + r)   # uniform field part: (l+r)/2 per term, times 2 because B enters as B(s1+s2)/2

cases = {
 "Fibonacci ZZ+(Z1+Z2)": X(np.diag([0, 0, 1]), (0, 0, 2)),
 "XZ point": X([[0, 0, 1], [0, 0, 0], [0, 0, 0]]),
 "FM Heisenberg": X(-np.eye(3)),
 "AF Heisenberg": X(np.eye(3)),
 "XX": X(np.diag([1, 1, 0])),
 "TFIM critical": X(np.diag([-1, 0, 0]), (0, 0, -1)),
 "XX saturation (gamma=0,h=1)": X(np.diag([0.5, 0.5, 0]), (0, 0, -1)),
 "XXZ D=2 at h=0 (gapped Neel)": X(np.diag([1, 1, 2])),
 "Ising FM (gapped, 2-fold)": X(np.diag([0, 0, -1])),
 "BG real psi p=0.2": proj_x([np.sqrt(0.8), 0, 0, np.sqrt(0.2)]),
}
L1, L2 = int(sys.argv[1]), int(sys.argv[2])
for nm, x in cases.items():
    t = time.time()
    g = objective(x / np.linalg.norm(x), L1, L2)
    print(f"{nm:34s} S_L1={softcount(x/np.linalg.norm(x), L1):9.3f} S_L2={softcount(x/np.linalg.norm(x), L2):9.3f}  g={g:7.3f}  ({time.time()-t:.2f}s)")
