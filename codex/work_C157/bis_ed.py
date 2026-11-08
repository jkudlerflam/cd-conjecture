"""Ring ED of the bisector-field chain H = sum X_j Z_{j+1} + B sum (X_j + Z_j)/sqrt2 (Kitaev chain in a [110] field,
B = 2h/K). Prints the lowest levels, (E - E0) L, and E0/L, for L = 8..LMAX."""
import sys
import numpy as np
from xzlib import term, ring_operator
import scipy.sparse.linalg as spla
Bs = [float(x) for x in sys.argv[1].split(",")]; Ls = [int(x) for x in sys.argv[2].split(",")]
for B in Bs:
    h = term("xz") + (B / np.sqrt(2)) * (term("Bx") + term("Bz"))
    for L in Ls:
        op = ring_operator(h, L)
        w = np.sort(spla.eigsh(op, k=8, which="SA", tol=1e-10, return_eigenvectors=False, ncv=40).real)
        print(f"B = {B} L = {L}: E0/L = {w[0]/L:.8f}  (E-E0)L = {np.round((w[1:]-w[0])*L, 4).tolist()}", flush=True)
