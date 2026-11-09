#!/usr/bin/env python3
"""C170: is a uniform field (the relevant momentum-0, charge-0 operator J^z + Jbar^z of the mu-frame chain) generated
at second order by U(1)-breaking harmonics at hz = 0?  Ring L = 12, a = 2 pi 4/12 (aL in 2 pi Z), Dz = 0.9 (liquid side).
Measure m = <sum sigma^z>/L in the ground state (sigma^z = mu^z)."""
import numpy as np, scipy.sparse.linalg as sla
from qed import ham, op, P
from lr2 import parent, harm
L = 12; a = 2 * np.pi * 4 / 12
Sz = sum(op(P[2], j, L) for j in range(L))
for kw in (dict(), dict(hx=0.15), dict(G=0.15), dict(hx=0.15, G=0.15), dict(hx=0.075, G=0.075), dict(hx=0.15, Dx=0.15), dict(g=0.15, hx=0.15)):
    dJ, dB = harm(**kw)
    H = ham(parent(a, 0.9) + dJ, dB, L, pbc=True)[0]
    w, v = sla.eigsh(H, k=2, which='SA', tol=1e-12)
    i = np.argmin(w); psi = v[:, i]
    m = np.real(psi.conj() @ (Sz @ psi)) / L
    print(f'{str(kw):32s} E0/L = {w[i]/L:.6f}  gap = {abs(w[1]-w[0]):.4f}  m = {m:+.5f}')
