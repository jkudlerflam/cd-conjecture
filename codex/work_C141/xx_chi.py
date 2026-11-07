"""Test of the per-sector CFT weights of the energy density (PROOF.md section 1.2).
XX chain H = sum (SxSx + SySy) = (1/2) sum (c^dag_x c_{x+1} + h.c.), eps(k) = cos k,
half filling, c = 1, v = 1, A = 1.  Bond energy eps_x = (1/2)(c^dag_x c_{x+1} + h.c.).
Predicted small-q forms per unit length:
  S_eps(q)  = c v^2 |q|^3 /(24 pi)
  m1(q)     = c v^3 q^4   /(24 pi)
  chi_eps(q)= 2 sum |<n|eps_q|0>|^2/omega_n /L = c v q^2/(12 pi)
"""
import numpy as np
for L in [4000, 16000]:
    k = 2*np.pi*(np.arange(L)+0.5)/L   # antiperiodic grid avoids zero modes
    occ = np.cos(k) < 0
    print(f"L={L}")
    for m in [1, 2, 4, 8, 16]:
        q = 2*np.pi*m/L
        kq = k + q
        # transition c^dag_k c_{k+q}: needs k+q occupied, k empty
        sel = (np.cos(kq) < 0) & (~occ)
        w = np.cos(k[sel]) - np.cos(kq[sel])
        me2 = np.cos(k[sel] + q/2)**2
        S = me2.sum()/L
        m1 = (w*me2).sum()/L
        chi = 2*(me2/w).sum()/L
        print(f"  q={q:.5f}  S/(q^3/24pi)={S/(q**3/(24*np.pi)):.5f}  "
              f"m1/(q^4/24pi)={m1/(q**4/(24*np.pi)):.5f}  chi/(q^2/12pi)={chi/(q**2/(12*np.pi)):.5f}")
