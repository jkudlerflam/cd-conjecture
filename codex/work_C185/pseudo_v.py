"""C185: dilute two-valley criterion (C175/C180) for the Fibonacci-cusp model with a repulsive shift of V2,
V2 = -2 + v (units |t2| = 1), V3 = -2/3, at phi = pi/2. v = 16 lam^2/Delta^2 for a transverse DM lam (PROOF 3).
Uses C175's scat.py unmodified. Columns: v, a_s, A_d, |A_o|, Delta_sigma = (a_s - A_d)/|A_o|, verdict, and
min over K of E2 - 2 eps_min (negative = two-body bound state). Binding theorem (contrib): bound iff v < u^2/4."""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, "../work_C175")
from scat import A_matrix, a_like, valleys, rel_coeffs

def E0(w, K, V2, V3, Lr=1500):
    c = rel_coeffs(K, w); n = Lr; pot = np.zeros(n); pot[0] = V2; pot[1] = V3
    ab = np.zeros((3, n), complex); ab[0, :] = c[0]+pot; ab[1, :-1] = c[-1]; ab[2, :-2] = c[-2]
    return sl.eig_banded(ab, lower=True, select='i', select_range=(0, 0), eigvals_only=True)[0]

print("#  u     v     a_s      A_d     |A_o|   Dsig    verdict        E2(K=ka+kb)-2eps_min  thr u^2/4")
for u in (1.0, 2.0, 0.5, 1.5):
    w = 1j*u; em = valleys(w)[0][1]
    for v in (0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        V2, V3 = -2+v, -2/3
        A, ka, kb, K = A_matrix(w, V2=V2, V3=V3); AH = (A+A.conj().T)/2
        a_s = a_like(w, 0, V2=V2, V3=V3)[0].real; Ad = AH[0, 0].real; Ao = abs(AH[0, 1])
        D = (a_s-Ad)/Ao
        verd = 'c=2' if abs(D) < 1 else ('one valley' if D <= -1 else 'valley crystal')
        Lr = 3000 if u < 1 else 1500
        eb = E0(w, K, V2, V3, 1200) - 2*em   # inter-valley channel K = k_a + k_b only; binding boundary is proved (contrib)
        print(f"{u:4.1f} {v:5.2f} {a_s:8.3f} {Ad:8.3f} {Ao:8.3f} {D:7.3f}  {verd:14s} {eb:+.5f}   {u*u/4:.4f}", flush=True)
