"""Extrapolate the ED columns of xxz_ed_check.out in 1/L^2 (fit x(L) = x_inf + a/L^2 + b/L^4 over
the four largest sizes) and compare with the closed-form targets.  Also evaluates the
path-integrated distinguishability  I = int_{-1}^{1} D_2(Delta) dDelta  (finite, since
D_2 ~ (pi sqrt2/6) (1+Delta)^{-1/2} near the ferromagnetic point), and the identity
A(Delta) = A(-1) + int_{-1}^{Delta} A' dDelta' with the trace-norm majorant.
"""
import re, numpy as np
from scipy.integrate import quad
from xxz_closed_form import D2_closed, v_of, vprime_of

txt = open("xxz_ed_check.out").read()
blocks = txt.split("== XXZ Delta = ")[1:]
print("Delta    L^2 dg_zz(extrap)  target (pi/6)A'   L^2 d(2g_xx)(extrap)  target   D_2(extrap)  D_2 closed   A(extrap)  A exact")
for blk in blocks:
    Delta = float(blk.split(":")[0])
    rows = [list(map(float, l.split())) for l in blk.splitlines() if re.match(r"\s+\d+\s+-", l)]
    rows = np.array(rows)
    L = rows[:,0]; Acol = rows[:,2]; gzz = rows[:,3]; gxx = rows[:,4]; D2 = rows[:,5]
    X = np.vstack([np.ones_like(L), L**-2., L**-4.]).T
    def ex(y): return np.linalg.lstsq(X[-4:], y[-4:], rcond=None)[0][0]
    D2c, a, b = D2_closed(Delta)
    print(f"{Delta:6.2f}  {ex(gzz):12.5f}      {a:10.5f}     {ex(gxx):12.5f}      {b:8.5f}   {ex(D2):10.5f}  {D2c:10.5f}   {ex(Acol):9.5f}  {v_of(Delta):8.5f}")

print()
I, err = quad(lambda D: D2_closed(D)[0], -1, 1, limit=400)
print(f"int_{{-1}}^{{1}} D_2(Delta) dDelta = {I:.5f}   (times ||S^zS^z|| = 1/4: {I/4:.5f})")
print(f"check: A(1) - A(-1) = {v_of(0.999999) - v_of(-0.999999):.5f} = int A' dDelta;  pi/2 = {np.pi/2:.5f}")
print("path majorant: dA/dDelta = (6/pi) lim L^2 Tr[S^zS^z X] <= (6/pi) ||S^zS^z|| D_2 = (3/(2 pi)) D_2, so")
print(f"  A(1) - A(-1) = pi/2 = {np.pi/2:.5f}  <=  (3/(2 pi)) int D_2 dDelta = {3/(2*np.pi)*I:.5f}")
