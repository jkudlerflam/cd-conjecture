# Check the three conditions on C(R) (m = 1) for free fields:
#  s >= 0      <=> 2C - R C' >= 0
#  C_V >= 0    <=> (C/R)'' >= 0
#  C_V >= s    <=> d^2 C / d q^2 >= 0, q = R^2   (real-tilt convexity)
import numpy as np
from scipy.special import kv
def C(R, kind, K=4000):
    k = np.arange(1, K+1)[:, None]
    sgn = (-1.0)**(k+1) if kind == 'NS' else 1.0
    S = np.sum(sgn*kv(1, k*R)/k, axis=0)
    c = 6*R/np.pi**2*S
    return -c if kind == 'Ramond' else c
R = np.exp(np.linspace(np.log(0.02), np.log(8), 4001)); q = R**2
for kind in ['boson', 'NS', 'Ramond']:
    c = C(R, kind)
    dc = np.gradient(c, R)
    s_cond = 2*c - R*dc
    cv_cond = np.gradient(np.gradient(c/R, R), R)
    tilt = np.gradient(np.gradient(c, q), q)
    sl = slice(5, -5)
    print(f"{kind:7s} min(2C-RC')={s_cond[sl].min(): .3e}  min((C/R)'')={cv_cond[sl].min(): .3e}  "
          f"min(C_qq)={tilt[sl].min(): .3e}  max(C_qq)={tilt[sl].max(): .3e}")
