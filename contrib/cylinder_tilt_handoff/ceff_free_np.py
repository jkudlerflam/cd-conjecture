# C(r) = -(6R/pi)(E0(R) - e_inf R) with m = 1, r = mR, free massive fields (numpy only)
# boson periodic:  E0 - eR = (1/2pi) int dth cosh th log(1 - e^{-r cosh th}) = -(1/pi) sum_n K1(n r)/n
# Majorana NS:     E0 - eR = -(1/2pi) int dth cosh th log(1 + e^{-r cosh th})
# Majorana R:      E0 - eR = -(boson periodic), mode by mode (E0 = -1/2 sum omega vs +1/2 sum omega, same momenta)
import numpy as np
def C(r, kind):
    tmax = np.log(2*60/r) + 2
    th = np.linspace(-tmax, tmax, 400001); x = r*np.cosh(th)
    if kind == 'B':  f = -np.log1p(-np.exp(-x))
    if kind == 'NS': f = np.log1p(np.exp(-x))
    I = np.trapezoid(np.cosh(th)*f, th)
    return 3*r/np.pi**2 * I
for r in [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 5, 10, 20]:
    b = C(r, 'B'); ns = C(r, 'NS')
    print(f"r={r:8.4g}  boson C={b: .6f}  Majorana NS C={ns: .6f}  Majorana R C={-b: .6f}")
