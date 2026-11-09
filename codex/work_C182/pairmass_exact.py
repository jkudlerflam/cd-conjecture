"""C182: exact ab-pair inverse mass at phi = pi/2 from the collaborator's closed form
(contrib/cusp_scattering_20261109/DERIVATION.md eqs 1, 5-8), evaluated for u = |w| up to 4, and checked against
the banded numerics of pairband.py at a few u."""
import numpy as np
from scipy.optimize import brentq
from pairband import lowest, Eb_exact
def inv_mass(u):
    p = brentq(lambda p: (1-p)*(3-p)*(1+p)**2 - u*u*p**3, 1e-14, 1-1e-14)
    s = -u*p/(1+p); z = (s+1j*np.sqrt(4*p-s*s))/2
    C = np.conj(1+u*z+z*z)
    Nn = 0.5*(abs(C)**2*p**2/(1-p) - np.real(C**2*z**4/(1-z**2)))
    Cm = lambda m: (abs(C)**2*np.real(z**m)*p**2/(1-p) - np.real(C**2*z**(4+m)/(1-z**2)))/Nn
    return -u/2*Cm(1) - 2*Cm(2), p
print("# u  1/m_pair(exact)  free=(8-u^2/2)/2  ratio  p  E_b")
for u in [0.5,1,1.5,2,2.5,3,3.25,3.5,3.6,3.7,3.8,3.9,3.95,3.99,3.999]:
    im, p = inv_mass(u)
    fr = (8-u*u/2)/2
    print(f"{u:6.3f} {im:+.6f} {fr:.6f} {im/fr:.4f} {p:.6f} {Eb_exact(u):.6f}")
# numerical check of the curvature at u = 3.5 and 3.9 (R large: pair size grows near u = 4)
for u, R in ((3.5, 3000), (3.9, 6000)):
    h = 2e-3; f = lambda K: lowest(K, 1j*u, R, k=1)[0]
    print("numeric curvature u=%g R=%d: %.5f" % (u, R, (f(np.pi+h)-2*f(np.pi)+f(np.pi-h))/h**2))
