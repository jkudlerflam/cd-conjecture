# XX chain as free fermions: eps(k) = -cos k (v = 1, c = 1, A = 1), energy current j(k) = eps*eps' = -cos k sin k.
# J commutes with H here, so d^2 G/da^2 at a = 0 is Var(J)/L = int dk/2pi j^2 f(1-f). IR prediction (pi/3) A T^3.
import numpy as np
k = np.linspace(-np.pi, np.pi, 2_000_001)[:-1]
eps = -np.cos(k); j = -np.cos(k)*np.sin(k)
for T in [0.005, 0.01, 0.02, 0.05]:
    f = 0.5*(1 - np.tanh(eps/(2*T)))
    var = np.mean(j**2 * f*(1-f))          # = int dk/2pi
    print(f"T={T:6.3f}  Var(J)/L={var:.6e}  (pi/3)T^3={np.pi/3*T**3:.6e}  ratio={var/(np.pi/3*T**3):.5f}")
