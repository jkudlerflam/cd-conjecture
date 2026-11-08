"""C162. 1D zero-energy scattering lengths (even channel, hard core at R = 0 is irrelevant when the potential is
strongly repulsive) for psi'' = k2 * V(R) psi, k2 = M / hbar^2, normalized psi -> R - a at large R.
Checks: (1) single exponential V = g e^{-lam R}: a = (2/lam)[ln(sqrt(k2 g)/lam) + gamma_E] (K_0 solution), so
a_ij for V_ij = g q_i q_j e^{-lam R} is additive (hard rods) and a_11 + a_22 - a_12 - a_21 = 0.
(2) two channels V_ij = g1 r_i l_j e^{-lam1 R} + g2 x_i y_j e^{-lam2 R}, lam2 > lam1:
J := a_11 + a_22 - a_12 - a_21 has the sign of g2 (x_1 - x_2)(y_1 - y_2) as k2 -> infinity.
(3) monotonicity da/dV >= 0 (Wronskian identity da/dlam = k2 int psi^2 dV/dlam), checked by finite differences."""
import numpy as np
from scipy.integrate import solve_ivp

def scatlen(V, k2, R1=None):
    """Riccati form u = psi'/psi, u' = k2 V - u^2, started on the WKB branch u = sqrt(k2 V) where k2 V = 400
    (deep enough under the barrier that the start is irrelevant), integrated out to R1, then a = R1 - 1/u."""
    from scipy.optimize import brentq
    R0 = brentq(lambda R: k2 * V(R) - 400.0, -20, 60)
    R1 = R0 + 60.0 if R1 is None else R1
    sol = solve_ivp(lambda R, u: [k2 * V(R) - u[0] ** 2], (R0, R1), [20.0], rtol=1e-11, atol=1e-13, method='LSODA')
    return R1 - 1.0 / sol.y[0, -1]

lam1, lam2 = 1.0, 2.5
print('(1) single exponential, g = 1, lam = 1')
for k2 in (1e2, 1e4, 1e6):
    a = scatlen(lambda R: np.exp(-lam1 * R), k2)
    print(f'   k2={k2:.0e}: a = {a:.6f}, K_0 formula = {2*(np.log(np.sqrt(k2))+np.euler_gamma):.6f}', flush=True)
print('(2) two channels plus a common core e^{-5R}; valley amplitudes r, l (slow, repulsive) and x, y (fast), J = a11 + a22 - a12 - a21')
# physical (Manton) convention for a decoupled massive scalar channel: V_ij = -c x_i y_j e^{-lam R}, c > 0,
# x_i = right-tail amplitude of the left particle (valley i), y_j = left-tail amplitude of the right particle (valley j).
# prediction: sign J = -sign((x_1 - x_2)(y_1 - y_2)); J < 0 is immiscible (Ising-ferro pseudospin).
c2 = 1.0
cases = [('K-like: fast tilt channel, x = +-B_R, y = +-B_L, B > 0', (1, 1), (1, 1), (0.6, -0.6), (0.4, -0.4)),
         ('U-like: x_- = -y_+, y_- = -x_+ in both channels', (1, 0.8), (0.8, 1), (0.6, -0.3), (0.3, -0.6)),
         ('no symmetry, self-repulsive fast channel (x_i y_i < 0), valley-dependent', (1, 1), (1, 1), (0.6, 0.2), (-0.5, -0.1))]
for name, r, l, x, y in cases:
    for k2 in (1e3, 1e5):
        a = {}
        for i in (0, 1):
            for j in (0, 1):
                V = lambda R, i=i, j=j: abs(r[i] * l[j]) * np.exp(-lam1 * R) - c2 * x[i] * y[j] * np.exp(-lam2 * R) + np.exp(-5.0 * R)
                a[i, j] = scatlen(V, k2)
        J = a[0, 0] + a[1, 1] - a[0, 1] - a[1, 0]
        pred = -np.sign((x[0] - x[1]) * (y[0] - y[1]))
        print(f'   {name}: k2={k2:.0e}  J={J:+.3e}  predicted sign {pred:+.0f}', flush=True)
print('(3) monotonicity: V_lam = e^{-R} + lam * 0.3 e^{-2.5 R}, k2 = 1e4')
k2 = 1e4
for lam in (0.0, 0.5, 1.0):
    V = lambda R: np.exp(-R) + lam * 0.3 * np.exp(-2.5 * R)
    Vp = lambda R: np.exp(-R) + (lam + 1e-3) * 0.3 * np.exp(-2.5 * R)
    print(f'   lam={lam}: a={scatlen(V, k2):.6f}, da/dlam (FD) = {(scatlen(Vp, k2) - scatlen(V, k2))/1e-3:+.4e}')
