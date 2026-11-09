"""C181: explicit qubit bond for effective-model parameters (C173 PROOF 0.1 and 3.1).
w = J_x + J_y + 2iD, Delta = J_x - J_y, V = 4 J_z, mu = 4 J_z + 2 h_z = |t2| (mu_eff + 6), |t2| = Delta^2/(2V).
Bond h = J_x XX + J_y YY + D (XY - YX) + J_z ZZ + (h_z/2)(Z(x)1 + 1(x)Z).
Also the size of the first number-changing vertices relative to |t2| (C179 claim 1, collaborator's cusp_validation):
separation-2 pair creation 2 Delta Re(w)/V, separation-3 pair creation 2 Delta^3/(3 V^2).
Usage: python3 bond.py wa phi mu_eff V_over_Delta"""
import sys, numpy as np
wa, phi, mu_eff, r = (float(x) for x in sys.argv[1:5])
Delta = 1.0; V = r * Delta; t2 = Delta ** 2 / (2 * V)
w = wa * t2 * np.exp(1j * phi)
s = w.real; D = w.imag / 2
Jx = (Delta + s) / 2; Jy = (s - Delta) / 2; Jz = V / 4
mu = t2 * (mu_eff + 6); hz = (mu - 4 * Jz) / 2
print(f"J_x = {Jx:.9f}, J_y = {Jy:.9f}, D = {D:.9f}, J_z = {Jz:.6f}, h_z = {hz:.9f}  (|t2| = {t2:.6g}, J_x+J_y = {s:.4e})")
print(f"sep-2 pair creation / |t2| = {2*Delta*s/V/t2:.3e}, sep-3 pair creation / |t2| = {2*Delta**3/(3*V**2)/t2:.3e}")
