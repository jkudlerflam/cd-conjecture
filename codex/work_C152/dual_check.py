"""Check the exact duality of the Ashkin-Teller plane (cos phi = 0, h_z = 0): (t, Delta, V) -> (Delta, t, -V),
implemented by gamma^2_j -> (-1)^j gamma^2_j. Compares full ring spectra (both Pi sectors) at L = 8, 10."""
import numpy as np
from model import bond, spin_from_fermion, ham, sector_states
def spec(p, L):
    h = bond(*p); ev = []
    for par in (0, 1):
        H = ham(h, L, par, True).toarray(); ev += list(np.linalg.eigvalsh(H))
    return np.sort(ev)
for (s, V) in [(0.9, -0.8), (0.5, 0.6), (0.3, -1.3)]:
    p1 = spin_from_fermion(1.0, np.pi / 2, s, V, V)        # t = 1, Delta = s, h_z = 0 (mu = V)
    p2 = spin_from_fermion(s, np.pi / 2, 1.0, -V, -V)      # dual: t = s, Delta = 1, V -> -V
    for L in (8, 10):
        a, b = spec(p1, L), spec(p2, L)
        print(f"s={s} V={V} L={L}: max|E - E_dual| = {np.max(np.abs(a - b)):.2e}, E0 = {a[0]:.6f}, {b[0]:.6f}")
