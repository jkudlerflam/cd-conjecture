"""Bravyi-Gosset gapless frustration-free bonds on rings.
psi = sqrt(1-p)|00> + sqrt(p)|11> (T_psi eigenvalues +-i sqrt(p(1-p)), gapless for 0<p<1 on open chains),
and a complex example with a non-normal T. Ring spectra: ground degeneracy and gap scaling."""
import numpy as np
from ed import spectrum, all_levels, gauge_M

def proj(psi):
    psi = np.asarray(psi, complex); psi = psi / np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def T_of(psi):
    # T = [[<psi|01>, <psi|11>], [-<psi|00>, -<psi|10>]], basis index 2*b1+b2
    c = np.conj(np.asarray(psi, complex))
    return np.array([[c[1], c[3]], [-c[0], -c[2]]])

def report(name, psi, Ls):
    h = proj(psi)
    mu = np.linalg.eigvals(T_of(psi))
    print(f"{name}: |mu| = {np.abs(mu)}, arg ratio = {np.angle(mu[1]/mu[0]):.4f}, M = {gauge_M(h)[0]:.4f}")
    for L in Ls:
        w, k = all_levels(spectrum(h, L))
        w = np.sort(w)
        deg = np.sum(w < 1e-9)
        gap = w[deg]
        print(f"   L={L}: E0={w[0]:.2e}, ground deg={deg}, gap={gap:.5f}, L*gap={L*gap:.4f}, L^2*gap={L*L*gap:.4f}")

p = 0.2
report("real psi p=0.2", [np.sqrt(1-p), 0, 0, np.sqrt(p)], [6, 7, 8, 9, 10, 12, 14])
report("singlet (FM Heisenberg)", [0, 1, -1, 0], [6, 8, 10, 12, 14])
# non-normal example: A (x) A applied to twisted singlet with twist pi/3
A = np.array([[1.0, 0.7], [0.0, 1.0]])
phi = np.pi / 3
psi0 = np.array([0, 1, -np.exp(1j * phi), 0])
psi = np.kron(np.linalg.inv(A).conj().T, np.linalg.inv(A).conj().T) @ psi0
report("deformed twisted singlet (twist pi/3, A=[[1,.7],[0,1]])", psi, [6, 8, 9, 10, 12])
