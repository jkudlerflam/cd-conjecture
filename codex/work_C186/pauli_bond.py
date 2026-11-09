"""C186: dictionary between the boson bond (PROOF eq. 0.1) and Pauli couplings (PROOF eq. 0.3), checked on rings.

Given cusp data (Delta, V, u, phi, mu_eff, v) for the transverse-DM bond Gamma (X Z - Z X), return the Pauli bond
    h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z1 + Z2) + Gxz XZ + Gzx ZX
whose second-order model is (3.1) with |w_eff| = u|t2|, arg w_eff = phi, V2 = -2 + v, r = sqrt(v), lambda_eff = 0.
"""
import numpy as np
from sw_general import bond_matrix, site_matrix, build_full

X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
# basis order in sw_general: index 0 = empty (spin up), 1 = occupied (spin down); Z = diag(1, -1) matches n = (1 - Z)/2


def pauli_bond(c):
    k = np.kron
    h = (c.get('Jx', 0) * k(X, X) + c.get('Jy', 0) * k(Y, Y) + c.get('D', 0) * (k(X, Y) - k(Y, X))
         + c.get('K', 0) * (k(X, Y) + k(Y, X)) + c.get('Jz', 0) * k(Z, Z)
         + c.get('hz', 0) / 2 * (k(Z, np.eye(2)) + k(np.eye(2), Z))
         + c.get('hx', 0) / 2 * (k(X, np.eye(2)) + k(np.eye(2), X))
         + c.get('hy', 0) / 2 * (k(Y, np.eye(2)) + k(np.eye(2), Y))
         + c.get('Gxz', 0) * k(X, Z) + c.get('Gzx', 0) * k(Z, X) + c.get('Gyz', 0) * k(Y, Z) + c.get('Gzy', 0) * k(Z, Y))
    return h


def boson_from_pauli(c):
    """(V, mu, w, p, beta, betap, lam) of PROOF eq. 0.1, 0.3, with lam the full unconditional creation amplitude."""
    V = 4 * c.get('Jz', 0)
    mu = 4 * c.get('Jz', 0) + 2 * c.get('hz', 0)
    w = c.get('Jx', 0) + c.get('Jy', 0) + 2j * c.get('D', 0)
    p = c.get('Jx', 0) - c.get('Jy', 0) + 2j * c.get('K', 0)
    gR = c.get('Gxz', 0) + 1j * c.get('Gyz', 0)    # multiplies b_j^dag (1 - 2 n_{j+1})
    gL = c.get('Gzx', 0) + 1j * c.get('Gzy', 0)    # multiplies b_j^dag (1 - 2 n_{j-1})
    lam = c.get('hx', 0) + 1j * c.get('hy', 0) + gR + gL
    beta, betap = -2 * gR, -2 * gL
    return V, mu, w, p, beta, betap, lam


def cusp_bond(Delta, V, u, phi, mu_eff, v):
    t2 = Delta ** 2 / (2 * V)
    d = Delta * np.sqrt(v)              # f_L - f_R, real
    G = d / 4                           # Gxz = G, Gzx = -G  ->  f_R = -d/2, f_L = +d/2, s = 0, lam = 0
    fL, fR = d / 2, -d / 2
    w_eff = u * t2 * np.exp(1j * phi)
    w = w_eff + fL * np.conj(fR) / V    # w_eff = w - fL fR*/V
    # bulk on-site: -mu + 3 Delta^2/V - (fL^2 + fR^2)/V = -mu_eff t2
    mu = mu_eff * t2 + 3 * Delta ** 2 / V - (fL ** 2 + fR ** 2) / V
    Jz = V / 4; hz = (mu - 4 * Jz) / 2
    Jx = (w.real + Delta) / 2; Jy = (w.real - Delta) / 2; D = w.imag / 2
    return dict(Jx=Jx, Jy=Jy, D=D, Jz=Jz, hz=hz, Gxz=G, Gzx=-G)


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    L = 8
    for trial in range(3):
        c = {k: rng.normal() for k in ('Jx', 'Jy', 'D', 'K', 'Jz', 'hz', 'hx', 'hy', 'Gxz', 'Gzx', 'Gyz', 'Gzy')}
        V, mu, w, p, beta, betap, lam = boson_from_pauli(c)
        Hb = build_full(L, bond_matrix(V, w, p, beta, betap), site_matrix(mu, lam), ring=True)
        Hp = build_full(L, pauli_bond(c), np.zeros((2, 2)), ring=True)
        D_ = Hp - Hb
        const = np.trace(D_).real / D_.shape[0]
        print(f"ring L={L}: max |H_pauli - H_boson - const| = {np.abs(D_ - const * np.eye(D_.shape[0])).max():.2e}")
    b = cusp_bond(1.0, 40.0, 1.0, np.pi / 2, -2.2, 0.25)
    print("example cusp bond (Delta=1, V=40, u=1, phi=pi/2, mu_eff=-2.2, v=0.25):", {k: round(x, 8) for k, x in b.items()})
