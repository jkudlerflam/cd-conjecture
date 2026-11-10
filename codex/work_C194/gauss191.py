"""C191: scaling dimensions in a multicomponent chiral Gaussian theory.

Chiral bosons phi_P (P = aR, aL, bR, bL), chirality s_P = +1 (R), -1 (L), psi_P ~ e^{i phi_P}.
Zero-mode energy E = (pi/L) N^T W N, with W_PP = |v_P| and W_PP' = f_PP'/(2 pi) (Landau parameters f).
Field Hamiltonian (1/4 pi) dphi^T Wf dphi with Wf = S W S (since rho_P = s_P dphi_P / 2 pi).
Normal modes: rows l_m of R are left eigenvectors of S Wf, normalized l_m S l_m^T = +-1.
A vertex e^{i alpha.phi} has scaling dimension x = (1/2) |R^{-T} alpha|^2.
"""
import numpy as np

S4 = np.diag([1.0, -1.0, 1.0, -1.0])        # basis (aR, aL, bR, bL)


def metric(W, S=S4):
    # W is the charge-energy matrix, E = (pi/L) N^T W N.  With rho_P = s_P d(phi_P)/(2 pi) the field
    # Hamiltonian is (1/4 pi) dphi^T (S W S) dphi, so the normal modes use Wf = S W S.
    Wf = S @ W @ S
    lam, Lv = np.linalg.eig((S @ Wf).T)        # columns: left eigenvectors of S W
    if np.max(np.abs(lam.imag)) > 1e-10:
        raise ValueError("unstable quadratic form (complex velocities)")
    lam = lam.real
    Lv = Lv.real
    R = []
    sig = []
    for m in range(len(lam)):
        l = Lv[:, m]
        nrm = l @ S @ l
        R.append(l / np.sqrt(abs(nrm)))
        sig.append(np.sign(nrm))
    R = np.array(R)
    G = np.linalg.inv(R.T @ R)                # x(alpha) = 1/2 alpha^T G alpha
    return G, lam, np.array(sig)


def dim(G, alpha):
    a = np.asarray(alpha, float)
    return 0.5 * a @ G @ a


# charge vectors (aR, aL, bR, bL) for e^{i(a theta + b phi)} per valley: Q_R = (a-b)/2, Q_L = (a+b)/2
def charges(ta=0, fa=0, tb=0, fb=0):
    return np.array([(ta - fa) / 2, (ta + fa) / 2, (tb - fb) / 2, (tb + fb) / 2])


OPS = {
    'O1 = cos 2(phi_a - phi_b)': charges(fa=2, fb=-2),          # x = 2 K_s
    'pair Psi_a Psi_b = e^{i(th_a+th_b)}': charges(ta=1, tb=1),   # x = 1/(2 K_c)
    'boson Psi_a = e^{i th_a}': charges(ta=1),                    # x = (1/K_c + 1/K_s)/8
    'fermion e^{i th_a + i(phi_a - phi_b)}': charges(ta=1, fa=1, fb=-1),
    'fermion e^{i th_a + i(phi_a + phi_b)}': charges(ta=1, fa=1, fb=1),
    'charge 2kF: cos 2(phi_a + phi_b)': charges(fa=2, fb=2),      # x = 2 K_c
}


def W_from(vR, vL, f_par, f_pp, f_mm, f_4):
    """vR > 0 speed of aR (= bL), vL > 0 speed of aL (= bR).
    f_par: (aR,aL) = (bR,bL); f_pp: (aR,bL); f_mm: (aL,bR); f_4: (aR,bR) = (aL,bL)."""
    W = np.diag([vR, vL, vL, vR]).astype(float)
    def put(i, j, f):
        W[i, j] = W[j, i] = f / (2 * np.pi)
    put(0, 1, f_par); put(2, 3, f_par)
    put(0, 3, f_pp); put(1, 2, f_mm)
    put(0, 2, f_4); put(1, 3, f_4)
    return W


if __name__ == "__main__":
    # tests
    W = W_from(1.0, 1.0, 0, 0, 0, 0)
    G, lam, sig = metric(W)
    for k, a in OPS.items():
        print("free", k, dim(G, a))
    # one-valley LL with g2: K = sqrt((1 - y/2)/(1 + y/2)), y = g2/(pi v): check O1-like and theta dims for decoupled valleys
    v, g2 = 1.0, -0.8
    W = W_from(v, v, g2, 0, 0, 0)
    G, lam, sig = metric(W)
    y = g2 / (np.pi * v)
    K = np.sqrt((1 - y / 2) / (1 + y / 2))
    print("decoupled valleys, K_formula", K, " x(e^{i th_a}) =", dim(G, charges(ta=1)), "expected", 1 / (4 * K),
          " x(e^{2i phi_a}) =", dim(G, charges(fa=2)), "expected", K)
    # Galilean drift: vR != vL free -> dims unchanged
    G, lam, sig = metric(W_from(0.8, 1.2, 0, 0, 0, 0))
    print("free with drift: x(O1) =", dim(G, OPS['O1 = cos 2(phi_a - phi_b)']), " velocities", lam)
    # SU(2) GY-like: f_par = 0, f_pp = f_mm = f_4 = g
    g = -0.3
    G, lam, sig = metric(W_from(1, 1, 0, g, g, g))
    print("GY first order g =", g, ": K_s from O1 =", dim(G, OPS['O1 = cos 2(phi_a - phi_b)']) / 2, "expected 1+g/(2 pi v) =", 1 + g / (2 * np.pi))
