"""Referee rebuild of h9 and h10 from the README text alone.

Written before opening bond9.py / bond10.py.  Letter order is the order in
which the README lists the letters:

    0 L_0, 1 L_1, 2 R_0, 3 R_1, 4 W, 5 A, 6 B, 7 X_A, 8 X_B, (9 V)

Interpretation choices made where the README text leaves room (all listed in
REFEREE.md):
  * "+ h.c." at the end of the display for w is applied to every line of w.
  * the on-site part is (1/2)(o x 1 + 1 x o).
  * h10: J*(w + K Pi_bad) is embedded in the material-material block only;
    the on-site operator (J*o - mu*n_occ) is applied as a genuine on-site term
    (1/2)(O x 1 + 1 x O) on all 100 pairs, including pairs with a vacancy;
    vacancy hop -t (|a V><V a| + h.c.) for every material a.
"""
import numpy as np

L0, L1, R0, R1, W, A, B, XA, XB, V = range(10)
NAMES = ["L0", "L1", "R0", "R1", "W", "A", "B", "XA", "XB", "V"]
Lsym = (L0, L1)
Rsym = (R0, R1)
HEADS = (W, A, B, XA, XB)
ACT = {"A": A, "B": B}
BLIND = {"A": XA, "B": XB}
CHARGE = np.array([0, 2, 0, 2, 0, -1, +1, -1, +1, 0], dtype=float)
N_L = (R0, R1, W, A, B, XA, XB)   # allowed left neighbours for the rail switch
N_R = (L0, L1, W, A, B, XA, XB)   # allowed right neighbours for the closing link


def _w(d, par):
    """Off-diagonal two-site operator w (including + h.c.) on C^d x C^d."""
    phi = {"A": par["phiA"], "B": 0.0}
    psi = {"A": par["psiA"], "B": par["psiB"]}
    z = {"A": lambda a: (-1) ** a, "B": lambda a: 1}
    half = np.zeros((d * d, d * d), dtype=complex)

    def add(l1, l2, r1, r2, c):
        # c |l1 l2><r1 r2|
        half[l1 * d + l2, r1 * d + r2] += c

    for a in (0, 1):
        La, Ra = Lsym[a], Rsym[a]
        add(La, W, W, Ra, -1.0)
        for s in "AB":
            add(La, ACT[s], ACT[s], Ra, -z[s](a))
            add(La, BLIND[s], BLIND[s], Ra, -z[s](a) * np.exp(1j * phi[s]))
    for s in "AB":
        for y in N_L:
            add(y, BLIND[s], y, ACT[s], -1.0)
        for y in N_R:
            add(BLIND[s], y, ACT[s], y, -par["td"] * np.exp(1j * psi[s]))
    add(B, R0, A, R1, par["g"])
    add(XA, B, W, W, par["lam"] * np.exp(-1j * np.pi / 4))
    add(XB, A, W, W, par["lam"] * np.exp(+1j * np.pi / 4))
    return half + half.conj().T


def _o(d, par):
    o = np.zeros(d)
    for h in HEADS:
        o[h] += par["nu"]
    for s in "AB":
        o[ACT[s]] += par["Delta"]
        o[BLIND[s]] += par["Delta"] + par["D"]
    o[:d] += 0.5 * par["f"] * CHARGE[:d]
    return o


def _bad(d):
    P = np.zeros((d * d, d * d))
    for a in Lsym:
        for b in Rsym:
            P[a * d + b, a * d + b] = 1.0
    return P


def h9(par):
    d = 9
    o = np.diag(_o(d, par))
    I = np.eye(d)
    return _w(d, par) + 0.5 * (np.kron(o, I) + np.kron(I, o)) + par["K"] * _bad(d)


def h10(par):
    d = 10
    mat9 = _w(9, par) + par["K"] * _bad(9)          # material-material only
    H = np.zeros((100, 100), dtype=complex)
    for i in range(9):
        for j in range(9):
            for k in range(9):
                for l in range(9):
                    H[i * 10 + j, k * 10 + l] = par["J"] * mat9[i * 9 + j, k * 9 + l]
    o10 = np.zeros(10)
    o10[:9] = par["J"] * _o(9, par) - par["mu"]
    O = np.diag(o10)
    I = np.eye(10)
    H += 0.5 * (np.kron(O, I) + np.kron(I, O))
    for a in range(9):
        H[a * 10 + V, V * 10 + a] += -par["t"]
        H[V * 10 + a, a * 10 + V] += -par["t"]
    return H


def random_par(rng, mobile=False):
    par = dict(
        nu=rng.normal(), Delta=rng.normal(), D=rng.normal(), f=rng.normal(),
        K=rng.normal(), td=rng.normal(), g=rng.normal(), lam=rng.normal(),
        phiA=rng.normal(), psiA=rng.normal(), psiB=rng.normal(),
    )
    if mobile:
        par.update(J=rng.normal(), t=rng.normal(), mu=rng.normal())
    return par


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    p = random_par(rng, True)
    H9 = h9(p)
    H10 = h10(p)
    print("h9 hermitian:", np.abs(H9 - H9.conj().T).max())
    print("h10 hermitian:", np.abs(H10 - H10.conj().T).max())
    off = H9 - np.diag(np.diag(H9))
    print("h9 nonzero offdiag:", int((np.abs(off) > 1e-14).sum()))
