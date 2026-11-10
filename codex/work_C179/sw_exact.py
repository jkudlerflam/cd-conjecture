"""Exact (all-orders) Schrieffer-Wolff effective Hamiltonian on the blockaded space, from the full qubit chain.

Independent of C173's code. The full open chain is built from Pauli matrices in the bond convention
    h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z(x)1 + 1(x)Z),
summed over bonds (so the bulk field is hz Z_j). Particle = spin down, n = (1 - Z)/2.
P = configurations with no two adjacent particles. The effective Hamiltonian on P is the des Cloizeaux /
direct-rotation form H_eff = U E U^dag, U = A (A^dag A)^(-1/2), A = <P basis | low eigenvectors>,
which equals the canonical Schrieffer-Wolff effective Hamiltonian (Bravyi-DiVincenzo-Loss 2011).
Matrix elements are in the spin (hard-core boson) basis, with no Jordan-Wigner signs.

Predictions (PROOF 1, spin language, a^dag = sigma^-):
  NNN hop  <j+2|H|j>          = -Delta^2/(2V)
  on-site  E(j) - E(0)        = -mu + 3 Delta^2/V            (bulk j)
  V2       E(j,j+2) - E(j) - E(j+2) + E(0) = -Delta^2/V
  V3       ... distance 3                 = -Delta^2/(3V)
  V4       ... distance 4                 = 0
  dN=2, separation 2: <b,b+2|H|0> = -2 Delta Re(w)/V
  dN=2, separation 3: <b,b+3|H|0> = +2 Delta^3/(3 V^2)
"""
import numpy as np, sys

X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2)


def op(L, ops):
    """ops: dict site -> 2x2. Site 0 is the most significant tensor factor; basis index bit (L-1-j) <-> site j."""
    M = np.array([[1.0 + 0j]])
    for j in range(L):
        M = np.kron(M, ops.get(j, I2))
    return M


def full_chain(L, Jx, Jy, D, Jz, hz):
    H = np.zeros((2 ** L, 2 ** L), complex)
    for j in range(L - 1):
        H += Jx * op(L, {j: X, j + 1: X}) + Jy * op(L, {j: Y, j + 1: Y}) + Jz * op(L, {j: Z, j + 1: Z})
        H += D * (op(L, {j: X, j + 1: Y}) - op(L, {j: Y, j + 1: X}))
        H += (hz / 2) * (op(L, {j: Z}) + op(L, {j + 1: Z}))
    return H


def occ(idx, L):
    """occupation list (particle = down = basis state |1> in Z eigenbasis ordering [up, down])."""
    return [(idx >> (L - 1 - j)) & 1 for j in range(L)]


def index_of(sites, L):
    s = 0
    for j in sites:
        s |= 1 << (L - 1 - j)
    return s


def heff(L, Delta, V, w, mu):
    Jz = V / 4; hz = (mu - 4 * Jz) / 2
    Jx = (w.real + Delta) / 2; Jy = (w.real - Delta) / 2; D = w.imag / 2
    H = full_chain(L, Jx, Jy, D, Jz, hz)
    # sanity: hz field convention. bulk field on site j is hz Z_j (two half-fields), edges get hz/2: correct the edges
    H += (hz / 2) * (op(L, {0: Z}) + op(L, {L - 1: Z}))
    # edge sites touch one ZZ bond only: add Jz Z on each edge so that every site has chemical potential mu
    H += Jz * (op(L, {0: Z}) + op(L, {L - 1: Z}))
    e, v = np.linalg.eigh(H)
    Pidx = [s for s in range(2 ** L) if not any(a & b for a, b in zip(occ(s, L)[:-1], occ(s, L)[1:]))]
    m = len(Pidx)
    gap = e[m] - e[m - 1]
    A = v[Pidx, :m]
    u, sv, vh = np.linalg.svd(A)
    U = u @ vh
    He = U @ np.diag(e[:m]) @ U.conj().T
    pos = {s: i for i, s in enumerate(Pidx)}
    return He, pos, gap, sv.min()


def elem(He, pos, L, f, i):
    return He[pos[index_of(f, L)], pos[index_of(i, L)]]


if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    Delta = 1.0
    print(f"L = {L}, Delta = 1. Units: t2 = Delta^2/(2V). Columns: measured, predicted, ratio.")
    for wt in (1j, 0.6 + 0.8j):
        for V in (20.0, 40.0, 80.0, 160.0):
            t2 = Delta ** 2 / (2 * V)
            w = wt * t2
            mu = -2.2 * t2 + 3 * Delta ** 2 / V
            He, pos, gap, smin = heff(L, Delta, V, w, mu)
            c = L // 2 - 1
            E0 = elem(He, pos, L, [], [])
            Ej = lambda *s: elem(He, pos, L, list(s), list(s))
            hop2 = elem(He, pos, L, [c + 2], [c])
            hop1 = elem(He, pos, L, [c + 1], [c])
            ons = Ej(c) - E0
            V2 = Ej(c, c + 2) - Ej(c) - Ej(c + 2) + E0
            V3 = Ej(c - 1, c + 2) - Ej(c - 1) - Ej(c + 2) + E0
            V4 = Ej(c - 1, c + 3) - Ej(c - 1) - Ej(c + 3) + E0
            p2 = elem(He, pos, L, [c, c + 2], [])
            p3 = elem(He, pos, L, [c - 1, c + 2], [])
            rows = [("hop1/w", hop1 / w, 1.0),
                    ("hop2", hop2, -Delta ** 2 / (2 * V)),
                    ("onsite+mu", ons.real + mu, 3 * Delta ** 2 / V),
                    ("V2", V2.real, -Delta ** 2 / V),
                    ("V3", V3.real, -Delta ** 2 / (3 * V)),
                    ("V4", V4.real, 0.0),
                    ("dN2 sep2", p2, -2 * Delta * w.real / V),
                    ("dN2 sep3", p3, 2 * Delta ** 3 / (3 * V ** 2))]
            print(f"\nw/t2 = {wt}, V = {V}: band gap = {gap:.3f}, min singular value of A = {smin:.5f}")
            for name, meas, pred in rows:
                r = meas / pred if abs(pred) > 0 else float('nan')
                print(f"   {name:10s} meas = {complex(meas):.6e}   pred = {pred:.6e}   ratio = {complex(r):.5f}")
