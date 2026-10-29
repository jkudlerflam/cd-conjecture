"""Kohn stiffness versus the bond: exact diagonalization checks.

For a nearest-neighbour chain H = sum_x h_{x,x+1} on a ring with a conserved
on-site charge q (diagonal in the computational basis), this script computes

  E_0(L)                 ring ground energy, Casimir fit for A = sum_a c_a v_a
  D_L(q)                 Kohn stiffness (1/L) d^2 E_0 / d theta^2 at theta = 0,
                         theta a uniform twist per bond (one bond carries L theta)
  <h''>                  ground-state expectation of -[q,[q,h]] on one bond
  ||h''||                operator norm of [q,[q,h~]] for the centred bond h~
  M                      bare half-spread of the bond (gauge does not lower it
                         for the chains used here, see fable_bound_proof)

and prints the chain of inequalities  D <= <h''> <= ||h''|| <= sigma_q^2 M,
together with the inferred bound on A in each case.

Only numpy and scipy. Run:  python stiffness_ed.py [quick]
"""
import sys
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ============================================================================
# bonds
# ============================================================================

def spin_matrices(S):
    dim = int(round(2 * S + 1))
    m = np.arange(S, -S - 1, -1)
    Sz = np.diag(m)
    Sp = np.zeros((dim, dim))
    for i in range(1, dim):
        Sp[i - 1, i] = np.sqrt(S * (S + 1) - m[i] * (m[i] + 1))
    Sm = Sp.T
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2j)
    return Sx, Sy, Sz


def bond_xxz(Delta):
    Sx, Sy, Sz = spin_matrices(0.5)
    h = np.kron(Sx, Sx) + np.kron(Sy, Sy) + Delta * np.kron(Sz, Sz)
    return h.real, np.diag(Sz).real


def bond_swap(N):
    d = N
    h = np.zeros((d * d, d * d))
    for a in range(d):
        for b in range(d):
            h[a * d + b, b * d + a] = 1.0
    return h, None


def bond_bilinear_biquadratic(beta):
    """spin-1 bond S.S + beta (S.S)^2 ; beta = -1 is Takhtajan-Babujian
    (SU(2)_2, c = 3/2), beta = +1 is Uimin-Lai-Sutherland (SU(3)_1, c = 2)."""
    Sx, Sy, Sz = spin_matrices(1.0)
    SS = np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)
    SS = SS.real
    h = SS + beta * SS @ SS
    return h, np.diag(Sz).real


# ============================================================================
# sector basis and sparse Hamiltonian
# ============================================================================

def sector_states(L, d, charge, target):
    """all basis states (as integer arrays of digits) with sum charge = target"""
    # enumerate digit strings via recursion on counts to avoid d^L blowup
    states = []
    # dynamic programming over positions with running charge
    # charges are small integers here (multiplied to integers)
    for digits in itertools.product(range(d), repeat=L):
        if abs(sum(charge[x] for x in digits) - target) < 1e-9:
            states.append(digits)
    return np.array(states, dtype=np.int8)


def sector_states_fast(L, d, charge, target):
    """vectorized enumeration: build in blocks by site"""
    # represent states as integers base d; keep running charge
    codes = np.zeros(1, dtype=np.int64)
    chg = np.zeros(1)
    for x in range(L):
        codes = (codes[:, None] * d + np.arange(d)[None, :]).ravel()
        chg = (chg[:, None] + charge[None, :]).ravel()
        # prune: remaining sites can change charge by at most (L-x-1)*range
        rem = L - x - 1
        lo = chg + rem * charge.min()
        hi = chg + rem * charge.max()
        keep = (lo <= target + 1e-9) & (hi >= target - 1e-9)
        codes = codes[keep]
        chg = chg[keep]
    keep = np.abs(chg - target) < 1e-9
    return np.sort(codes[keep])


def digits_of(codes, L, d):
    out = np.zeros((len(codes), L), dtype=np.int64)
    c = codes.copy()
    for x in range(L - 1, -1, -1):
        out[:, x] = c % d
        c //= d
    return out


def build_H(codes, L, d, h, charge, Phi=0.0, bond_ops=None):
    """Ring Hamiltonian in the sector spanned by codes. The twist Phi is put on
    the bond (L-1, 0): h -> e^{i Phi q_{L-1}} h e^{-i Phi q_{L-1}}.
    Optionally bond_ops = list of (bond_index, matrix) to replace h on given bonds."""
    n = len(codes)
    dig = digits_of(codes, L, d)
    rows, cols, vals = [], [], []
    hh = h.reshape(d, d, d, d)  # (a', b', a, b)
    nz = np.argwhere(np.abs(hh) > 1e-14)
    for x in range(L):
        y = (x + 1) % L
        a = dig[:, x]
        b = dig[:, y]
        hx = hh
        if bond_ops is not None:
            for bi, mat in bond_ops:
                if bi == x:
                    hx = mat.reshape(d, d, d, d)
        nzx = np.argwhere(np.abs(hx) > 1e-14)
        for (ap, bp, aa, bb) in nzx:
            sel = np.nonzero((a == aa) & (b == bb))[0]
            if len(sel) == 0:
                continue
            val = hx[ap, bp, aa, bb]
            if x == L - 1 and Phi != 0.0 and charge is not None:
                val = val * np.exp(1j * Phi * (charge[ap] - charge[aa]))
            newcodes = codes[sel] + (ap - aa) * d ** (L - 1 - x) + (bp - bb) * d ** (L - 1 - y)
            idx = np.searchsorted(codes, newcodes)
            rows.append(idx)
            cols.append(sel)
            vals.append(np.full(len(sel), val))
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    dtype = complex if np.iscomplexobj(vals) else float
    H = sp.csr_matrix((vals.astype(dtype), (rows, cols)), shape=(n, n))
    return H


def ground(H, k=1):
    if H.shape[0] <= 400:
        w, v = np.linalg.eigh(H.toarray())
        return w[0], v[:, 0]
    w, v = spla.eigsh(H, k=k, which='SA', tol=1e-11, maxiter=20000)
    i = np.argmin(w)
    return w[i], v[:, i]


def bond_expectation(codes, L, d, psi, O):
    """<psi| O_{0,1} |psi> for a two-site operator O on bond (0,1)"""
    H01 = build_H(codes, L, d, np.zeros_like(O), None, 0.0, bond_ops=[(0, O)])
    return np.vdot(psi, H01 @ psi).real


# ============================================================================
# analysis helpers
# ============================================================================

def half_spread(h):
    w = np.linalg.eigvalsh(h)
    return (w[-1] - w[0]) / 2


def centre(h):
    w = np.linalg.eigvalsh(h)
    return h - (w[-1] + w[0]) / 2 * np.eye(len(h))


def double_commutator(h, charge):
    d = len(charge)
    Q = np.kron(np.diag(charge), np.eye(d))
    return Q @ Q @ h + h @ Q @ Q - 2 * Q @ h @ Q  # = [Q,[Q,h]]


def casimir_fit(Ls, E0s):
    """fit E_0/L = e_0 - a/L^2 + b/L^4 ; A = 6 a / pi"""
    Ls = np.array(Ls, float)
    y = np.array(E0s) / Ls
    X = np.vstack([np.ones_like(Ls), -1 / Ls ** 2, 1 / Ls ** 4]).T
    if len(Ls) < 3:
        X = X[:, :2]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 6 * coef[1] / np.pi, coef[0]


def stiffness(codes, L, d, h, charge, E0, delta=0.02):
    """D_L = L * d^2 E_0 / dPhi^2 with Phi the twist on one bond (= (1/L) d^2E/dtheta^2
    for a uniform twist theta per bond). Symmetric finite difference, two step sizes
    for a Richardson check."""
    out = []
    for dl in (delta, 2 * delta):
        Ep, _ = ground(build_H(codes, L, d, h, charge, Phi=dl))
        Em, _ = ground(build_H(codes, L, d, h, charge, Phi=-dl))
        out.append(L * (Ep + Em - 2 * E0) / dl ** 2)
    # Richardson: D(delta) has O(delta^2) error
    return (4 * out[0] - out[1]) / 3, out


# ============================================================================
# main table
# ============================================================================

def run_chain(name, h, charge, d, Ls, target_fn, K_or_note):
    print("=" * 78)
    print(name)
    M = half_spread(h)
    ht = centre(h)
    C = double_commutator(ht, charge)
    normC = np.linalg.norm(C, 2)
    sigma = charge.max() - charge.min()
    print(f"  d = {d}, M (bare half-spread) = {M:.6f}, ||[q,[q,h~]]|| = {normC:.6f}, "
          f"sigma_q^2 M = {sigma**2 * M:.6f}")
    E0s = []
    rows = []
    for L in Ls:
        codes = sector_states_fast(L, d, charge, target_fn(L))
        H = build_H(codes, L, d, h, charge)
        E0, psi = ground(H)
        E0s.append(E0)
        hpp = -bond_expectation(codes, L, d, psi, C)  # <h''> = -<[q,[q,h]]>
        D, raw = stiffness(codes, L, d, h, charge, E0)
        rows.append((L, E0, D, hpp))
        print(f"  L = {L:2d}  E0/L = {E0/L:+.8f}  D_L = {D:.6f}  (raw {raw[0]:.6f}, {raw[1]:.6f})"
              f"  <h''>_L = {hpp:.6f}   D/<h''> = {D/hpp:.4f}")
    if len(Ls) >= 3:
        A, e0 = casimir_fit(Ls[-4:], E0s[-4:])
        print(f"  Casimir fit over L = {Ls[-4:]}: A = {A:.5f}, e_0 = {e0:.6f}")
    else:
        A = float('nan')
    print(f"  note: {K_or_note}")
    return rows, A, M, normC


if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    Lmax2 = 14 if quick else 18
    even = lambda L: 0.0

    # XXZ family, charge S^z, sector S^z_tot = 0
    for Delta in (1.0, 0.5, 0.0, -0.5):
        h, q = bond_xxz(Delta)
        gamma = np.arccos(Delta)
        v = np.pi * np.sin(gamma) / (2 * gamma) if gamma > 0 else np.pi / 2
        K = np.pi / (2 * (np.pi - gamma))
        Dexact = v * K / np.pi
        run_chain(f"XXZ Delta = {Delta}", h, q, 2, list(range(8, Lmax2 + 1, 2)), even,
                  f"Bethe: v = {v:.5f}, K = {K:.5f}, D = vK/pi = {Dexact:.5f}; "
                  f"stiffness bound A = pi D/K <= pi <h''>/K ; norm bound (pi/K) ||h''|| = {np.pi/K*np.linalg.norm(double_commutator(centre(h), q), 2):.4f}")

    # SU(3) swap chain, Cartan charge lambda = (1,-1,0)/sqrt2 and (1,1,-2)/sqrt6
    h, _ = bond_swap(3)
    for lam, lname in (((1, -1, 0) / np.sqrt(2), "lambda_3-type"), ((1, 1, -2) / np.sqrt(6), "lambda_8-type")):
        q = np.array(lam)
        run_chain(f"SU(3) swap, charge {lname} {np.round(q,4)}", h, q, 3, [6, 9, 12], even,
                  "SU(3)_1: v = 2pi/3, c = 2, D_alpha = v/pi = 0.66667 per orthonormal Cartan direction")

    # spin-1 bilinear-biquadratic: TB point (beta=-1) and ULS point (beta=+1)
    Ls1 = [6, 8, 10] if quick else [6, 8, 10, 12]
    Ls_uls = [6, 9, 12]
    for beta, note, Lset in ((-1.0, "Takhtajan-Babujian SU(2)_2, c = 3/2, v = 2 pi in these units (A = 3 pi)", Ls1),
                             (1.0, "ULS = SU(3)_1 in spin-1 language (h = P + 1), c = 2, v = 2pi/3; S^z = sqrt2 x orthonormal Cartan charge", Ls_uls)):
        h, q = bond_bilinear_biquadratic(beta)
        run_chain(f"spin-1 S.S + ({beta}) (S.S)^2, charge S^z", h, q, 3, Lset, even, note)
