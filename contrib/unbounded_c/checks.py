"""Independent checks for contrib/unbounded_c (numpy + scipy only).

(a) Exact diagonalization of the d = 4 period selector (C7 Theorem 4, rederived
    independently here with letters M, L, W, R = #, L, X, R) on a ring of L = 10
    at a coupling where the selected period is 5 > d.
(b) Reduced exact solution: selected period, gap to other cell lengths, internal gap.
(c) Siegel-Weil average count of Narain lattice vectors with Delta <= 2 (and <= 1),
    to locate the smallest c at which a Narain CFT with no relevant or marginal
    non-current primaries is guaranteed to exist.
(d) Second-order rank raising: a rank-2 port bond produces a second-order
    inter-block coupling of connected rank N - 1.

Run:  python3 checks.py      (about a minute; (a) dominates)
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.special import gammaln, zeta
from scipy.integrate import quad
from fractions import Fraction

# ---------------------------------------------------------------- (a) + (b)
M, Lc, W, R = 0, 1, 2, 3
ALLOWED = {(M, Lc), (M, W), (Lc, Lc), (Lc, W), (W, R), (W, M), (R, R), (R, M)}


def selector_hamiltonian(L, mu, t=1.0, K=8.0):
    dim = 4 ** L
    idx = np.arange(dim, dtype=np.int64)
    digits = np.stack([(idx // 4 ** i) % 4 for i in range(L)], axis=1)
    diag = mu * (digits == M).sum(axis=1).astype(float)
    rows, cols = [], []
    for i in range(L):
        j = (i + 1) % L
        a, b = digits[:, i], digits[:, j]
        bad = np.ones(dim, dtype=bool)
        for (x, y) in ALLOWED:
            bad &= ~((a == x) & (b == y))
        diag += K * bad
        src = np.nonzero((a == Lc) & (b == W))[0]          # L W -> W R
        dst = src + (W - Lc) * 4 ** i + (R - W) * 4 ** j
        rows += [dst, src]
        cols += [src, dst]
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    H = sp.coo_matrix((-t * np.ones(len(rows)), (rows, cols)), shape=(dim, dim))
    return (H + sp.diags(diag)).tocsr()


def cell_energy(w, mu, t=1.0):
    """Lowest energy of a cell: one marker plus w non-marker sites, one W walking."""
    return mu - 2 * t * np.cos(np.pi / (w + 1))


def reduced_spectrum(L, mu, t=1.0, nmax=6):
    """Lowest energies over all legal marker tilings of a ring of L sites
    (cells of period p = w + 1 >= 2), plus the two zero-energy marker-free states.
    Only ground states of each cell are used for the tiling minimum; internal
    excitations are reported separately."""
    best = {}
    # compositions of L into parts p >= 2 (unordered multiset suffices for energies)
    def rec(rem, minp, parts):
        if rem == 0:
            e = sum(cell_energy(p - 1, mu, t) for p in parts)
            best[tuple(parts)] = e
            return
        for p in range(minp, rem + 1):
            rec(rem - p, p, parts + [p])
    rec(L, 2, [])
    best[("no markers",)] = 0.0
    return sorted(best.items(), key=lambda kv: kv[1])[:nmax]


def selector_table(mus, t=1.0, wmax=4000):
    out = []
    w = np.arange(1, wmax)
    for mu in mus:
        dens = cell_energy(w, mu, t) / (w + 1)
        k = np.argmin(dens)
        p = w[k] + 1
        estar = dens[k]
        excess = cell_energy(w, mu, t) - (w + 1) * estar
        excess[k] = np.inf
        gamma = excess.min()
        omega = 2 * t * (np.cos(np.pi / p) - np.cos(2 * np.pi / p))
        out.append((mu, p, estar, gamma, omega, np.pi * np.sqrt(3 * t / (2 * t - mu))))
    return out


# ---------------------------------------------------------------- (c)
def narain_expected_count(c, D=2.0):
    """Siegel-Weil average over Narain moduli of the number of nonzero lattice
    vectors with Delta = h + hbar <= D.  Density of primaries:
    (2 pi)^c / Gamma(c/2)^2 (h hbar)^(c/2-1) a_j,  j = h - hbar in Z,
    a_0 = zeta(c-1)/zeta(c),  a_j = sigma_{1-c}(|j|)/zeta(c)."""
    pref = np.exp(c * np.log(2 * np.pi) - 2 * gammaln(c / 2))
    total = 0.0
    jmax = int(np.floor(D))
    for j in range(-jmax, jmax + 1):
        aj = zeta(c - 1) / zeta(c) if j == 0 else \
            sum(dd ** (1 - c) for dd in range(1, abs(j) + 1) if abs(j) % dd == 0) / zeta(c)
        J = abs(j)
        hmin, hmax = J, (D + J) / 2.0          # h >= |j| (taking h >= hbar), Delta = 2h - |j|
        if hmax <= hmin:
            continue
        val, _ = quad(lambda h: (h * (h - J)) ** (c / 2 - 1), hmin, hmax)
        total += aj * val
    return pref * total


def exact_connected_rank(N, eps):
    """Exact rank (rational arithmetic) of the connected coefficient matrix of
    sum_ab C_ab |a><a| x |b><b|, C_ab = 1/(2 + eps (a + b)): double-centre C,
    which removes the one-block parts, and row-reduce.  Cauchy matrices are
    badly conditioned, so floating-point SVD undercounts the rank for N >~ 8."""
    C = [[Fraction(1) / (2 + eps * (a + b)) for b in range(1, N + 1)] for a in range(1, N + 1)]
    rm = [sum(r) / N for r in C]
    cm = [sum(C[a][b] for a in range(N)) / N for b in range(N)]
    tm = sum(rm) / N
    A = [[C[a][b] - rm[a] - cm[b] + tm for b in range(N)] for a in range(N)]
    rank, col = 0, 0
    rows = len(A)
    for col in range(N):
        piv = next((r for r in range(rank, rows) if A[r][col] != 0), None)
        if piv is None:
            continue
        A[rank], A[piv] = A[piv], A[rank]
        for r in range(rows):
            if r != rank and A[r][col] != 0:
                f = A[r][col] / A[rank][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[rank])]
        rank += 1
    return rank


# ---------------------------------------------------------------- (d)
def second_order_rank_example(N, eps=0.3, Delta=1.0):
    """Block = port qubit (x) register C^N.  H_A = Delta |1><1|_port (x) (1 + eps E),
    E = diag(1..N).  Bond h = sigma^+ sigma^+ + h.c. on the two ports (connected
    rank 2).  First order vanishes on the ground space |0>|psi>.  Second order:
    H2 = - sum_{ab} |ab><ab| / (Delta (2 + eps(a + b)))."""
    E = np.diag(np.arange(1, N + 1, dtype=float))
    C = np.array([[1.0 / (Delta * (2 + eps * (E[a, a] + E[b, b]))) for b in range(N)]
                  for a in range(N)])
    # connected rank = rank of C after removing row/column means (one-block parts)
    rank = exact_connected_rank(N, Fraction(eps).limit_denominator(1000))
    # exact check by diagonalization at small lambda
    lam = 1e-3
    dimA = 2 * N
    P1 = np.diag([0.0, 1.0])
    HA = Delta * np.kron(P1, np.eye(N) + eps * E)
    H0 = np.kron(HA, np.eye(dimA)) + np.kron(np.eye(dimA), HA)
    sp_ = np.array([[0, 1], [0, 0]], dtype=float)   # |0><1|  (sigma^+ in the 0/1 labels below)
    up = sp_.T                                         # |1><0|
    hbond = np.kron(up, up) + np.kron(sp_, sp_)        # creates/annihilates a port pair
    # assemble V on (port_A, reg_A, port_B, reg_B): reorder via explicit loops
    dimtot = dimA * dimA
    Vfull = np.zeros((dimtot, dimtot))
    for pa in range(2):
        for pb in range(2):
            for qa in range(2):
                for qb in range(2):
                    amp = hbond[pa * 2 + pb, qa * 2 + qb]
                    if amp == 0:
                        continue
                    for ra in range(N):
                        for rb in range(N):
                            i = (pa * N + ra) * dimA + (pb * N + rb)
                            j = (qa * N + ra) * dimA + (qb * N + rb)
                            Vfull[i, j] += amp
    ev = np.linalg.eigvalsh(H0 + lam * Vfull)[: N * N]
    pred = np.sort(-lam ** 2 * C.ravel())
    return rank, np.max(np.abs(ev - pred)) / lam ** 2


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    print("(a) ED of the d = 4 selector, ring L = 10, t = 1, K = 8, mu = 0.9")
    H = selector_hamiltonian(10, 0.9)
    vals = np.sort(sla.eigsh(H, k=14, which="SA", return_eigenvectors=False, tol=1e-10))
    print("   lowest 14 eigenvalues:", vals)
    print("   reduced prediction (tilings):")
    for parts, e in reduced_spectrum(10, 0.9):
        print("     ", parts, round(e, 6))
    print("   predicted ground energy 2(mu - 2cos(pi/5)) =", 2 * (0.9 - 2 * np.cos(np.pi / 5)))

    print("\n(b) selected period versus mu (t = 1):  mu, period p, e*, gamma_p, omega_p, pi*sqrt(3/(2-mu))")
    for row in selector_table([0.9, 1.5, 1.8, 1.9, 1.95, 1.99, 1.999]):
        print("   ", "  ".join(f"{x:.6g}" for x in row))

    print("\n(c) Siegel-Weil expected number of nonzero Narain vectors")
    for c in [13, 15, 20, 25, 28, 29, 30, 31, 32, 33, 35, 40]:
        print(f"    c = {c:2d}:  E[#(Delta<=2)] = {narain_expected_count(c, 2.0):9.4f}   "
              f"E[#(Delta<=1)] = {narain_expected_count(c, 1.0):9.4g}")

    print("\n(d) second-order rank raising through a rank-2 port bond")
    for N in [3, 5, 8, 12]:
        r, err = second_order_rank_example(N)
        print(f"    N = {N:2d}: connected rank of H2 = {r} (= N - 1); max |exact - formula| / lambda^2 = {err:.2e}")
