"""C186: second-order Schrieffer-Wolff reduction of the most general nearest-neighbour bond near the Fibonacci point,
checked against the full qubit chain (R1 of RESULT.md).

Spin language: particle = spin down, b^dag = sigma^-, bit j of the integer state = occupation of site j.
Bond (sites 1,2 = j, j+1), basis index 2*n1 + n2:
    h = V n1 n2 + [w b2^dag b1 + p b1^dag b2^dag + beta b1^dag n2 + betap b2^dag n1 + h.c.]
Site terms: -mu n + lam b^dag + conj(lam) b.
So the creation amplitude at site j is c_j = lam + betap n_{j-1} + beta n_{j+1}:
f_L = lam + betap (left neighbour occupied), f_R = lam + beta (right neighbour occupied).

Part A: exact second-order matrix  K = P H1 P - P H1 Q H0^-1 Q H1 P  (H0 = V sum n n), at w = lam = 0,
        against the closed form of PROOF 2 (rules), entry by entry, open chains and rings, L = 8..12.
Part B: des Cloizeaux (all orders) effective Hamiltonian of the full chain at V -> infinity with
        w, lam, mu of order eps^2/V, bulk coefficients against the closed form, relative deviations vs V.
"""
import numpy as np, itertools, sys

bd = np.array([[0, 0], [1, 0]], complex)   # b^dag : |0> -> |1>
b = bd.T.copy()
nn = np.diag([0.0, 1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def bond_matrix(V, w, p, beta, betap):
    h = V * np.kron(nn, nn)
    h = h + w * np.kron(b, bd) + p * np.kron(bd, bd) + beta * np.kron(bd, nn) + betap * np.kron(nn, bd)
    h = h + (w * np.kron(b, bd) + p * np.kron(bd, bd) + beta * np.kron(bd, nn) + betap * np.kron(nn, bd)).conj().T
    return h


def site_matrix(mu, lam):
    return -mu * nn + lam * bd + np.conj(lam) * b


def bit(s, j):
    return (s >> j) & 1


def build_full(L, hb, hs, ring=False):
    D = 2 ** L
    H = np.zeros((D, D), complex)
    bonds = [(j, j + 1) for j in range(L - 1)] + ([(L - 1, 0)] if ring else [])
    for s in range(D):
        for (i, j) in bonds:
            loc = 2 * bit(s, i) + bit(s, j)
            base = s & ~((1 << i) | (1 << j))
            for out in range(4):
                a = hb[out, loc]
                if a != 0:
                    t = base | ((out >> 1) << i) | ((out & 1) << j)
                    H[t, s] += a
        for i in range(L):
            loc = bit(s, i)
            base = s & ~(1 << i)
            for out in range(2):
                a = hs[out, loc]
                if a != 0:
                    H[base | (out << i), s] += a
    return H


def adj_count(s, L, ring):
    c = 0
    for j in range(L - 1 + (1 if ring else 0)):
        c += bit(s, j) & bit(s, (j + 1) % L)
    return c


def pspace(L, ring):
    return [s for s in range(2 ** L) if adj_count(s, L, ring) == 0]


def second_order_exact(L, V, w, p, beta, betap, mu, lam, ring=False):
    """K = P H1 P - P H1 Q H0^{-1} Q H1 P, exactly, as a matrix on P."""
    H0d = np.array([V * adj_count(s, L, ring) for s in range(2 ** L)], float)
    Hfull = build_full(L, bond_matrix(V, w, p, beta, betap), site_matrix(mu, lam), ring)
    H1 = Hfull - np.diag(H0d)
    P = pspace(L, ring)
    Qm = H0d > 0
    A = H1[np.ix_(Qm, P)]                       # Q H1 P
    K = H1[np.ix_(P, P)] - A.conj().T @ (A / H0d[Qm][:, None])
    return K, P


def closed_form(L, V, p, beta, betap, mu, ring=False):
    """Closed form of PROOF 2 at w = lam = 0, built by the local rules."""
    P = pspace(L, ring)
    pos = {s: i for i, s in enumerate(P)}
    K = np.zeros((len(P), len(P)), complex)

    def n(s, j):
        if ring:
            return bit(s, j % L)
        return bit(s, j) if 0 <= j < L else 0

    def exists(j):
        return ring or (0 <= j < L)

    pp = abs(p) ** 2
    for s in P:
        i = pos[s]
        N = bin(s).count('1')
        e = -mu * N
        # pairing, diagonal
        nb = L if ring else L - 1
        for k in range(nb):
            k1 = (k + 1) % L
            if n(s, k) == 0 and n(s, k1) == 0:
                cnt = 1 + n(s, k - 1) + n(s, k + 2)
                e += -pp / (V * cnt)
        # facilitated, diagonal
        for j in range(L):
            if n(s, j) == 0:
                m = n(s, j - 1) + n(s, j + 1)
                if m > 0:
                    amp = betap * n(s, j - 1) + beta * n(s, j + 1)
                    e += -abs(amp) ** 2 / (V * m)
        K[i, i] += e
        for j in range(L):
            jm, jp, jp2, jp3 = (j - 1) % L, (j + 1) % L, (j + 2) % L, (j + 3) % L
            # NNN hop j -> j+2 (and the reverse via hermiticity below)
            if exists(j + 2) and n(s, j) == 1 and n(s, j + 1) == 0 and n(s, j + 2) == 0 and n(s, j + 3) == 0:
                t = (s & ~(1 << j)) | (1 << jp2)
                K[pos[t], i] += -pp / (2 * V)
                K[i, pos[t]] += -pp / (2 * V)
            # NN hop j+1 -> j : create j (right nbr occupied: beta), annihilate j+1 (left nbr occ: conj(betap))
            if exists(j + 1) and n(s, j + 1) == 1 and n(s, j) == 0 and n(s, j - 1) == 0:
                t = (s & ~(1 << jp)) | (1 << j)
                amp = -beta * np.conj(betap) / V
                K[pos[t], i] += amp
                K[i, pos[t]] += np.conj(amp)
            # dN = +1 at j
            if n(s, j - 1) == 0 and n(s, j) == 0 and n(s, j + 1) == 0:
                t = s | (1 << j)
                amp = 0j
                if exists(j + 1):   # pair (j, j+1), then annihilate j+1
                    amp += -p * np.conj(betap + beta * n(s, j + 2)) / (V * (1 + n(s, j + 2)))
                if exists(j - 1):   # pair (j-1, j), then annihilate j-1
                    amp += -p * np.conj(beta + betap * n(s, j - 2)) / (V * (1 + n(s, j - 2)))
                K[pos[t], i] += amp
                K[i, pos[t]] += np.conj(amp)
    return K, P


def partA(rng):
    print("Part A: exact second-order matrix vs closed form (w = lam = 0), random complex couplings")
    worst = 0.0
    for ring in (False, True):
        for L in range(8, 13):
            for trial in range(2 if L < 12 else 1):
                V = rng.uniform(3, 30)
                p, beta, betap = (rng.normal(size=3) + 1j * rng.normal(size=3))
                mu = rng.normal()
                Kx, P = second_order_exact(L, V, 0.0, p, beta, betap, mu, 0.0, ring)
                Kc, _ = closed_form(L, V, p, beta, betap, mu, ring)
                err = np.abs(Kx - Kc).max()
                worst = max(worst, err)
                print(f"  ring={ring!s:5} L={L:2d} V={V:6.2f} dimP={len(P):4d} max|K_exact - K_closed| = {err:.2e}"
                      f"   (max|K| = {np.abs(Kx).max():.2e})", flush=True)
    print(f"Part A worst deviation: {worst:.2e}")
    return worst


# ---------------- Part B: des Cloizeaux at finite V ----------------

def des_cloizeaux(L, V, w, p, beta, betap, mu, lam, ring):
    H = build_full(L, bond_matrix(V, w, p, beta, betap), site_matrix(mu, lam), ring)
    if not ring:
        # edge correction: give the edge sites the bulk 'missing-bond' environment is not possible exactly;
        # we only read bulk coefficients from the centre, so edges are left as they are.
        pass
    import scipy.linalg as sl
    e, U = sl.eigh(H, driver='evr')
    P = pspace(L, ring)
    m = len(P)
    A = U[P, :m]
    u_, sv, vh = np.linalg.svd(A)
    Uo = u_ @ vh
    He = Uo @ np.diag(e[:m]) @ Uo.conj().T
    return He, {s: i for i, s in enumerate(P)}, e[m] - e[m - 1], sv.min()


def predicted(V, w, p, beta, betap, mu, lam):
    fL, fR = lam + betap, lam + beta
    d, ssum = fL - fR, fL + fR
    pp = abs(p) ** 2
    return {
        "hop2": -pp / (2 * V),
        "V2": -pp / V + abs(d) ** 2 / (2 * V),
        "V3": -pp / (3 * V),
        "V4": 0.0,
        "hop1 (b+_{j+1} b_j)": w - fL * np.conj(fR) / V,
        "onsite+mu": 3 * pp / V - (abs(fL) ** 2 + abs(fR) ** 2) / V,
        "dN1 uniform": lam - p * np.conj(ssum) / V,
        "dN1 n_{j+2}": p * np.conj(d) / (2 * V),
        "dN1 n_{j-2}": -p * np.conj(d) / (2 * V),
        "dN2 sep2": 0.0,
        "dN2 sep3": 0.0,
    }


def measure(He, pos, L, ring):
    def idx(*sites):
        s = 0
        for j in sites:
            s |= 1 << (j % L)
        return pos[s]
    E = lambda *s: He[idx(*s), idx(*s)].real
    c = L // 2 - 1 if not ring else 2
    E0 = E()
    m = {}
    m["hop2"] = He[idx(c + 2), idx(c)]
    m["V2"] = E(c, c + 2) - E(c) - E(c + 2) + E0
    m["V3"] = E(c - 1, c + 2) - E(c - 1) - E(c + 2) + E0
    m["V4"] = E(c - 1, c + 3) - E(c - 1) - E(c + 3) + E0
    m["hop1 (b+_{j+1} b_j)"] = He[idx(c + 1), idx(c)]
    m["onsite+mu"] = E(c) - E0
    m["dN1 uniform"] = He[idx(c), idx()]
    # conditioned: create at c with a particle at c+2 (or c-2), minus the uniform part
    m["dN1 n_{j+2}"] = He[idx(c, c + 2), idx(c + 2)] - m["dN1 uniform"]
    m["dN1 n_{j-2}"] = He[idx(c - 2, c), idx(c - 2)] - m["dN1 uniform"]
    m["dN2 sep2"] = He[idx(c, c + 2), idx()]
    m["dN2 sep3"] = He[idx(c - 1, c + 2), idx()]
    return m


def partB(L=11, ring=False):
    print(f"\nPart B: des Cloizeaux effective Hamiltonian of the full chain, L = {L}, ring = {ring}")
    # a generic bond: p, beta, betap order eps = 1; w, lam, mu of order eps^2/V
    p0 = 1.0 * np.exp(0.3j)
    beta0, betap0 = 0.35 + 0.2j, -0.25 + 0.4j
    for V in (20.0, 40.0, 80.0, 160.0):
        t2 = abs(p0) ** 2 / (2 * V)
        w = (0.4 + 0.9j) * t2
        lam = (0.3 - 0.5j) * t2
        mu_eff = -2.2
        # mu in the bond is the bare chemical potential; onsite+mu absorbs the generated shift
        mu = mu_eff * t2 + 3 * abs(p0) ** 2 / V
        He, pos, gap, smin = des_cloizeaux(L, V, w, p0, beta0, betap0, mu, lam, ring)
        meas = measure(He, pos, L, ring)
        if not ring:
            pass
        meas["onsite+mu"] = meas["onsite+mu"] + mu
        pred = predicted(V, w, p0, beta0, betap0, mu, lam)
        print(f"  V = {V:6.1f}  t2 = {t2:.3e}  gap = {gap:.2f}  min sv = {smin:.4f}")
        for k in pred:
            pv, mv = pred[k], meas[k]
            dev = abs(mv - pv) / t2
            print(f"     {k:22s} meas = {complex(mv):.6e}  pred = {complex(pv):.6e}  |meas-pred|/t2 = {dev:.2e}")


if __name__ == "__main__":
    rng = np.random.default_rng(186)
    what = sys.argv[1] if len(sys.argv) > 1 else "AB"
    if "A" in what:
        partA(rng)
    if "B" in what:
        if "R" not in what:
            partB(11, False)
        partB(12, True)
