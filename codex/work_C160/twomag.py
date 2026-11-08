"""C160 (Claude Opus agent, 2026-11-08). Exact two-magnon problem near h_F (PROOF section 4).

Number-conserving part of a bond in the frame of the polarization axis n = z:
   h = Jp (XX + YY) + Jz ZZ + Dz (XY - YX) + Bz (Z1 + Z2)/2,   H = sum_j h_{j,j+1} on a ring.
Above the all-up state a magnon is a down spin. Prediction (PROOF eq. (4.2)-(4.4)):
   one magnon     e(k) = c0 + 2|a| cos(k - phi),  c0 = -4 Jz - 2 Bz,  |a| = 2 sqrt(Jp^2 + Dz^2), tan phi = Dz/Jp up to sign
   bound pair     E_b(K) = 2 c0 + V + tau^2 / V,  V = 4 Jz,  tau = 2|a| cos(K/2 - phi),  exists iff |V| > |tau|
   continuum      [2 c0 - 2|tau|, 2 c0 + 2|tau|]
Check by exact diagonalization of the two-magnon sector of an L-site ring, momentum-resolved.
"""
import numpy as np, itertools, sys

X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2)


def bond(Jp, Jz, Dz, Bz):
    return (Jp * (np.kron(X, X) + np.kron(Y, Y)) + Jz * np.kron(Z, Z) + Dz * (np.kron(X, Y) - np.kron(Y, X))
            + Bz * (np.kron(Z, I2) + np.kron(I2, Z)) / 2)


def two_magnon(h, L):
    basis = list(itertools.combinations(range(L), 2))
    idx = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    H = np.zeros((n, n), complex)
    E0 = L * h[0, 0].real
    for i, b in enumerate(basis):
        down = set(b)
        for j in range(L):
            k = (j + 1) % L
            s = (1 if j in down else 0) * 2 + (1 if k in down else 0)   # 0 = up up, ..., 3 = down down
            for s2 in range(4):
                amp = h[s2, s]
                if abs(amp) < 1e-15:
                    continue
                nd = set(down); nd.discard(j); nd.discard(k)
                if s2 & 2: nd.add(j)
                if s2 & 1: nd.add(k)
                if len(nd) != 2:
                    continue
                H[idx[tuple(sorted(nd))], i] += amp
    H -= E0 * np.eye(n)
    T = np.zeros((n, n))
    for i, b in enumerate(basis):
        T[idx[tuple(sorted(((b[0] + 1) % L, (b[1] + 1) % L)))], i] = 1
    return H, T, basis


def one_magnon(h, L):
    """single magnon dispersion from the bond (exact, number conserving)."""
    E0 = L * h[0, 0].real
    ks = 2 * np.pi * np.arange(L) / L
    # diagonal: magnon at j touches bonds (j-1,j) and (j,j+1)
    diag = (h[1, 1] + h[2, 2] - 2 * h[0, 0]).real
    hop = h[2, 1]          # |down at j, up at j+1> <- |up at j, down at j+1|: magnon hops j+1 -> j
    return ks, diag + hop * np.exp(1j * ks) + np.conj(hop) * np.exp(-1j * ks)


if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    cases = {'Heisenberg FM': (-1, -1, 0, 0), 'easy-axis eps=0.1': (-1, -1.1, 0, 0),
             'easy-axis + DM': (-1, -1.1, 0.15, 0.05), 'easy-plane eps=0.1': (-1.1, -1, 0, 0),
             'easy-plane + DM': (-1.1, -1, 0.2, 0)}
    for name, (Jp, Jz, Dz, Bz) in cases.items():
        h = bond(Jp, Jz, Dz, Bz)
        H, T, basis = two_magnon(h, L)
        Tp = [np.eye(len(basis))]
        for m in range(1, L):
            Tp.append(T @ Tp[-1])
        def sector(K, Hm, Tpow):
            P = sum(np.exp(-1j * K * m) * Tpow[m] for m in range(L)) / L
            U, sv, _ = np.linalg.svd(P)
            Q = U[:, sv > 0.5]
            return np.linalg.eigvalsh(Q.conj().T @ Hm @ Q)
        # one-magnon dispersion from the same translation convention (L x L problem)
        H1 = np.zeros((L, L), complex); T1 = np.zeros((L, L))
        E0 = L * h[0, 0].real
        for j in range(L):
            H1[j, j] = (h[1, 1] + h[2, 2] - 2 * h[0, 0]).real
            H1[j, (j + 1) % L] += h[2, 1]; H1[(j + 1) % L, j] += h[1, 2]
            T1[(j + 1) % L, j] = 1
        T1p = [np.eye(L)]
        for m in range(1, L):
            T1p.append(T1 @ T1p[-1])
        K1 = 2 * np.pi * np.arange(L) / L
        w1 = np.array([sector(K, H1, T1p)[0] for K in K1])
        a = np.sum(w1 * np.exp(1j * K1)) / L          # e1(K) = c0 + 2|a| cos(K - phi), phi = arg a
        c0 = w1.mean(); phi = np.angle(a)
        fit = np.max(np.abs(w1 - c0 - 2 * abs(a) * np.cos(K1 - phi)))
        V = 4 * Jz
        maxerr = 0.0; nb = 0
        for m in range(L):
            K = 2 * np.pi * m / L
            Es = sector(K, H, Tp)
            # pair energy: e1(k1)+e1(k2) with k1+k2=K; the relative-motion hopping is tau
            tau = 2 * abs(a) * np.cos(K / 2 - phi)
            if abs(V) > abs(tau) + 1e-9:
                Eb = 2 * c0 + V + tau ** 2 / V
                # the bound state is the lowest level below the continuum if V < 0
                cand = Es[0] if V < 0 else Es[-1]
                xi = 1 / np.log(abs(V / tau)) if abs(tau) > 1e-12 else 0.0   # bound-state size
                if xi < L / 8:
                    maxerr = max(maxerr, abs(cand - Eb)); nb += 1
        print(f'{name:22s} c0={c0:+.4f} |a|={abs(a):.4f} phi={phi:+.4f} 1-magnon fit {fit:.1e}  well-bound (xi < L/8) states at {nb}/{L} momenta, max |ED - formula| = {maxerr:.2e}')
