"""Level statistics of the full qubit ring for the candidate bond family, resolved by momentum K and Z-parity.
Bond: h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z1 + 1Z) on every bond of a ring of L sites.
Family: Jx = -Jy = Delta/2 = 1/2, Jz = V/4, D = |w| t2 / 2 with t2 = Delta^2/(2V) (|w| = 1, phi = pi/2),
hz from mu_eff = -2.2 t2. Sectors K = 0, pi, +-pi/2 are skipped (U and inversion-type symmetries act there).
Control: D = 0 (C153's hidden staggered U(1), Q = sum (-1)^j Z_j, not resolved here) should depress <r>."""
import numpy as np, scipy.sparse as sp, sys
from ringk import rstat
L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
dim = 1 << L
s = np.arange(dim, dtype=np.int64)
bits = lambda j: (s >> (j % L)) & 1   # 1 = spin down = particle
def zval(j): return 1 - 2 * bits(j)
def build(Jx, Jy, D, Jz, hz):
    diag = np.zeros(dim)
    rows = [s]; cols = [s]; vals = None
    off_r, off_c, off_v = [], [], []
    for j in range(L):
        k = (j + 1) % L
        diag += Jz * zval(j) * zval(k) + hz * zval(j)
        # XX, YY, XY, YX all flip both spins j, k. Matrix elements <s'|..|s>, s' = s ^ (bit j | bit k)
        t = s ^ ((1 << j) | (1 << k))
        zj, zk = zval(j), zval(k)   # Z eigenvalues of the initial state
        # X|z> = |-z>, Y|z> = i z |-z>  (Y|up> = i|down>, Y|down> = -i|up>; z=+1 up)
        amp = Jx * 1 + Jy * (1j * zj) * (1j * zk) + D * ((1j * zk) - (1j * zj))
        off_r.append(t); off_c.append(s); off_v.append(amp)
    H = sp.csr_matrix((np.concatenate([diag.astype(complex)] + off_v), (np.concatenate([s] + off_r), np.concatenate([s] + off_c))), shape=(dim, dim))
    return H
def rotl(x, d):
    d %= L; m = dim - 1
    return ((x << d) | (x >> (L - d))) & m
allrot = np.stack([rotl(s, d) for d in range(L)])
rep = allrot.min(axis=0); dmin = allrot.argmin(axis=0); shift = (-dmin) % L
par = np.array([bin(x).count('1') & 1 for x in range(dim)])
ureps = np.unique(rep)
per_of = {}
for r in ureps:
    per_of[int(r)] = next(p for p in range(1, L + 1) if rotl(r, p) == r)
def sector(H, m, p):
    K = 2 * np.pi * m / L
    R = [int(r) for r in ureps if par[r] == p and abs(np.exp(-1j * K * per_of[int(r)]) - 1) < 1e-9]
    col = {r: i for i, r in enumerate(R)}
    sel = np.array([int(r) in col for r in rep])
    idx = np.nonzero(sel)[0]
    c = np.array([col[int(rep[i])] for i in idx])
    per = np.array([per_of[int(rep[i])] for i in idx])
    # |K,r> = sqrt(per)/L sum_d e^{-iKd} T^d |r>; config x = T^{shift} rep appears L/per times with phase e^{-iK shift}
    v = np.exp(-1j * K * shift[idx]) * np.sqrt(per) / L * (L / per)
    B = sp.csr_matrix((v, (idx, c)), shape=(dim, len(R)))
    return (B.conj().T @ H @ B).toarray()
# translation convention check: T^d maps bit j -> bit j+d, consistent with rotl.
for VD in [float(x) for x in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['4', '8'])]:
    Delta = 1.0; V = VD * Delta; t2 = Delta ** 2 / (2 * V)
    Jz = V / 4; mu = -2.2 * t2 + 3 * Delta ** 2 / V; hz = (mu - 4 * Jz) / 2
    for label, D in (("candidate family D = t2/2", t2 / 2), ("control D = 0 (C153 class)", 0.0)):
        H = build(0.5, -0.5, D, Jz, hz)
        rs = []
        for m in (1, 2, 3, 5, 6, 7):
            for p in (0, 1):
                A = sector(H, m, p)
                assert np.abs(A - A.conj().T).max() < 1e-10
                e = np.linalg.eigvalsh(A)
                r, nr, nd = rstat(e)
                rs.append(r)
        print(f"L={L} V/Delta={VD:g} {label:28s}: <r> over 12 sectors = {np.mean(rs):.4f} +- {np.std(rs)/np.sqrt(len(rs)):.4f} (sector dim ~{len(e)})", flush=True)
