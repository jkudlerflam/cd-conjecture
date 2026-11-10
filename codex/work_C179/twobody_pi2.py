"""Two-body problem of (3.2) at phi = pi/2 on a ring of L sites (C173's twobody.py was run at phi = 1.3 only).
Relative coordinate r = x2 - x1 in [2, L-2] (blockade), |K, r> = sum_x e^{iKx} |x, x + r>.
Particle 2 hopping by d: r -> r + d, amplitude a_d. Particle 1 hopping by d: r -> r - d, amplitude a_d e^{-iKd}.
a_{+1} = w, a_{-1} = w*, a_{+-2} = t2 = -1; V2 = -2 at r = 2, L-2; V3 = -2/3 at r = 3, L-3.
Bound state declared if E_2(K) < min_k [eps(k) + eps(K-k)] - 1e-3, eps(k) = 2|w| sin k - 2 cos 2k (pre-registered R3)."""
import numpy as np, sys
def H2(K, w, L, t2=-1.0, V2=-2.0, V3=-2.0/3.0):
    rs = np.arange(2, L - 1); n = len(rs); H = np.zeros((n, n), complex)
    amp = {1: w, -1: np.conj(w), 2: t2, -2: t2}
    for i, r in enumerate(rs):
        H[i, i] += V2 * ((r == 2) + (r == L - 2)) + V3 * ((r == 3) + (r == L - 3))
        for d, a in amp.items():
            if 2 <= r + d <= L - 2: H[r + d - 2, i] += a
            if 2 <= r - d <= L - 2: H[r - d - 2, i] += a * np.exp(-1j * K * d)
    assert np.abs(H - H.conj().T).max() < 1e-12
    return np.linalg.eigvalsh(H)[0]
ks = np.linspace(-np.pi, np.pi, 40001)
for wa in (0.5, 1.0, 2.0):
    eps = lambda k: 2 * wa * np.sin(k) - 2 * np.cos(2 * k)
    emin = eps(ks).min()
    for L in (400, 800):
        rows = []
        for K in np.linspace(-np.pi, np.pi, 97):
            rows.append((K, H2(K, 1j * wa, L), (eps(ks) + eps(K - ks)).min()))
        d = np.array([r[1] - r[2] for r in rows]); Ks = np.array([r[0] for r in rows])
        g = min(r[1] for r in rows)
        Kb = Ks[d < -1e-3]
        print(f"|w|={wa} L={L}: 2eps_min={2*emin:.5f}, min_K E_2 - 2eps_min = {g-2*emin:+.5f}; "
              f"min_K [E_2(K) - continuum(K)] = {d.min():+.5f} at K/pi = {Ks[d.argmin()]/np.pi:+.3f}; bound at {len(Kb)}/97 K"
              + (f", K/pi in [{Kb.min()/np.pi:+.3f}, {Kb.max()/np.pi:+.3f}]" if len(Kb) else ""), flush=True)
        if L == 800:
            for K in (0.0, np.pi / 2, np.pi, -np.pi / 2):
                c = (eps(ks) + eps(K - ks)).min(); e = H2(K, 1j * wa, L)
                print(f"     K/pi={K/np.pi:+.2f}: E_2 - continuum = {e-c:+.5f}")
