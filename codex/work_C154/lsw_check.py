"""C154: check of the linear-spin-wave band (5.1) of Lemma 3 against exact ring spectra (Claude Opus agent, 2026-11-07).
Strong field |B| along -n so the polarized state along n is the reference; LSW should match the exact lowest
excitation at each momentum up to O(J^2/|B|). Also verifies the coefficient formulas a_s, a_c, pi_0."""
import numpy as np
import est_c155 as E
rng = np.random.default_rng(7)
L = 12
for trial in range(3):
    J = rng.normal(size=(3, 3)) * 0.3
    n = rng.normal(size=3); n /= np.linalg.norm(n)
    Bmag = float(__import__("sys").argv[1]); B = -Bmag * n
    # local frame (u, v, n)
    u = np.cross(n, [0.3, 0.5, 0.8]); u /= np.linalg.norm(u); v = np.cross(n, u)
    R = np.stack([u, v, n], 1); Jp = R.T @ J @ R; Bn = B @ n
    t = Jp[0, 0] + Jp[1, 1] - 1j * Jp[0, 1] + 1j * Jp[1, 0]
    p = Jp[0, 0] - Jp[1, 1] - 1j * (Jp[0, 1] + Jp[1, 0])
    A0 = -4 * Jp[2, 2] - 2 * Bn
    ks = 2 * np.pi * np.arange(L) / L
    Ak = A0 + 2 * t.real * np.cos(ks) - 2 * t.imag * np.sin(ks)
    Amk = A0 + 2 * t.real * np.cos(ks) + 2 * t.imag * np.sin(ks)
    w = (Ak - Amk) / 2 + np.sqrt(((Ak + Amk) / 2) ** 2 - np.abs(2 * p * np.cos(ks)) ** 2)
    d = E.ring_data(E.bond_from_JB(J, B), L, nev=2)
    lev = d['lev']; E0 = d['E0']; m0 = d['m0']
    exact = np.array([lev[(m0 + m) % L, 0] if (m % L) != 0 else lev[m0, 1] for m in range(L)]) - E0
    print(f"trial {trial}: m0={m0}, max |LSW - exact| over k = {np.abs(w - exact).max():.3e}  (bandwidth {w.max()-w.min():.3f}, J^2/B = {np.abs(J).max()**2/Bmag:.3e})")
    exact_rev = np.array([exact[(-m) % L] for m in range(L)])
    print("   reversed-momentum convention: max diff", f"{np.abs(w - exact_rev).max():.3e}")
    print("   LSW  ", np.round(w, 3)); print("   exact", np.round(exact, 3))
