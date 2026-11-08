"""Translation-invariant Hartree-Fock (Wick) for interacting Majorana chains
    H = t sum_j b_j + sum_j sum_{r>=2} U_r b_j b_{j+r},   b_j = -i g_j g_{j+1},
in the thermodynamic limit, with G_r = <i g_j g_{j+r}> (G_{-r} = -G_r, G_0 = 0).
Wick: <b_j> = -G_1,  <b_j b_{j+r}> = G_1^2 - G_r^2 + G_{r+1} G_{r-1}.
Energy per site e(G) = -t G_1 + sum_r U_r (G_1^2 - G_r^2 + G_{r+1} G_{r-1}).
Mean-field Hamiltonian sum_j sum_{r>0} h_r i g_j g_{j+r} with h_r = de/dG_r, single-particle
eps(k) = 4 sum_r h_r sin(k r), ground state G_r = -(1/pi) int_0^pi sign(eps(k)) sin(k r) dk.
(Check: h_1 = t > 0 only gives G_1 = -2/pi, G_3 = -2/(3 pi).)
Zeros: 2 automatic (k = 0, pi) + 2 x (sign changes of eps on (0, pi)); c_MF = zeros/4.
This ansatz excludes translation breaking (dimerization); it tests only the Lifshitz route to c = 3/2.
"""
import numpy as np

R = 200
k = (np.arange(8000) + 0.5) * np.pi / 8000
SIN = np.sin(np.outer(np.arange(R + 3), k))  # SIN[r] = sin(k r)


def corr(h):
    eps = 4 * (h[1:, None] * SIN[1:R + 1]).sum(0)
    s = np.sign(eps)
    G = np.zeros(R + 3)
    G[1:R + 2] = -(s[None, :] * SIN[1:R + 2]).mean(1)  # (1/pi) int dk = mean over grid
    nz = int(np.sum(s[1:] != s[:-1]))
    return G, eps, nz


def grad(G, t, U):
    h = np.zeros(R + 1)
    h[1] = -t
    for r, u in U.items():
        h[1] += 2 * u * G[1]
        h[r] += -2 * u * G[r]
        h[r + 1] += u * G[r - 1]
        if r - 1 >= 1:
            h[r - 1] += u * G[r + 1]
    return h


def energy(G, t, U):
    return -t * G[1] + sum(u * (G[1] ** 2 - G[r] ** 2 + G[r + 1] * G[r - 1]) for r, u in U.items())


def solve(t, U, iters=1500, mix=0.4, tol=1e-12):
    h = np.zeros(R + 1); h[1] = -t if t != 0 else -1.0
    G, eps, nz = corr(h)
    for it in range(iters):
        h = grad(G, t, U)
        Gn, eps, nz = corr(h)
        err = np.abs(Gn - G).max()
        G = (1 - mix) * G + mix * Gn
        if err < tol:
            break
    h = grad(G, t, U)
    G2, eps, nz = corr(h)
    return dict(e=energy(G, t, U), G=G, h=h, zeros=2 + 2 * nz, c=(2 + 2 * nz) / 4, err=err, it=it + 1)


if __name__ == "__main__":
    print("RZFA calibration: H = sum (i t g g + g gggg) = -t sum b - g sum b_j b_{j+2}  (t = 1): t_B = -1, U_2 = -g")
    for g in (1.0, 0.3, 0.0, -0.2, -0.26, -0.28, -0.29, -0.30, -0.32, -0.5, -1.0, -2.0):
        r = solve(-1.0, {2: -g})
        h = r["h"]
        print(f"  g/t = {g:+.2f}: G1 = {r['G'][1]:+.4f} G3 = {r['G'][3]:+.4f}  h1 = {h[1]:+.4f} h3 = {h[3]:+.4f} "
              f"v ~ h1 + 3 h3 = {h[1] + 3*h[3]:+.4f}  zeros = {r['zeros']}  c_MF = {r['c']}  (it {r['it']}, err {r['err']:.0e})")
