"""C194 exact pair-breaking gaps of the two controls, and the extrapolation fit applied to them.
Two-leg XX (blockade off, w = 0, t2 = -1, V = 0): two open chains of L/2 sites, hard-core bosons = free fermions,
  levels -2 cos(pi k/(l+1)); E(N) = min over N_a + N_b = N.
Hard rods (blockade on, w = 1, t2 = 0, V = 0): N rods on L sites <-> N free fermions on L - N + 1 sites
  (y_i = x_i - (i - 1)), levels -2 cos(pi k/(L - N + 2)).
Delta_pb = E(N+1) + E(N-1) - 2E(N), N even nearest n L.  Fit Delta(L) = D_inf + a/L + b/L^2 (registered form)."""
import numpy as np, sys
from fitgap import fit_gap

def ff(l, M):
    k = np.arange(1, M + 1); return float(np.sum(-2 * np.cos(np.pi * k / (l + 1)))) if M > 0 else 0.0

def E_xx(L, N):
    l = L // 2
    return min(ff(l, a) + ff(l, N - a) for a in range(max(0, N - l), min(N, l) + 1))

def E_hr(L, N):
    return ff(L - N + 1, N)

def dpb(E, L, N): return E(L, N + 1) + E(L, N - 1) - 2 * E(L, N)

if __name__ == '__main__':
    for name, E, ns in (('two-leg XX', E_xx, (1/32, 1/16, 1/8, 3/16, 1/4)), ('hard rods', E_hr, (1/32, 1/16, 1/8, 1/4))):
        for n in ns:
            Ls = [64, 128, 192, 256]
            Ns = [2 * int(round(n * L / 2)) for L in Ls]
            d = [dpb(E, L, N) for L, N in zip(Ls, Ns)]
            f = fit_gap(Ls, d)
            print(f"{name:10s} n={n:.4f} N={Ns} L*Dpb={[round(L*x,4) for L,x in zip(Ls,d)]} D_inf={f['D']:+.2e} +- {f['err']:.1e}  ({f['D']/f['err']:+.1f} sigma)  a={f['a']:.3f}")
