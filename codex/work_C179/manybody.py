"""Ring ED of (3.2) (hard-core bosons, periodic): parity gap and total-density Luttinger parameter K_rho at n = 1/6.
Delta_pb = E(N+1) + E(N-1) - 2E(N) (N even). K_rho (boson convention, Tonks = 1, two decoupled Tonks gases = 2):
  v K = pi L E''(Phi)  (Phi = total twist),   v/K = (L/pi) [E(N+2) + E(N-2) - 2E(N)] / 4,   K = sqrt(vK / (v/K)).
The two-pair vertex O_4 has dimension x_4 = 4/K_rho and is relevant iff K_rho > 2 (R4)."""
import numpy as np, sys, time
from ringv import e0, basis
for (L, N) in [tuple(int(x) for x in a.split(',')) for a in sys.argv[2:]]:
    for label, w in (("phi=pi/2", 1j * float(sys.argv[1])), ("phi=1.3", float(sys.argv[1]) * np.exp(1.3j))):
        t = time.time()
        E = {M: e0(L, M, w)[0] for M in (N - 2, N - 1, N, N + 1, N + 2)}
        dP = 0.05
        Ep = e0(L, N, w, Phi=dP)[0]; Em = e0(L, N, w, Phi=-dP)[0]
        curvP = (Ep + Em - 2 * E[N]) / dP ** 2
        vK = np.pi * L * curvP
        v_over_K = (L / np.pi) * (E[N + 2] + E[N - 2] - 2 * E[N]) / 4
        K = np.sqrt(vK / v_over_K) if vK > 0 and v_over_K > 0 else float('nan')
        dpb = E[N + 1] + E[N - 1] - 2 * E[N]
        print(f"L={L} N={N} {label} |w|={sys.argv[1]}: Delta_pb = {dpb:+.5f}, L*Delta_pb = {L*dpb:+.4f}; "
              f"vK = {vK:.4f}, v/K = {v_over_K:.4f}, K_rho = {K:.3f}, v = {np.sqrt(vK*v_over_K) if vK>0 and v_over_K>0 else float('nan'):.3f}, "
              f"mu window for N: [{E[N]-E[N-1]:.4f}, {E[N+1]-E[N]:.4f}]  ({time.time()-t:.0f}s)", flush=True)
