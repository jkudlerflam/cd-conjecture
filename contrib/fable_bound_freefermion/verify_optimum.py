"""Re-find the m=3 optimum and verify A by a brute-force band scan (sign changes on a fine grid,
slopes by finite differences), independent of the companion-matrix root finder."""
import numpy as np
from scipy.optimize import minimize
from ff_bound import casimir_A, trace_norm, bloch_nc, bloch_mj
from params import unpack_nc, unpack_mj, neg_ratio

def brute_A(kind, A0, B, N=200001):
    h = bloch_nc(A0, B) if kind == "nc" else bloch_mj(A0, B)
    ks = np.linspace(-np.pi, np.pi, N, endpoint=False)
    # continuity-sorted bands are not needed: count sign changes of each sorted eigenvalue branch;
    # a crossing of two bands through zero at the same k is counted by both sorted branches.
    E = np.array([np.linalg.eigvalsh(h(k)) for k in ks])
    dk = ks[1] - ks[0]
    tot = 0.0
    for j in range(E.shape[1]):
        e = E[:, j]
        s = np.sign(e)
        idx = np.where(s * np.roll(s, -1) < 0)[0]
        for i in idx:
            i2 = (i + 1) % N
            tot += abs(e[i2] - e[i]) / dk
    return tot / (2.0 if kind == "nc" else 4.0)

if __name__ == "__main__":
  rng = np.random.default_rng(11)
  for kind, n, trials in (("nc", 3, 40), ("mj", 3, 20), ("mj", 6, 30), ("nc", 2, 20), ("mj", 4, 20)):
      best, bestx = 0, None
      dim = 4 * n * n if kind == "nc" else 2 * n * n
      for _ in range(trials):
          r = minimize(neg_ratio, rng.standard_normal(dim), args=(kind, n), method="Nelder-Mead",
                       options={"maxiter": 6000, "xatol": 1e-9, "fatol": 1e-11})
          if -r.fun > best:
              best, bestx = -r.fun, r.x
      A0, B = (unpack_nc if kind == "nc" else unpack_mj)(bestx, n)
      norm = trace_norm(B) if kind == "nc" else trace_norm(B) / 2
      Ab = brute_A(kind, A0, B)
      print(f"{kind} n={n}: optimizer ratio {best:.5f}; brute-force A/norm = {Ab/norm:.5f}")
      sv = np.linalg.svd(B, compute_uv=False)
      print("   singular values/max:", np.round(sv / sv.max(), 4), " ||A0||/sigma_max:", round(np.linalg.norm(A0, 2) / sv.max(), 4))
      np.savez(f"opt_{kind}_{n}.npz", A0=A0, B=B)
