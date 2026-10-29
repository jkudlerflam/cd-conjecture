"""Adversarial search for sup A/||B||_1 (NC) and sup A/(||B||_1/2) (MJ), with every candidate
re-checked by the brute-force band scan. Usage: python3 search_ratio.py kind n trials seed"""
import sys, numpy as np
from scipy.optimize import minimize
from ff_bound import casimir_A, trace_norm
from params import unpack_nc, unpack_mj
from verify_optimum import brute_A

kind, n, trials, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
unpack = unpack_nc if kind == "nc" else unpack_mj
dim = 4 * n * n if kind == "nc" else 2 * n * n

def nrm(B):
    return trace_norm(B) if kind == "nc" else trace_norm(B) / 2

def neg(x):
    try:
        A0, B = unpack(x, n)
        return -casimir_A(kind, A0, B) / nrm(B)
    except Exception:
        return 0.0

best = 0.0
for t in range(trials):
    r = minimize(neg, rng.standard_normal(dim), method="Nelder-Mead",
                 options={"maxiter": 5000, "xatol": 1e-9, "fatol": 1e-11})
    A0, B = unpack(r.x, n)
    rb = brute_A(kind, A0, B, N=20001) / nrm(B)
    if rb > best + 1e-9:
        best = rb
        np.savez(f"best_{kind}_{n}_seed{seed}.npz", A0=A0, B=B)
        sv = np.linalg.svd(B, compute_uv=False)
        print(f"trial {t}: optimizer {-r.fun:.5f}  brute {rb:.5f}  sv/max {np.round(sv/sv.max(),3)}  |A0|/smax {np.linalg.norm(A0,2)/sv.max():.3f}", flush=True)
print(f"FINAL {kind} n={n}: best brute-force ratio = {best:.5f}", flush=True)
