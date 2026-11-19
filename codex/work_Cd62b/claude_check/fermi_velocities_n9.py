"""Fermi velocities of the d = 9 chain at n = 9, from the repository's own data.

Uses contrib/markerless_d9_d10/certificate/full_loop_whitened_n9.npz and
full_loop_bands_n9.json, and the velocity formula of full_loop_integer.py
analyze(): v = |u^dag h'(k) u| / u^dag W u at each Fermi momentum k, with
h(k) = B + T e^{ik} + T^T e^{-ik} in whitened coordinates, B = mu + field,
and W = field * 2^field_power the image of the orthonormal site metric.

It also computes the eigenvalues of the Hermitian part of the hopping in the
orthonormal site basis (generalized eigenvalues of (sym T, W)), in the same
units, to compare the two ladders.

Run from the repository root:
  python3 codex/work_Cd62b/claude_check/fermi_velocities_n9.py
(needs numpy and scipy).
"""
import json, math, os, sys
import numpy as np
from scipy.linalg import eigh

root = sys.argv[1] if len(sys.argv) > 1 else "."
cdir = os.path.join(root, "contrib/markerless_d9_d10/certificate")
a = np.load(os.path.join(cdir, "full_loop_whitened_n9.npz"), allow_pickle=True)
bands = json.load(open(os.path.join(cdir, "full_loop_bands_n9.json")))

def cplx(x):
    return np.array([[complex(v) for v in row] for row in x], dtype=complex)

T = cplx(a["T"]); mu = cplx(a["mu"]); field = cplx(a["field"])
fp = int(a["field_power"])
W = (field * 2.0 ** fp).real
B = mu + field
n = T.shape[0]

# Hermitian part of the hopping in the orthonormal site basis, same units:
# x = R^T u, so h_orig is congruent to h_w with metric W; its eigenvalues are
# the generalized eigenvalues of (sym T_w, W). The generalized eigenvectors
# Y (with Y^dag W Y = 1) give the orthonormal-frame weights of a zero mode u:
# weight_a = |y_a^dag W u|^2 / (u^dag W u).
Sw = ((T + T.conj().T) / 2)
Sw = (Sw + Sw.conj().T) / 2
s, Y = eigh(Sw, W.astype(complex))
order = np.argsort(abs(s))[::-1]
s, Y = s[order], Y[:, order]
ls = np.log10(abs(s))

vel, slopes, deep2 = [], [], []
for k in bands["roots"]:
    h = B + T * np.exp(1j * k) + T.T * np.exp(-1j * k)
    w, U = np.linalg.eigh(h)
    u = U[:, np.argmin(abs(w))]
    hp = 1j * T * np.exp(1j * k) - 1j * T.T * np.exp(-1j * k)
    slope = abs(np.vdot(u, hp @ u).real)
    norm = np.vdot(u, W @ u).real
    slopes.append(slope)
    vel.append(slope / norm)
    wts = abs(Y.conj().T @ (W @ u)) ** 2 / norm
    deep2.append(float(wts[-2:].sum()))
vel = np.array(vel)
lv = np.log10(vel)

out = dict(
    n=n,
    fermi_points=len(vel),
    log10_velocity_max=float(lv.max()),
    log10_velocity_min=float(lv.min()),
    log10_velocity_span=float(lv.max() - lv.min()),
    certificate_log10_min_speed_normalized=bands["log10_min_speed_normalized"],
    log10_velocities_sorted=sorted(map(float, lv), reverse=True),
    log10_abs_sym_hopping_eigenvalues=list(map(float, ls)),
    log10_hopping_span=float(ls[0] - ls[-1]),
    whitened_slope_min=float(min(slopes)),
    whitened_slope_max=float(max(slopes)),
    zero_mode_weight_on_two_deepest_orbitals_min=min(deep2),
)
assert abs(out["log10_velocity_min"] - out["certificate_log10_min_speed_normalized"]) < 1e-6
assert ls[-1] < out["log10_velocity_min"] and out["log10_velocity_max"] < ls[-2]
print(json.dumps(out, indent=1))
