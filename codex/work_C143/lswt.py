"""Linear spin-wave (Holstein-Primakoff) data of a bond around its uniform classical minimum n:
transverse block K = E^T J E (E = two unit vectors orthogonal to n), longitudinal energy, and the spin-wave
dispersion in the free-fermion form of PROOF (2.1): f(k) = -a sin k + sqrt((c cos k + Beff)^2 + g^2 sin^2 k)."""
import json, sys
import numpy as np
from scipy.optimize import minimize
from ed import params_to_JB
from classical_landscape import energy
fn, idx = sys.argv[1], int(sys.argv[2])
d = [json.loads(l) for l in open(fn)][idx]
J, B = params_to_JB(np.array(d["x"]))
rng = np.random.default_rng(1)
best = min((minimize(energy, rng.uniform(0, 2 * np.pi, 2), args=(J, B, 1), method="BFGS") for _ in range(50)), key=lambda r: r.fun)
th, ph = best.x; n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
e1 = np.cross(n, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
E = np.stack([e1, e2], 1)
K = E.T @ J @ E
Jnn = n @ J @ n; Bn = B @ n
# spin-wave: flipping costs -(2 Jnn + Bn)*... in Pauli units the classical local field on a site is h_loc = B + (J + J^T) n
hloc = B + (J + J.T) @ n
a = K[0, 1] - K[1, 0]; c = K[0, 0] + K[1, 1]; g = abs((K[0, 0] - K[1, 1]) - 1j * (K[0, 1] + K[1, 0]))
Beff = -(hloc @ n)   # spin anti-aligned with the local field: excitation cost |h_loc|, written with the sign of PROOF (2.1)
print("n =", np.round(n, 4), " |h_loc| =", round(float(np.linalg.norm(hloc)), 4), " h_loc parallel to n:", round(float(abs(hloc @ n) / np.linalg.norm(hloc)), 6))
print("transverse block K =", np.round(K, 4).tolist(), " det K =", round(float(np.linalg.det(K)), 5))
print(f"a={a:.4f} c={c:.4f} g={g:.4f} Beff={Beff:.4f}")
ks = np.linspace(-np.pi, np.pi, 2001)
f = -a * np.sin(ks) + np.sqrt((c * np.cos(ks) + Beff) ** 2 + g ** 2 * np.sin(ks) ** 2)
print("spin-wave f(k): min", round(float(f.min()), 4), "max", round(float(f.max()), 4), "(bandwidth relative to |h_loc|:", round(float((f.max() - f.min()) / np.linalg.norm(hloc)), 4), ")")
