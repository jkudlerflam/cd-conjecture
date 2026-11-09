"""C167: finite-entanglement exponent of the staggered moment, m_s(chi) ~ xi(chi)^(-x_pi), and c_pair, along the
hi_K_1 path.  x_pi from successive chi pairs (two-site VUMPS, fesx.py output).  In a half-floating liquid x_pi is
the dimension of the commensurate vertex V_c at momentum pi; the liquid is stable iff x_pi > 1/2 (C159 table 4.1,
V_c^2 of dimension 4 x_pi), and an independent Ising mode at its critical point would give x_pi = 1/8.
usage: python3 xpi_fes.py loc/path/*.json"""
import sys, json, re, numpy as np
rows = []
for fn in sys.argv[1:]:
    lam = float(re.findall(r'_([0-9.]+)\.json', fn)[0])
    d = json.load(open(fn))['results']
    ms = [np.linalg.norm((np.array(r['mag'])[:, 0] - np.array(r['mag'])[:, 1]) / 2) for r in d]
    xi = [r['xi'] for r in d]; S = [np.mean(r['S']) for r in d]
    xs = [-np.log(ms[k + 1] / ms[k]) / np.log(xi[k + 1] / xi[k]) for k in range(len(d) - 1) if ms[k] > 1e-6]
    cs = [6 * (S[k + 1] - S[k]) / np.log(xi[k + 1] / xi[k]) for k in range(len(d) - 1)]
    rows.append((lam, fn.split('/')[-1], d[-1]['chi'], d[-1]['e'], ms[-1], xi[-1], xs[-2:], cs[-2:]))
for r in sorted(rows):
    print(f"{r[1]:22s} lam={r[0]:.3f} chi_max={r[2]:3d} e={r[3]:.9f} m_s={r[4]:.4f} xi={r[5]:7.2f} "
          f"x_pi(last pairs)={np.round(r[6], 3).tolist()} c_pair(last)={np.round(r[7], 2).tolist()}")
