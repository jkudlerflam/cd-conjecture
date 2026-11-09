#!/usr/bin/env python3
"""C177: entanglement-profile estimators robust to incommensurate oscillations.
c_osc: least squares S(l) = (c/6) X(l) + s0 + sum_q [A_q cos(q l) + B_q sin(q l)] X'(l)  over L/8 <= l <= 7L/8, with
X = log[(2L/pi) sin(pi l/L)], oscillation envelope X'(l) = 1/sin(pi l/L)^p, p = 1 fixed (heuristic), and one or two
oscillation wavevectors chosen on a grid by least residual.  Sw: mean S over the central quarter (window average).
usage: sfit.py file.json ..."""
import sys, json, numpy as np

def fit(S, L, nq=2, p=1.0):
    l = np.arange(1, L); X = np.log(2 * L / np.pi * np.sin(np.pi * l / L))
    w = (l >= L // 8) & (l <= L - L // 8)
    env = 1 / np.sin(np.pi * l / L) ** p
    qs = np.linspace(0.05, np.pi, 316)
    def solve(qsel):
        cols = [X[w] / 6, np.ones(w.sum())]
        for q in qsel:
            cols += [np.cos(q * l[w]) * env[w], np.sin(q * l[w]) * env[w]]
        A = np.stack(cols, 1)
        c, res, *_ = np.linalg.lstsq(A, S[w], rcond=None)
        r = float(np.sum((A @ c - S[w]) ** 2))
        return c[0], r
    best = (None, np.inf, None)
    q1s = []
    for q in qs:
        c, r = solve([q])
        if r < best[1]: best = (c, r, [q])
    if nq == 2:
        b2 = best
        for q in qs:
            if abs(q - best[2][0]) < 0.05: continue
            c, r = solve([best[2][0], q])
            if r < b2[1]: b2 = (c, r, [best[2][0], q])
        best = b2
    return best

if __name__ == '__main__':
    for f in sys.argv[1:]:
        d = json.load(open(f)); S = np.array(d['S']); L = d['L']
        c1, r1, q1 = fit(S, L, 1); c2, r2, q2 = fit(S, L, 2)
        Sw = float(np.mean(S[3 * L // 8 - 1:5 * L // 8]))
        print(f"{d['tag']:6s} L={L:4d} Sbar={d['Sbar']:.4f} Sw={Sw:.4f} c_prof(pi)={d['c_prof']:.3f} c_osc1={c1:.3f} (q={q1[0]:.3f}) "
              f"c_osc2={c2:.3f} (q={np.round(q2, 3).tolist()}) rms={np.sqrt(r2 / (L * 3 / 4)):.1e} trunc={d['trunc']:.1e}")
