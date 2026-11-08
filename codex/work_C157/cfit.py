"""Central-charge fits of open-chain S(l) from a dmrg_run jsonl: parity-resolved fits of S = (c/6) log chord + const
over windows l in [w L, (1-w) L] for w = 0.05, 0.125, 0.25."""
import sys, json
import numpy as np
def fits(S, L):
    S = np.asarray(S); l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    out = {}
    for w in (0.05, 0.125, 0.25):
        m = (l >= w * L) & (l <= (1 - w) * L)
        out[w] = [round(float(6 * np.polyfit(x[m & (l % 2 == p)], S[m & (l % 2 == p)], 1)[0]), 3) for p in (0, 1)]
    return out
if __name__ == "__main__":
    for f in sys.argv[1:]:
        for line in open(f):
            d = json.loads(line)
            print(f"{d['label']} L={d['L']} alpha={d['alpha']} beta={d['beta']} {d['extra']}: b_mid={d['b_mid']:+.3f} b_stag={d['b_stag']:+.3f} a_stag={d['a_stag']:+.3f} S_mid={d['S_mid']:.3f} c(w)={fits(d['S'], d['L'])}")


def avg_fit(S, L, w=0.125):
    """Fit of the parity-averaged entropy (S(l) + S(l+1))/2 against the chord at l + 1/2."""
    S = np.asarray(S); l = np.arange(1, L - 1) + 0.5
    Sa = 0.5 * (S[:-1] + S[1:])
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= w * L) & (l <= (1 - w) * L)
    return round(float(6 * np.polyfit(x[m], Sa[m], 1)[0]), 3)
