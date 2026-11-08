"""Size-scaling central-charge estimators robust to incommensurate oscillations (open chains).
For each record: Sbar = mean of S(l) over the central window l/L in [3/8, 5/8] (cuts as stored), and the parity-averaged
profile fit (C144 cfit.avg_fit) for reference. Between sizes L1 < L2 of the same label: c_pair = 6 (Sbar2 - Sbar1) / log(L2/L1),
the open-chain Calabrese-Cardy slope (S = (c/6) log L + const at fixed l/L).
Usage: python3 cscale.py key file1.jsonl file2.jsonl ...   (key: json field used to group, e.g. label or g).
Records may store S per cut of L-1 sites (qubit chains) or per block cut (rzfa_dmrg.py, L = spins, cuts at 2m)."""
import sys, json
import numpy as np
key = sys.argv[1]
recs = {}
for f in sys.argv[2:]:
    for line in open(f):
        d = json.loads(line)
        S = np.asarray(d["S"]); L = d["L"]
        if len(S) == L - 1:
            l = np.arange(1, L)
        else:
            l = 2 * np.arange(1, len(S) + 1)
        m = (l >= 3 * L / 8) & (l <= 5 * L / 8)
        Sbar = float(S[m].mean())
        x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
        w = (l >= L / 8) & (l <= 7 * L / 8)
        # smoothed fit: running mean over 8 consecutive stored cuts, then linear fit in the chord variable
        k = 8 if len(S) == L - 1 else 4
        Ss = np.convolve(S, np.ones(k) / k, mode="valid"); xs = np.convolve(x, np.ones(k) / k, mode="valid")
        ls = np.convolve(l, np.ones(k) / k, mode="valid")
        ws = (ls >= L / 8) & (ls <= 7 * L / 8)
        csm = float(6 * np.polyfit(xs[ws], Ss[ws], 1)[0])
        recs.setdefault(str(d[key]), []).append((L, Sbar, csm, d.get("trunc", d.get("hist", [[0, 0, 0]])[-1][2])))
for k, v in recs.items():
    v.sort()
    print(f"{key} = {k}")
    for i, (L, Sb, cs, tr) in enumerate(v):
        s = f"   L = {L:4d}: Sbar = {Sb:.4f}  smoothed profile fit c = {cs:+.3f}  trunc = {tr:.1e}"
        if i > 0:
            L0, S0 = v[i - 1][0], v[i - 1][1]
            s += f"   c_pair({L0},{L}) = {6 * (Sb - S0) / np.log(L / L0):+.3f}"
        print(s)
