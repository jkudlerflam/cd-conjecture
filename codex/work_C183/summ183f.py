#!/usr/bin/env python3
"""C183 stage F summary (Claude Opus 5.5 agent, 2026-11-09), written before any output existed.
usage: summ183f.py DIR   (DIR holds oc/ and ctrl/ of job 17604545)
Per point S1, S2: S-bar, c_prof, truncation, n_mid at L = 64, 128, 256; c2 = 6 dS-bar / dln L for 64 -> 128 and 128 -> 256;
the parity gap L * (E_odd - E_even) at 64 and 128 (and N_odd - N_even); n(k) proxy from the one-body row
G(r) = <s+_{i0} s-_{i0+r}> (i0 = L/4, r = 1 .. L/2): n(k) ~ Re sum_r G(r) e^{-ikr} (Hann window), with the weights on
the two valleys (k in (0, pi/2) and (pi/2, pi), the valleys sitting at k_a and pi - k_a on the stratum).
Controls: c2 of XX and the two-leg chain.  Evaluates F1 / F2 and P1 / P4 as registered in RESULT.md."""
import sys, os, json, math
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else '.'


def load(p):
    try:
        return json.load(open(os.path.join(D, p)))
    except Exception:
        return None


def nk(d):
    G = np.array([complex(*z) for z in d['G']])
    m = min(len(G), d['L'] // 2)
    G = G[:m]
    r = np.arange(1, m + 1)
    win = 0.5 * (1 + np.cos(np.pi * r / (m + 1)))
    k = np.linspace(-np.pi, np.pi, 721)
    n = np.array([2 * np.real(np.sum(win * G * np.exp(-1j * kk * r))) for kk in k])
    q = [float(n[(k > lo) & (k < lo + np.pi / 2)].clip(min=0).sum() * (k[1] - k[0])) for lo in (-np.pi, -np.pi / 2, 0, np.pi / 2)]
    return float(k[np.argmax(n)]), q


for kind in ('xx', 'leg2'):
    S = {L: load(f'ctrl/{kind}_{L}.json') for L in (64, 128, 256)}
    line = f"control {kind}: " + ', '.join(f"L={L} Sbar={d['Sbar']:.4f} cprof={d['c_prof']:.3f}" for L, d in S.items() if d)
    for L1, L2 in ((64, 128), (128, 256)):
        if S[L1] and S[L2]:
            line += f"; c2({L1}->{L2}) = {6 * (S[L2]['Sbar'] - S[L1]['Sbar']) / math.log(2):.3f}"
    print(line)
for tag in ('S1', 'S2'):
    S = {L: load(f'oc/{tag}_{L}.json') for L in (64, 128, 256)}
    for L, d in S.items():
        if d:
            kpk, q = nk(d) if 'G' in d else (np.nan, [np.nan] * 4)
            print(f"{tag} L={L}: Sbar={d['Sbar']:.4f} c_prof={d['c_prof']:.3f} trunc={d['trunc']:.1e} n_mid={d['n_mid']:.4f} "
                  f"E={d['E']:.10f} chi={d['chis'][-1]} n(k) peak k={kpk:+.3f}, quadrant weights (-pi..pi) {np.round(q, 3).tolist()}")
    c2 = {}
    for L1, L2 in ((64, 128), (128, 256)):
        if S[L1] and S[L2]:
            c2[(L1, L2)] = 6 * (S[L2]['Sbar'] - S[L1]['Sbar']) / math.log(2)
            print(f"{tag} c2({L1}->{L2}) = {c2[(L1, L2)]:.3f}")
    g = {}
    for L in (64, 128):
        o, e = load(f'oc/{tag}odd_{L}.json'), load(f'oc/{tag}even_{L}.json')
        if o and e:
            g[L] = L * (o['E'] - e['E'])
            print(f"{tag} L={L}: E_odd - E_even = {o['E'] - e['E']:+.3e}, L*Dp = {g[L]:+.4f}, N_odd - N_even = {o['N'] - e['N']:+.3f}")
    if len(c2) == 2 and len(g) == 2:
        F1 = min(c2.values()) >= 1.6 and abs(g[128]) <= 1.2 * abs(g[64])
        F2 = c2[(128, 256)] < 1.4 and abs(g[128]) > abs(g[64])
        P1 = (c2[(128, 256)] > 1.2 and S[128]['c_prof'] > 1.2 and S[256]['c_prof'] > 1.2
              and S[256]['c_prof'] >= S[128]['c_prof'] - 0.03 and S[128]['trunc'] < 1e-7 and S[256]['trunc'] < 1e-7)
        P4 = abs(g[128]) < 1.5 * abs(g[64])
        print(f"{tag}: F1 {'yes' if F1 else 'no'}, F2 {'yes' if F2 else 'no'}, P1 {'pass' if P1 else 'FAIL'}, P4 {'pass' if P4 else 'FAIL'}")
