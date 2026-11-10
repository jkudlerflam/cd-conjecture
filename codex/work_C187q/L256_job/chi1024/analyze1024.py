"""C187d: the registered readings of README.md for the chi = 1024 continuation.
Usage, from this folder: python3 -E -s analyze1024.py > typhon_out/analyze1024.txt"""
import json, os, sys, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', '..', 'work_C181'))
import analyze as A
OUT = os.path.join(here, 'typhon_out')
Ws = (1, 6, 12); LN2 = np.log(2.0)
def recs(t):
    return [json.loads(x) for x in open(os.path.join(OUT, t + '.jsonl')).read().strip().split('\n') if x]
def stages(R):
    st = {}
    for x in R: st[x['chi']] = x          # last sweep of each chi stage
    return st
def sb(x, W): return A.Sbar(x['S'], x['cfg']['L'], W)
def extrap(st, W, chis):
    t = np.array([st[c]['trunc'] for c in chis]); s = np.array([sb(st[c], W) for c in chis])
    a, b = np.polyfit(t, s, 1); return b
def r2(st, W):
    e3 = extrap(st, W, (512, 768, 1024)); e2 = extrap(st, W, (768, 1024)); return e3, e2
names = [('ctl2 (c = 2 control)', 'bf_ctl2_L128_N24', 'ctl2_L256_N48'),
         ('phi = 1.53', 'bf_w1n316_p1.5300_L128_N24', 'w1n316_p1.5300_L256_N48'),
         ('phi = 1.555', 'bf_w1n316_p1.5550_L128_N24', 'w1n316_p1.5550_L256_N48')]
print("# stage table: tag, then for chi = 384, 512, 768, 1024 (last sweep of the stage): E, trunc, Sbar1, Sbar6, Sbar12")
alltags = [n[1] for n in names] + ['ctl2_L256_N48'] + [f'w1n316_p{p}_L256_N{N}' for p in ('1.5300', '1.5550') for N in (47, 48, 49)]
ST = {}
for t in alltags:
    R = recs(t); ST[t] = stages(R); r = R[-1]
    print(f"{t}  sweeps {len(R)}  final chi {r['chi']}  maxchi {r['maxchi']}  dE(last two sweeps) {r['E']-R[-2]['E']:+.2e}  dSbar6(last two) {sb(r,6)-sb(R[-2],6):+.2e}")
    for c in (384, 512, 768, 1024):
        if c in ST[t]:
            x = ST[t][c]
            print(f"    chi {c:4d}  E {x['E']:.10f}  trunc {x['trunc']:.2e}  Sbar {sb(x,1):.4f} {sb(x,6):.4f} {sb(x,12):.4f}")
print()
print("# convergence in chi: |Sbar6(1024) - Sbar6(768)|, bound 0.005 at L = 256 and 0.002 at L = 128")
conv = {}
for t in alltags:
    L = ST[t][1024]['cfg']['L']; d = sb(ST[t][1024], 6) - sb(ST[t][768], 6); ok = abs(d) <= (0.005 if L == 256 else 0.002)
    conv[t] = ok
    print(f"{t:30s} L = {L}  step {d:+.4f}  {'converged' if ok else 'NOT converged'}")
print()
print("# R1, matched chi = 1024 on both sizes: c_two(128 -> 256) at W = 1, 6, 12; then c_prof at each size")
R1 = {}; R2 = {}
for name, a, b in names:
    c = [6 * (sb(ST[b][1024], W) - sb(ST[a][1024], W)) / LN2 for W in Ws]; R1[name] = c
    cp = (A.c_prof(np.array(ST[a][1024]['S']), 128), A.c_prof(np.array(ST[b][1024]['S']), 256))
    print(f"{name:22s} " + ' '.join(f"{x:.3f}" for x in c) + f"   c_prof {cp[0]:.3f} -> {cp[1]:.3f}")
print("# for comparison, matched at chi = 512 and at chi = 768")
for cc in (512, 768):
    for name, a, b in names:
        c = [6 * (sb(ST[b][cc], W) - sb(ST[a][cc], W)) / LN2 for W in Ws]
        print(f"  chi {cc}  {name:22s} " + ' '.join(f"{x:.3f}" for x in c))
print()
print("# R2, linear extrapolation in discarded weight through chi = 512, 768, 1024 (three-point), with the two-point value through 768, 1024")
print("# per size: Sbar_W extrapolated (three-point, two-point), W = 1, 6, 12")
for name, a, b in names:
    c3 = []; c2 = []; unc = []
    for W in Ws:
        a3, a2 = r2(ST[a], W); b3, b2 = r2(ST[b], W)
        c3.append(6 * (b3 - a3) / LN2); c2.append(6 * (b2 - a2) / LN2); unc.append(6 * (abs(b3 - b2) + abs(a3 - a2)) / LN2)
    for t in (a, b):
        print(f"    {t:30s} " + '  '.join(f"{r2(ST[t],W)[0]:.4f} / {r2(ST[t],W)[1]:.4f}" for W in Ws))
    R2[name] = (c3, c2, unc)
    print(f"{name:22s} R2 (three-point) " + ' '.join(f"{x:.3f}" for x in c3) + "   two-point " + ' '.join(f"{x:.3f}" for x in c2)
          + "   uncertainty (sum of the two sizes) " + ' '.join(f"{x:.3f}" for x in unc))
print()
print("# calibration: ctl2 R2 within 0.25 of 2")
c3 = R2['ctl2 (c = 2 control)'][0]
print("ctl2 R2 " + ' '.join(f"{x:.3f}" for x in c3) + "   " + ("calibrated" if all(abs(x - 2) <= 0.25 for x in c3) else "NOT calibrated in every window"))
print()
print("# parity gap L * [E(N+1) + E(N-1) - 2 E(N)] at L = 256 by chi stage; L = 128 references 1.53 (phi = 1.53) and 1.50 (phi = 1.555)")
ref = {'1.5300': 1.53, '1.5550': 1.50}
for p in ('1.5300', '1.5550'):
    S = [ST[f'w1n316_p{p}_L256_N{n}'] for n in (47, 48, 49)]
    out = []
    for c in (384, 512, 768, 1024):
        if all(c in s for s in S): out.append((c, 256 * (S[0][c]['E'] + S[2][c]['E'] - 2 * S[1][c]['E'])))
    g = out[-1][1]
    print(f"phi = {p}: " + '  '.join(f"{c}: {v:+.3f}" for c, v in out) + f"   ratio to L = 128 at chi = 1024: {g/ref[p]:.3f}  ({'within' if g <= 1.5*ref[p] else 'ABOVE'} 1.5)")
print()
print("# valley polarization and mid-chain density at chi = 1024")
for name, a, b in names[1:]:
    for t in (a, b):
        r = recs(t)[-1]; m = r['meas']; L = r['cfg']['L']
        pol = abs(m['Wa_f'] - m['Wb_f']) / (m['Wa_f'] + m['Wb_f']); dn = np.array(m['dens'])
        print(f"{t:30s} p = {pol:.3f}  n_mid = {dn[L//4:3*L//4].mean():.4f}")
print()
print("# verdict rule")
cand = ['phi = 1.53', 'phi = 1.555']
r1pass = all(sum(x >= 1.6 for x in R1[n]) >= 2 for n in cand)
r2pass = all(sum(x >= 1.75 for x in R2[n][0]) >= 2 for n in cand)
fail = any(sum(x <= 1.3 for x in R1[n]) >= 2 for n in cand)
print(f"R1 >= 1.6 in two of three windows at both phi: {r1pass}")
print(f"R2 >= 1.75 in two of three windows at both phi: {r2pass}")
print(f"R1 <= 1.3 in two windows at either phi: {fail}")
print("verdict: " + ("E2" if fail else ("E1" if (r1pass and r2pass) else "E3")))
print()
print("# not registered: the same reading with the N = 47 and N = 49 runs at L = 256 against the N = 24 partner (chi = 1024), W = 1, 6, 12")
for p in ('1.5300', '1.5550'):
    a = ST[f'bf_w1n316_p{p}_L128_N24'][1024]
    for N in (47, 48, 49):
        b = ST[f'w1n316_p{p}_L256_N{N}'][1024]
        print(f"phi = {p}  N = {N}  Sbar6 {sb(b,6):.4f}  c_two " + ' '.join(f"{6*(sb(b,W)-sb(a,W))/LN2:.3f}" for W in Ws))
print()
print("# not registered: shell effect in N. d2 = Sbar_W(N) - [Sbar_W(N-1) + Sbar_W(N+1)]/2 at fixed L, W = 1, 6, 12.")
print("# L = 256 at chi = 1024 (this job), L = 128 at chi = 384 (C181 backfill, the only chi where N = 23 and 25 exist)")
BF = os.path.join(here, '..', '..', '..', 'work_C181', 'typhon_out_bf', 'out_bf')
def bf_stage(t):
    st = {}
    for x in open(os.path.join(BF, t + '.jsonl')).read().strip().split('\n'):
        if x:
            r = json.loads(x); st[r['chi']] = r
    return st
D2 = {}
for p in ('1.5300', '1.5550'):
    b = [ST[f'w1n316_p{p}_L256_N{N}'][1024] for N in (47, 48, 49)]
    a = [bf_stage(f'bf_w1n316_p{p}_L128_N{N}')[384] for N in (23, 24, 25)]
    d256 = [sb(b[1], W) - 0.5 * (sb(b[0], W) + sb(b[2], W)) for W in Ws]
    d128 = [sb(a[1], W) - 0.5 * (sb(a[0], W) + sb(a[2], W)) for W in Ws]
    D2[p] = (d128, d256)
    print(f"phi = {p}  L = 128: Sbar6(N = 23, 24, 25) " + ' '.join(f"{sb(x,6):.4f}" for x in a) + "   d2 " + ' '.join(f"{x:+.4f}" for x in d128))
    print(f"phi = {p}  L = 256: Sbar6(N = 47, 48, 49) " + ' '.join(f"{sb(x,6):.4f}" for x in b) + "   d2 " + ' '.join(f"{x:+.4f}" for x in d256))
    print(f"phi = {p}  shift of c_two from removing d2 at both sizes, 6 (d2_128 - d2_256) / ln 2: " + ' '.join(f"{6*(x-y)/LN2:+.3f}" for x, y in zip(d128, d256))
          + "   R1 with the shift: " + ' '.join(f"{c + 6*(x-y)/LN2:.3f}" for c, x, y in zip(R1['phi = 1.53' if p == '1.5300' else 'phi = 1.555'], d128, d256)))
