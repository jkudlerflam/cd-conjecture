"""C187b: the registered L = 128 -> 256 reading (README.md). Pairs C181's backfill L = 128 runs
(../../work_C181/typhon_out_bf/out_bf/bf_*) with typhon_out/*.jsonl.
Usage: python3 -E -s pair256.py   (run from this folder)"""
import json, os, sys, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'work_C181'))
import analyze as A
BF = os.path.join(here, '..', '..', 'work_C181', 'typhon_out_bf', 'out_bf')
OUT = os.path.join(here, 'typhon_out')
def recs(f):
    return [json.loads(x) for x in open(f).read().strip().split('\n') if x]
def last(f):
    return recs(f)[-1]
Ws = (1, 6, 12)
def ctwo(r1, r2):
    L1, L2 = r1['cfg']['L'], r2['cfg']['L']
    return [6 * (A.Sbar(r2['S'], L2, W) - A.Sbar(r1['S'], L1, W)) / np.log(L2 / L1) for W in Ws]
print("# convergence at L = 256: tag, sweeps, chi, E, dE(last two sweeps), trunc, maxchi, Sbar6, dSbar6(last two), Sbar6 at the end of each chi stage")
tags = ['ctl1fix_L256_N48', 'ctl1_L256_N48', 'ctl2_L256_N48'] + [f'w1n316_p{p}_L256_N{N}' for p in ('1.5300', '1.5550') for N in (47, 48, 49)]
for t in tags:
    R = recs(os.path.join(OUT, t + '.jsonl')); r = R[-1]; L = r['cfg']['L']
    st = {}
    for x in R: st[x['chi']] = A.Sbar(x['S'], L, 6)
    print(f"{t:26s} {r['sweep']:2d} {r['chi']:4d} {r['E']:.9f} {r['E']-R[-2]['E']:+.1e} {r['trunc']:.1e} {r['maxchi']:4d} "
          f"{A.Sbar(r['S'], L, 6):.4f} {A.Sbar(r['S'], L, 6)-A.Sbar(R[-2]['S'], L, 6):+.1e}  " + ' '.join(f"{c}:{s:.4f}" for c, s in sorted(st.items())))
print("# c_two(128 -> 256), W = 1, 6, 12, and c_prof at each size")
pairs = [('ctl1fix (c = 1 calibration)', 'bf_ctl1_L128_N24', 'ctl1fix_L256_N48'),
         ('ctl1 (config error, t2 = 0)', 'bf_ctl1_L128_N24', 'ctl1_L256_N48'),
         ('ctl2 (c = 2 calibration)', 'bf_ctl2_L128_N24', 'ctl2_L256_N48'),
         ('phi = 1.53', 'bf_w1n316_p1.5300_L128_N24', 'w1n316_p1.5300_L256_N48'),
         ('phi = 1.555', 'bf_w1n316_p1.5550_L128_N24', 'w1n316_p1.5550_L256_N48')]
for name, a, b in pairs:
    r1, r2 = last(os.path.join(BF, a + '.jsonl')), last(os.path.join(OUT, b + '.jsonl'))
    c = ctwo(r1, r2)
    print(f"{name:30s} " + ' '.join(f"{x:.3f}" for x in c) + f"   c_prof {A.c_prof(np.array(r1['S']), 128):.3f} -> {A.c_prof(np.array(r2['S']), 256):.3f}"
          f"   Sbar6 {A.Sbar(r1['S'],128,6):.4f} -> {A.Sbar(r2['S'],256,6):.4f}")
print("# parity gap E(N+1) + E(N-1) - 2 E(N), times L")
for p in ('1.5300', '1.5550'):
    out = []
    for L, N, d, pre in ((128, 24, BF, 'bf_'), (256, 48, OUT, '')):
        E = [last(os.path.join(d, f'{pre}w1n316_p{p}_L{L}_N{n}.jsonl'))['E'] for n in (N - 1, N, N + 1)]
        g = E[0] + E[2] - 2 * E[1]; out.append((L, g, L * g))
    print(f"phi = {p}: " + '  '.join(f"L={L} gap {g:.6f} L*gap {lg:.3f}" for L, g, lg in out) + f"   ratio {out[1][2]/out[0][2]:.3f}")
print("# valley polarization and mid-chain density at L = 256")
for p in ('1.5300', '1.5550'):
    for L, N, d, pre in ((128, 24, BF, 'bf_'), (256, 48, OUT, '')):
        m = last(os.path.join(d, f'{pre}w1n316_p{p}_L{L}_N{N}.jsonl'))['meas']
        pol = abs(m['Wa_f'] - m['Wb_f']) / (m['Wa_f'] + m['Wb_f']); dn = np.array(m['dens'])
        print(f"phi = {p} L={L}: p = {pol:.3f}, n_mid = {dn[L//4:3*L//4].mean():.4f}")
print("# chi dependence: L * parity gap from the last sweep of each chi stage (all three N at the same chi), and Sbar6 of the main run")
def stage(f):
    st = {}
    for x in recs(f): st[x['chi']] = x
    return st
for p in ('1.5300', '1.5550'):
    for L, N, d, pre in ((128, 24, BF, 'bf_'), (256, 48, OUT, '')):
        S = [stage(os.path.join(d, f'{pre}w1n316_p{p}_L{L}_N{n}.jsonl')) for n in (N - 1, N, N + 1)]
        chis = sorted(set(S[0]) & set(S[1]) & set(S[2]))
        print(f"phi = {p} L={L}: " + '  '.join(f"{c}: {L*(S[0][c]['E']+S[2][c]['E']-2*S[1][c]['E']):+.3f} / {A.Sbar(S[1][c]['S'], L, 6):.4f}" for c in chis if c >= 128))
print("# c_two(128 -> 256) at matched chi = 384 on both sizes, W = 1, 6, 12")
for name, a, b in pairs:
    s1, s2 = stage(os.path.join(BF, a + '.jsonl')), stage(os.path.join(OUT, b + '.jsonl'))
    print(f"{name:30s} " + ' '.join(f"{x:.3f}" for x in ctwo(s1[384], s2[384])))
print("# for the trend: c_two(64 -> 128) from the backfill files alone, W = 1, 6, 12 (not part of the registered rule)")
for name, a in (('ctl1', 'bf_ctl1'), ('ctl2', 'bf_ctl2'), ('phi = 1.53', 'bf_w1n316_p1.5300'), ('phi = 1.555', 'bf_w1n316_p1.5550')):
    r1, r2 = last(os.path.join(BF, a + '_L64_N12.jsonl')), last(os.path.join(BF, a + '_L128_N24.jsonl'))
    print(f"{name:30s} " + ' '.join(f"{x:.3f}" for x in ctwo(r1, r2)) + f"   chi {r1['chi']}, {r2['chi']}  trunc {r1['trunc']:.1e}, {r2['trunc']:.1e}")
