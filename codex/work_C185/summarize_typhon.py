"""C185 Typhon readout (coordinator, 2026-11-10): R6' two-size block readings, the dmv1 parity gap of R6, and R7's
full-chain S_mid against the lambda = 0 reference. Usage: python3 summarize_typhon.py (reads typhon_out, typhon_out_bf)."""
import json, math, glob, os
def last(f):
    return json.loads(open(f).read().strip().split('\n')[-1])
def sb(S):
    m = len(S) // 2; return sum(S[m - 1:m + 2]) / 3
k = 6 / math.log(2)
print("R6' block model (Smid / Sbar of three middle cuts, final truncation)")
bk = {}
for f in sorted(glob.glob('typhon_out/outbk/*.jsonl') + glob.glob('typhon_out_bf/outbk/*.jsonl')):
    r = last(f); t = os.path.basename(f)[:-6]; bk[t] = r
    print(f"  {t:24s} n_mid {r['n_mid']:.3f}  Smid {r['Smid']:.4f}  Sbar {sb(r['S']):.4f}  trunc {r['hist'][-1][2]:.1e}")
for base in ('bk_k0v0_p1571', 'bk_k1v1_p1571', 'bk_k0v1_p1571', 'bk_k1v1_p1625', 'bk_k1v1_p1640', 'bk_k1v1_p1655'):
    a, b = bk.get(base + '_Lb16'), bk.get(base + '_Lb32')
    print(f"  {base}: c2(48 -> 96) = " + (f"{k * (sb(b['S']) - sb(a['S'])):.3f} (Smid: {k * (b['Smid'] - a['Smid']):.3f})" if a and b else 'no L = 96 run'))
print('R6 parity gap at dmv1 (L Delta_P, Delta_P = E(N+1) + E(N-1) - 2E(N))')
for L in (64, 128):
    E = {s: last(f'typhon_out/out/dmv1_L{L}{s}.jsonl')['E'] for s in ('', '_Np', '_Nm')}
    d = E['_Np'] + E['_Nm'] - 2 * E['']
    print(f'  L = {L}: Delta_P = {d:.6f}, L Delta_P = {L * d:.4f}')
print('R7 full chain, V/Delta = 8, S_mid at L = 64 against the lambda = 0 reference')
oc = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob('typhon_out/oc/*.json') + glob.glob('typhon_out_bf/oc/*.json')}
for t, d in sorted(oc.items()):
    print(f"  {t:16s} " + ', '.join(f'{kk}={d[kk]:.4f}' if isinstance(d.get(kk), float) else f'{kk}={d.get(kk)}' for kk in ('Smid', 'c_prof', 'trunc', 'L') if kk in d))
