"""C175 Typhon job: two-size readings against the registration in RESULT.md.
Usage: python3 summarize_typhon.py [outdir]   (default typhon_out)
Estimator (registered): c = 6 [Sbar(192) - Sbar(96)] / ln 2, Sbar the mean of the three middle block cuts.
Also printed: the same with 1, 5 and 7 middle cuts and with the mean over the middle half of the chain,
the profile fit c_cc stored by block_dmrg.py, the final truncation, and the left-right asymmetry of the profile."""
import json, math, sys
d = sys.argv[1] if len(sys.argv) > 1 else 'typhon_out'
TAGS = ['W1_225', 'W1_224', 'W15_255', 'W2_290', 'W1M_215', 'W2M_250', 'W1R_220']
def last(f):
    t = open(f).read().strip()
    return json.loads(t.split('\n')[-1]) if t else None
def sb(S, k):
    m = len(S) // 2
    return sum(S[m - k:m + k + 1]) / (2 * k + 1)
def half(S):
    n = len(S); q = S[n // 4:3 * n // 4 + 1]
    return sum(q) / len(q)
print('tag       L   n_mid   Smid    Sbar3   c_cc   trunc    dE(last 2 sweeps)  asym     hours')
R = {}
for tag in TAGS:
    for Lb in (32, 64):
        r = last(f'{d}/bd_{tag}_{Lb}.jsonl')
        if r is None:
            print(f'{tag:8s} {3*Lb:3d}  empty file, the run hit its 1400 minute limit'); continue
        S = r['S']; n = len(S); h = r['hist']
        asym = max(abs(S[i] - S[n - 1 - i]) for i in range(n // 2))
        R[(tag, Lb)] = r
        print(f"{tag:8s} {3*Lb:3d}  {r['n_mid']:.4f}  {r['Smid']:.4f}  {sb(S,1):.4f}  {r['c_cc']:5.2f}  {h[-1][2]:.1e}  {h[-1][1]-h[-2][1]:+.1e}           {asym:.1e}  {r['secs']/3600:.1f}")
print()
print('tag       c(96->192): 1 cut, 3 cuts (registered), 5 cuts, 7 cuts, middle half;  reportable (trunc < 1e-6 at both sizes)')
for tag in TAGS:
    if (tag, 32) not in R or (tag, 64) not in R:
        print(f'{tag:8s} no pair'); continue
    a, b = R[(tag, 32)], R[(tag, 64)]
    cs = [6 * (sb(b['S'], k) - sb(a['S'], k)) / math.log(2) for k in (0, 1, 2, 3)] + [6 * (half(b['S']) - half(a['S'])) / math.log(2)]
    ok = max(a['hist'][-1][2], b['hist'][-1][2]) < 1e-6
    print(f"{tag:8s} " + '  '.join(f'{c:6.3f}' for c in cs) + f"   {'yes' if ok else 'no'}")
