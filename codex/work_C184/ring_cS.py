"""Ring entropy slope c_S(16,20) = 3[S(20) - S(16)]/ln(20/16) at the L = 20 points (ring20.jsonl) using ring_phi1.jsonl for L = 16.
Also the parity gap L*|E0(Pi=0)-E0(Pi=1)| and the lowest gap L*(E1-E0) at L = 12, 16, 20."""
import json, math, collections
d = collections.defaultdict(dict)
for f in ['ring_phi1.jsonl', 'ring20.jsonl']:
    for l in open(f):
        r = json.loads(l); d[(r['Delta'], r['V'])][r['L']] = r
for k in sorted(d):
    if 20 not in d[k]: continue
    r = d[k]
    cS = 3 * (r[20]['S'] - r[16]['S']) / math.log(20 / 16)
    cS1 = 3 * (r[16]['S'] - r[12]['S']) / math.log(16 / 12)
    gp = [L * abs(r[L]['E0_par'][0] - r[L]['E0_par'][1]) for L in (12, 16, 20)]
    g1 = [L * (r[L]['levels'][1] - r[L]['levels'][0]) for L in (12, 16, 20)]
    print(f"Delta={k[0]:.1f} V={k[1]:.1f}  S={r[12]['S']:.3f},{r[16]['S']:.3f},{r[20]['S']:.3f}  c_S(12,16)={cS1:.3f} c_S(16,20)={cS:.3f}  L*dpar={gp[0]:.3f},{gp[1]:.3f},{gp[2]:.3f}  L*g1={g1[0]:.3f},{g1[1]:.3f},{g1[2]:.3f}")
