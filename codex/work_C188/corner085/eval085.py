# Evaluate the corner-point check against K1 to K5 and W1 to W3 of README.md.
import json, math
def rd(f): return [json.loads(l) for l in open(f) if l.strip()]
S = {}
for L in (48, 96, 192):
    sw = rd('sweeps_L%d.jsonl' % L); st = dict(zip((160, 240, 320), rd('stages_L%d.jsonl' % L)))  # one record per chi stage, in order
    a, b = sw[-2], sw[-1]
    k1 = abs(b['E'] - a['E']); k2 = abs(b['S_mid4'] - a['S_mid4'])
    k3 = abs(st[320]['S_mid4'] - st[240]['S_mid4']); k4 = b['trunc']
    k5 = st[320]['mmax']; e = st[320]['e']; m = len(e)//2
    D = sum(abs(e[i] - e[i+1]) for i in range(m-4, m+4))/8   # mid-chain bond-energy alternation
    S[L] = {c: st[c]['S_mid4'] for c in (160, 240, 320)}
    ok = k1 <= 1e-7 and k2 <= 1e-3 and k3 <= 0.01 and k4 <= 1e-8 and k5 <= 0.05
    print('L=%3d E=%.10f S4=%.7f K1=%.1e K2=%.1e K3=%.5f K4=%.1e K5=%.1e %s D_mid=%.5f c_fit=%s'
          % (L, st[320]['E'], S[L][320], k1, k2, k3, k4, k5, 'pass' if ok else 'FAIL', D,
             [round(x, 3) for x in st[320].get('c_fit', [])]))
for c in (160, 240, 320):
    print('chi=%d c_pair(48,96)=%.3f c_pair(96,192)=%.3f' % (c,
          6*(S[96][c]-S[48][c])/math.log(2), 6*(S[192][c]-S[96][c])/math.log(2)))
