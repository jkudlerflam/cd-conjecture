"""C182: summarize the Typhon output (after collection into typhon_out/).
  python3 summarize.py clusters  -> per-particle cluster ladder table
  python3 summarize.py screen    -> DMRG screen table: c_two(48->96, W=1,6), trunc, n_mid, valley weights, parity gap
"""
import numpy as np, json, glob, os, sys
from analyze import Sbar, c_prof

D = sys.argv[2] if len(sys.argv) > 2 else 'typhon_out'

def clusters():
    print("# |w| arg phi delta  e_N = E_N/N - eps_a for N = 2..Nmax (R)   onset species   B3=E2+ea-E3  B(pair,trimer)=E2+E3-E5  B(trimer,atom)=E3+ea-E4")
    for f in sorted(glob.glob(f'{D}/out/cl_*.json')):
        try: d = json.load(open(f))
        except Exception: continue
        cl = d['clusters']; ea = d['ea']
        e = {int(N): cl[N]['E']/int(N) - ea for N in cl}
        e[1] = 0.0
        best = min(e, key=lambda N: e[N])
        E = {int(N): cl[N]['E'] for N in cl}
        s = f"{d['wa']:4.2f} {d['arg']:>6} {d['phi']:.5f} {d['delta']:.4f} " + ' '.join(f"e{N}={e[N]:+.5f}" for N in sorted(e) if N > 1)
        s += f"  onset N={best}"
        if 3 in E: s += f"  B3={E[2]+ea-E[3]:+.5f}"
        if 5 in E: s += f"  B23={E[2]+E[3]-E[5]:+.5f}"
        if 4 in E:
            # trimer+atom channel: lowest 4-body level near K = K3 + k_a (the aaab candidate)
            K3 = cl['3']['K']; Kt = np.angle(np.exp(1j*(K3 + d['ka'])))
            loc = cl['4']['local']
            near = min(loc, key=lambda r: abs(np.angle(np.exp(1j*(r['K']-Kt)))))
            s += f"  E4(K~K3+ka)={near['E']-4*ea:+.5f} B(3,1)={E[3]+ea-near['E']:+.5f}  B22=2E2-E4={2*E[2]-E[4]:+.5f}"
        print(s)
    for f in sorted(glob.glob(f'{D}/out/l6_*.json')):
        try: print(open(f).read().strip())
        except Exception: pass

def screen():
    runs = {}
    for f in glob.glob(f'{D}/out/*.jsonl'):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        if recs: runs[recs[0]['tag']] = recs
    sers = sorted(set(t.rsplit('_L', 1)[0] for t in runs if '_N' not in t.rsplit('_L', 1)[1]))
    print("# series  c_two(48->96) W=1 W=6  cprof96  trunc48 trunc96  n_mid96  Wmin/Wmax  L*DeltaP(96)")
    for s in sers:
        a, b = runs.get(f'{s}_L48'), runs.get(f'{s}_L96')
        if not a or not b: continue
        ra, rb = a[-1], b[-1]
        if rb['sweep'] != len(rb['cfg']['chis'])-1: continue
        c1 = 6*(Sbar(rb['S'], 96, 1)-Sbar(ra['S'], 48, 1))/np.log(2)
        c6 = 6*(Sbar(rb['S'], 96, 6)-Sbar(ra['S'], 48, 6))/np.log(2)
        m = rb.get('meas')
        line = f"{s:28s} {c1:+.3f} {c6:+.3f} {c_prof(rb['S'], 96):+.3f} {ra['trunc']:.0e} {rb['trunc']:.0e}"
        if m:
            nm = np.mean(m['dens'][24:72]); W = sorted([m['Wa_f'], m['Wb_f']])
            line += f" {nm:.3f} {W[0]/max(W[1],1e-12):.2f}"
        p, q = runs.get(f'{s}_L96_Np1'), runs.get(f'{s}_L96_Nm1')
        if p and q:
            dp = p[-1]['E'] + q[-1]['E'] - 2*rb['E']
            line += f" {96*dp:.3f}"
        print(line)

if __name__ == "__main__":
    {'clusters': clusters, 'screen': screen}[sys.argv[1]]()
