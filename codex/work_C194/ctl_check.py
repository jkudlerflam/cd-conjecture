"""C194: DMRG controls against the exact energies, and Delta_pb from DMRG vs exact."""
import json, glob, os
from controls_exact import E_xx, E_hr
E = {}
for f in sorted(glob.glob('out_ctl/*.jsonl')):
    recs = [json.loads(l) for l in open(f)]
    sw = [r for r in recs if r['kind'] == 'sweep']; fin = [r for r in recs if r['kind'] == 'final'][0]
    c = fin['cfg']; ex = (E_xx if '_xx_' in c['tag'] else E_hr)(c['L'], c['N'])
    E[c['tag']] = (fin['E'], ex)
    print(f"{c['tag']:18s} E={fin['E']:.10f} exact={ex:.10f} err={fin['E']-ex:+.1e} last-sweep change={sw[-1]['E']-sw[-2]['E']:+.1e} trunc={sw[-1]['trunc']:.1e}")
for nm in ('xx', 'hr'):
    for L, N in ((64, 8), (128, 4), (128, 16)):
        g = lambda k, i: E[f'ctl_{nm}_L{L}_N{N+k}'][i]
        d = [g(1, i) + g(-1, i) - 2 * g(0, i) for i in (0, 1)]
        print(f"{nm} L={L} N={N}: L*Dpb DMRG={L*d[0]:.6f} exact={L*d[1]:.6f} diff(Dpb)={d[0]-d[1]:+.1e}")
