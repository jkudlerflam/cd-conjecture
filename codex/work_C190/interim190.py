# interim reading: p153c, L=64 (local pilot or Typhon) against L=128 (Typhon partial), by chi stage
import json, sys, os, glob
import numpy as np
sys.path.insert(0, os.getcwd())
from analyze190 import sbar, c_prof
T2 = 1.0 / 80.0   # |t2| = Delta^2/(2V), Delta = 1, V = 40

def stages(f, L):
    r = json.load(open(f)); out = {}
    for x in r['hist']:
        if 'Sprof' in x:
            S = np.array(x['Sprof'])
            out[x['chi']] = dict(S1=sbar(S, L, 1), S6=sbar(S, L, 6), S12=sbar(S, L, 12), E=x['E_phys'],
                                 n=x['n_mid'], N=x['N_tot'], res=x['res_max'], tr=x['trunc'], sw=x['sweep'],
                                 cp=c_prof(S, L))
    return r.get('status'), out

runs = {}
for f in sys.argv[1:]:
    tag = os.path.basename(f)[:-5]
    L = int(tag.split('_L')[1].split('_')[0]); P = int(tag.split('_P')[1][0])
    st, sg = stages(f, L)
    runs[tag] = (L, P, st, sg)
    print(tag, st)
    for chi, d in sorted(sg.items()):
        print(f"   chi={chi:4d} sw={d['sw']:2d} E={d['E']:.10f} S1,6,12={d['S1']:.4f},{d['S6']:.4f},{d['S12']:.4f} "
              f"c_prof={d['cp']:.3f} n_mid={d['n']:.4f} N={d['N']:.3f} res={d['res']:.1e} trunc={d['tr']:.1e}")
print()
for chi in (192, 256, 384, 512):
    for ta, (La, Pa, _, sa) in runs.items():
        for tb, (Lb, Pb, _, sb) in runs.items():
            if Lb == 2 * La and chi in sb:
                ca = chi if chi in sa else max(sa) if sa else None
                if ca is None: continue
                a, b = sa[ca], sb[chi]
                print(f"chi {ca}->{chi}  {ta} -> {tb}: c_two W=1,6,12 = " + ', '.join(
                    f"{6 * (b[k] - a[k]) / np.log(2):.3f}" for k in ('S1', 'S6', 'S12')) +
                    f"   dn_mid={b['n'] - a['n']:+.4f}")
    for ta, (La, Pa, _, sa) in runs.items():
        for tb, (Lb, Pb, _, sb) in runs.items():
            if La == Lb and Pa == 0 and Pb == 1 and chi in sa and chi in sb:
                print(f"chi {chi}  L={La}  {ta} / {tb}: E_odd-E_even = {sb[chi]['E'] - sa[chi]['E']:+.3e}, "
                      f"L*dGC/|t2| = {La * (sb[chi]['E'] - sa[chi]['E']) / T2:+.3f}")
