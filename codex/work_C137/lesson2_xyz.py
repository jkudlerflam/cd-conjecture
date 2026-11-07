"""Lesson 2 on an interacting near-critical massive family at d = 2:
XYZ bond (1+gam) SxSx + (1-gam) SySy + Delta SzSz, Delta = 0.5, gam = c/n.
u = 0 is optimal by the Z2 x Z2 pi-rotation symmetry (any optimal u can be averaged to zero),
and the bond is reflection symmetric. Usage: python lesson2_xyz.py c nmin nmax"""
import sys, warnings, json
warnings.filterwarnings('ignore')
from common import *
c = float(sys.argv[1]); nmin, nmax = int(sys.argv[2]), int(sys.argv[3]); Delta = 0.5
path = HERE / 'results' / f'xyz_c{c}.json'
recs = json.load(open(path)) if path.exists() else []
done = {r['n'] for r in recs}
for n in range(nmin, nmax + 1):
    if n in done: continue
    gam = c / n
    h = xyz(1 + gam, 1 - gam, Delta)
    r = optimize_window(WindowFamily(h, n, udirs=(), reflect=True), maxit=300)
    rec = dict(c=c, n=n, gam=gam, Delta=Delta, opt=r['value'], copt=r['c'].tolist(), pred_gap=r['pred_gap'])
    for wn, cw in [('ssd', ssd_window(n)), ('sine', sine_window(n))]:
        rec[wn] = float(lowest(WindowFamily(h, n, udirs=(), cfix=cw).matrix(np.zeros(0)), 2)[0][0])
    recs.append(rec); recs.sort(key=lambda r: r['n']); write_json(path, recs)
    print(c, n, rec['opt'], rec['ssd'], rec['sine'], flush=True)
