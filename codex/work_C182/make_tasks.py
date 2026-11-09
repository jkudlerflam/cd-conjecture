"""C182: build the Typhon task list (stage 1: cluster ladders and the DMRG phase-diagram screen).
Stage 2 (converged runs at the best screen points) is generated on the node by select_stage2.py.
"""
import numpy as np, json, os
from pairband import valleys
from clusters import phi_for_delta

from pairband import Eb_exact

WS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
XS = [0.0, 0.5, 1.0, 2.0]               # delta / E_b
PHIS_FAR = [1.45, 1.30]
NS = [(1, 12), (1, 6), (1, 4), (1, 3), (5, 12)]
SCH48 = [16]*4 + [32]*4 + [64]*3 + [128]*3 + [192]*2
SCH96 = [16]*4 + [32]*4 + [64]*4 + [128]*3 + [192]*2 + [256]*2

def points():
    pts = []
    for wa in WS:
        Eb = Eb_exact(wa)
        for x in XS:
            ph = phi_for_delta(wa, x*Eb)
            pts.append((wa, f"x{x:g}", ph, x*Eb))
        for ph in PHIS_FAR:
            v = valleys(wa*np.exp(1j*ph))
            pts.append((wa, f"p{ph:g}", ph, (v[1][1]-v[0][1]) if len(v) > 1 else None))
    return pts

if __name__ == "__main__":
    os.makedirs('typhon/cfg', exist_ok=True)
    tasks_dmrg, tasks_cl = [], []
    meta = []
    for (wa, lab, ph, dl) in points():
        meta.append(dict(wa=wa, lab=lab, phi=ph, delta=dl, Eb=Eb_exact(wa)))
        for (a, b) in NS:
            for L, sch in ((48, SCH48), (96, SCH96)):
                N = L*a//b
                ser = f"s_w{wa:g}_{lab}_n{a}o{b}"
                cfgs = [dict(tag=f"{ser}_L{L}", L=L, N=N, wa=wa, phi=ph, chis=sch, meas_chi=10**9)]
                if L == 96:
                    for dN, s in ((1, 'p'), (-1, 'm')):
                        cfgs.append(dict(tag=f"{ser}_L{L}_N{s}1", L=L, N=N+dN, wa=wa, phi=ph, chis=sch, nomeas=True))
                for c in cfgs:
                    fn = f"typhon/cfg/{c['tag']}.json"
                    json.dump([c], open(fn, 'w'))
                    tasks_dmrg.append(f"timeout 600m python3 run.py {fn} out > logs/{c['tag']}.log 2>&1")
    # controls (C176's definitions) at n = 1/4 and 1/6
    for (a, b) in ((1, 4), (1, 6)):
        for L, sch in ((48, SCH48), (96, SCH96)):
            N = L*a//b
            for nm, extra in (('c1', dict(wa=1.0, phi=0.0, t2=0.0, V2=0.0, V3=0.0)),
                              ('c2', dict(wa=0.0, phi=0.0, V2=0.0, V3=0.0, blockade=False))):
                c = dict(tag=f"{nm}_n{a}o{b}_L{L}", L=L, N=N, chis=sch, meas_chi=10**9, **extra)
                fn = f"typhon/cfg/{c['tag']}.json"; json.dump([c], open(fn, 'w'))
                tasks_dmrg.insert(0, f"timeout 600m python3 run.py {fn} out > logs/{c['tag']}.log 2>&1")
    # cluster ladders: |w| in 1..3.5, delta/E_b grid (with a finer grid near the pair/atom crossing)
    RS = {1.0: (1500, 110, 80, 0), 1.5: (1200, 80, 48, 0), 2.0: (1000, 60, 40, 28),
          2.5: (1000, 50, 36, 26), 3.0: (1000, 46, 34, 24), 3.5: (1000, 46, 34, 24)}
    for wa, (R2, R3, R4, R5) in RS.items():
        Nmax = 5 if R5 else 4
        for x in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
            rs = [R2, R3, R4] + ([R5] if R5 else [])
            tag = f"cl_w{wa:g}_d{x:g}"
            tasks_cl.append(f"timeout 1400m python3 clusters.py {wa} d{x} {Nmax} {' '.join(map(str, rs))} > out/{tag}.json 2> logs/{tag}.err")
    # long cluster tasks first, then the dmrg screen
    open('typhon/tasks1.txt', 'w').write('\n'.join(tasks_cl + tasks_dmrg) + '\n')
    json.dump(meta, open('typhon/points.json', 'w'), indent=1)
    print(len(tasks_cl), "cluster tasks,", len(tasks_dmrg), "dmrg tasks")
    for m in meta:
        print(m)
