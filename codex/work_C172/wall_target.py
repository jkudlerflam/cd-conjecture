#!/usr/bin/env python3
"""C172: targeted Ising-limit test. Phase-aligned pitch (a near pi, where the charge-2 threshold |P| = |A| is sharp),
|P|/|A| = 0.8 .. 1.0, random charge-1 vectors of norm up to 2 and random hz. Can charge-1 terms lower the threshold?"""
import numpy as np
from ed_kink import fam
from wall_pt import configs_one, configs_two, projected
rng = np.random.default_rng(7)
L2, L1 = 28, 29
c1 = configs_one(L1); c2 = configs_two(L2)
for ratio in (0.8, 0.9, 0.95, 1.0, 1.05, 1.2):
    bs = []
    for it in range(12):
        a = np.pi + rng.uniform(-0.3, 0.3)
        amp1 = [0, 0.5, 1.0, 2.0][it % 4]
        v = rng.standard_normal(6); v *= amp1 / np.linalg.norm(v)
        G, Gp, Dx, Dy, hx, hy = v
        hz = 0.0 if it < 4 else rng.uniform(-0.5, 0.5)
        J, B = fam(a=a, g=ratio, G=G, Gp=Gp, Dx=Dx, Dy=Dy, hx=hx, hy=hy, hz=hz)
        e1 = np.linalg.eigvalsh(projected(L1, J, B, c1))[0]
        e2 = np.linalg.eigvalsh(projected(L2, J, B, c2))[0]
        bs.append((round(e2 - 2 * e1, 4), round(a - np.pi, 2), amp1, round(hz, 2)))
    print(f'|P|/|A| = {ratio}: min b = {min(b[0] for b in bs):+.4f}; (b, a-pi, |charge1|, hz):', bs, flush=True)
