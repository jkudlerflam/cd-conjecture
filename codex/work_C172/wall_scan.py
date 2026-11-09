#!/usr/bin/env python3
"""C172: random scan of the first-order wall problem (wall_pt.py) with all harmonic types.
Tests the derived statement: if |P| < |A| (in-plane det > 0) no two-wall bound state below 2 E1min."""
import sys, numpy as np
from ed_kink import fam
from wall_pt import configs_one, configs_two, projected
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
L2, L1 = 28, 29
c1 = configs_one(L1); c2 = configs_two(L2)
for it in range(N):
    a = rng.uniform(0, 2 * np.pi)
    ratio = rng.uniform(0.0, 1.4)                      # |P|/|A|
    phi = rng.uniform(0, 2 * np.pi)
    g, gxy = ratio * np.cos(phi), ratio * np.sin(phi)  # |P| = 2 sqrt(g^2+gxy^2), |A| = 2
    amp1 = rng.uniform(0, 1.5)
    v = rng.standard_normal(6); v *= amp1 / np.linalg.norm(v)
    G, Gp, Dx, Dy, hx, hy = v
    hz = rng.uniform(-1.0, 1.0)
    J, B = fam(a=a, g=g, gxy=gxy, G=G, Gp=Gp, Dx=Dx, Dy=Dy, hx=hx, hy=hy, hz=hz)
    e1 = np.linalg.eigvalsh(projected(L1, J, B, c1))[0]
    e2 = np.linalg.eigvalsh(projected(L2, J, B, c2))[0]
    print(f'{it:3d} a={a:.3f} P/A={ratio:.3f} amp1={amp1:.3f} hz={hz:+.3f}  E1={e1:.5f} E2={e2:.5f} b={e2-2*e1:+.5f}', flush=True)
