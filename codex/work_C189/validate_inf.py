"""C189 validation of idmrg.py against C176's fibdmrg.py on finite chains.
The first iDMRG iteration optimizes an open 2q-site window with trivial environments, i.e. the finite chain of
2q sites at N = 2p, with idmrg.py's own two-site code (not fibdmrg's). Its energy must match fibdmrg's converged
energy for the same chain. Usage: python3 validate_inf.py"""
import numpy as np, sys
from fibdmrg import Chain, DMRG
from idmrg import Model, IDMRG

cases = [
    ('sym n=1/8, 16 sites', dict(p=1, q=8, w=1j, gauge=1)),
    ('sym n=3/16, 32 sites', dict(p=3, q=16, w=1j, gauge=1)),
    ('c1 hard rods n=1/4, 8 sites', dict(p=1, q=4, w=1.0, t2=0.0, V2=0.0, V3=0.0)),
    ('sym phi=1.50 complex, 16 sites', dict(p=1, q=8, w=np.exp(1.5j), gauge=0)),
]
for name, c in cases:
    p, q = c.pop('p'), c.pop('q')
    w = c.pop('w')
    M = Model(p, q, w, **c)
    dm = IDMRG(M, chi0=8)          # fibdmrg only seeds the window at chi 8 (two sweeps at 8, 16)
    dm.sweep_window(256, nsw=6)    # idmrg's own sweeps at chi 256
    r = dm.iterate(256)
    ch = Chain(2 * q, 2 * p, w, **c)
    ref = DMRG(ch, seed=3)
    for chi in (16, 32, 64, 128, 256, 256, 256):
        Eref, tr = ref.sweep(chi)
    S_ref = ref.S[q - 1]
    print(f"{name}: idmrg E = {r['E']:.12f}, fibdmrg E = {Eref:.12f}, diff = {r['E'] - Eref:.1e}; "
          f"S(mid) idmrg = {r['S']:.10f}, fibdmrg = {S_ref:.10f}")
