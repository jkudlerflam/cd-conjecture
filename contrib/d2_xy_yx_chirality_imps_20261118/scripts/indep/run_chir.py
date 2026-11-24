#!/usr/bin/env python3
"""Independent-code check of the four-term chain (2026-11-24).  Runs codex/work_C190's Z2 DMRG
(z2dmrg.py / z2dmrg_xp.py, run190g.py driver, not TeNPy) on the pi/4-rotated chain of REPORT 9.213:
    H = sum [ (XX - YY)/2 + (delta/2)(XY - YX) + ZZ ] - h sum Z,
i.e. the C190 bond family with Jx = 1/2, Jy = -1/2, Jz = 1, D = delta/2, field hz = -h on every site.
Model name 'chir<delta>' (e.g. chir0.1875); config 'hz' = -h.  Everything else is run190g.py unchanged.
usage: python3 run_chir.py <config.json> <outdir>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models as Mo
_orig = Mo.bond_family


def bond_family(name):
    if name.startswith('chir'):
        d = float(name[4:])
        return dict(Jx=0.5, Jy=-0.5, Jz=1.0, D=d / 2, hz0=None)
    return _orig(name)


Mo.bond_family = bond_family
import run190g
run190g.main()
