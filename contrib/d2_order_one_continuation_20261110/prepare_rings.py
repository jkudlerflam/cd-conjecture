#!/usr/bin/env python3
"""Freeze the registered attractive-ring input from the branch's exact points."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
p=json.loads((HERE/'attractive_scan/points.json').read_text())
cfg=dict(sizes=[14,16,20],sectors=[0,1],periodic=True,k=4,tol=1e-11,
         points=[dict(name=name,couplings=p[name]) for name in ['pref_ref','pref_rat']])
(HERE/'preferred_ring_grid.json').write_text(json.dumps(cfg,indent=2)+'\n')
