"""Fixed V=4 physical model and unchanged C190 parity solver."""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'contrib/simple_qubit_candidates_20261117'
sys.path.insert(0,str(ROOT/'codex/work_C190'))
import z2dmrg as Z
import analyze190 as Ana


def build(cfg):
    L=cfg['L'];hz=-7/4;shift=3/4
    h=Z.pauli_bond(Jx=.5,Jy=-.5,Jz=1.,D=1/16)
    gp=Z.gauge_phases(L,'stag')
    hs=[Z.rotate_bond(h,gp[j],gp[j+1]) for j in range(L-1)]
    h1=[hz*Z.ZZ_+shift*Z.I2 for _ in range(L)]
    if cfg.get('boundary','full')=='half_end_fields':
        for j in (0,L-1):h1[j]=(hz/2)*Z.ZZ_+shift*Z.I2
    else:assert cfg.get('boundary','full')=='full'
    return Z.MPO(L,hs,h1),gp,shift,dict(Jx=.5,Jy=-.5,Jz=1.,D=1/16,hz=hz,V=4,kappa=1,a=1)


def write_json(path,obj):
    import json,os
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
