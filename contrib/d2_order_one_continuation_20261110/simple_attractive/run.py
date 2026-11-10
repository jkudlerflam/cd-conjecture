#!/usr/bin/env python3
"""One sequential worker for the preregistered small-integer point."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'d2_order_one_20261110'))
import dmrg_screen as D
import ed_screen as E
import json
point=dict(name='simple_attr_7_m1_5_m7',couplings=dict(Jx=7/12,Jy=-1/12,Dz=5/12,Jz=-7/12,hz=0.))
out=HERE/'results';out.mkdir(exist_ok=True)
checkpoint=Path('/private/tmp/d2-order-one-continuation-checkpoints/simple_attractive')
checkpoint.mkdir(parents=True,exist_ok=True)

if __name__=='__main__':
    D.atomic_json(HERE/'point.json',dict(**point,metrics=E.coupling_metrics(point['couplings'])))
    for L in (48,96):
        for P in (0,1):
            path=out/f'{point["name"]}_L{L}_P{P}_s17.json'
            if path.exists():
                r=json.loads(path.read_text())
                if r.get('status')=='done':continue
                raise RuntimeError('Interrupted record: explicit resume registration required.')
            print('START',L,P,flush=True)
            D.run_point(point,L,P,chi=128,seed=17,stages=[16,32,64,128],
                max_sweeps=8,output_dir=out,checkpoint_dir=checkpoint)
            print('END',L,P,flush=True)
