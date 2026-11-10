#!/usr/bin/env python3
"""Registered intermediate chain length, without parameter retuning."""
from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'d2_order_one_20261110'))
import dmrg_screen as D
point=dict(name='b_r5_2_d1_10_q0_hm7_8',couplings=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875))
output=HERE/'independent_sizes';output.mkdir(exist_ok=True)
checkpoint=Path('/private/tmp/d2-order-one-continuation-checkpoints/independent_sizes')
checkpoint.mkdir(parents=True,exist_ok=True)
if __name__=='__main__':
    for P in (0,1):
        path=output/f'{point["name"]}_L64_P{P}_s17.json'
        if path.exists():
            r=json.loads(path.read_text())
            if r.get('status')=='done':continue
            raise RuntimeError('Interrupted record: explicit resume registration required.')
        D.run_point(point,64,P,chi=192,seed=17,stages=[16,32,64,128,192],
                    max_sweeps=8,output_dir=output,checkpoint_dir=checkpoint)
