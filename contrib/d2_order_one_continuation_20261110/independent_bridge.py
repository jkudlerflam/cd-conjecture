#!/usr/bin/env python3
"""Registered independent seed for the R=8.75 point."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'d2_order_one_20261110'))
import dmrg_screen as D
point=dict(name='b_r5_2_d1_10_q0_hm7_8',couplings=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875))
output=HERE/'independent_seed'
output.mkdir(exist_ok=True)
checkpoint=Path('/private/tmp/d2-order-one-continuation-checkpoints/independent_seed')
checkpoint.mkdir(parents=True,exist_ok=True)
if __name__=='__main__':
    path=output/'b_r5_2_d1_10_q0_hm7_8_L48_P0_s23.json'
    if path.exists():raise RuntimeError('Preserve completed or interrupted result; explicit resume registration required.')
    D.run_point(point,48,0,chi=192,seed=23,stages=[16,32,64,128,192],
                max_sweeps=8,output_dir=output,checkpoint_dir=checkpoint)
