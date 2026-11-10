#!/usr/bin/env python3
"""Global energy variance for completed, hash-stable continuation snapshots."""
from pathlib import Path
import hashlib
import json
import sys
import time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'d2_order_one_20261110'))
import dmrg_screen as D
import numpy as np

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    output=HERE/'variances.json'
    data=json.loads(output.read_text()) if output.exists() else {}
    wanted={'bridge_r4_u1_k4','b_r5_2_d1_10_q0_hm7_8','pref_ref','pref_rat','pref_zrat','attr_ref','simple_attr_7_m1_5_m7'}
    for folder in ('validation','bridge_scan/dmrg','attractive_scan/results','independent_seed','independent_sizes','simple_attractive/results'):
        for path in sorted((HERE/folder).glob('*.json')):
            raw=path.read_bytes();r=json.loads(raw)
            if r.get('status')!='done' or r.get('name') not in wanted:continue
            key=str(path.relative_to(HERE));rs=hashlib.sha256(raw).hexdigest()
            if data.get(key,{}).get('record_sha256')==rs:continue
            cp=Path(r['checkpoint']);cs=sha(cp)
            with np.load(cp,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
            if len(meta['history'])!=len(r['history']) or meta['couplings']!=r['couplings']:continue
            start=time.time();psi=D.Z.MPS.load(cp)
            mpo,_,shift=D.build_model(r['L'],r['couplings'])
            var,energy=D.Z.variance(psi,mpo)
            if sha(path)!=rs or sha(cp)!=cs:
                print('CHANGED; discarded',key,flush=True);continue
            error=abs(energy+mpo.const-r['L']*shift-r['E_physical'])
            assert error<1e-8,(key,error)
            data[key]=dict(name=r['name'],L=r['L'],parity=r['parity'],seed=r['seed'],chi=r['chi'],
                variance=var,variance_per_site=var/r['L'],energy_agreement=error,
                sweep_converged=r['convergence']['sweep_converged'],record_sha256=rs,
                checkpoint_sha256=cs,seconds=time.time()-start,
                caveat='Variance tests eigenstate quality; it alone certifies neither the ground state, parity ordering, nor a CFT. Tiny signed values are subject to floating-point cancellation.')
            D.atomic_json(output,data)
            print(key,'variance/site',var/r['L'],flush=True)

if __name__=='__main__':main()
