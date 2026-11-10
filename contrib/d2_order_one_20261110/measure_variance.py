#!/usr/bin/env python3
"""Postprocess completed bridge checkpoints; never optimize or read an active state."""
import hashlib
import json
from pathlib import Path
import time
from dmrg_screen import Z, build_model, atomic_json

HERE=Path(__file__).resolve().parent


def main():
    output=HERE/'bridge_variance.json'
    data=json.loads(output.read_text()) if output.exists() else {}
    for L in [24,48,96]:
        files=list((HERE/'dmrg').glob(f'bridge_r4_u1_k4_L{L}_*.json'))
        rr=[(p,json.loads(p.read_text())) for p in files]
        if not rr or not all(r['status']=='done' for p,r in rr):continue
        p,r=min(rr,key=lambda pr:pr[1]['E_physical'])
        sha=hashlib.sha256(p.read_bytes()).hexdigest()
        if data.get(str(L),{}).get('record_sha256')==sha:continue
        start=time.time()
        psi=Z.MPS.load(r['checkpoint']);mpo,_,shift=build_model(L,r['couplings'])
        variance,E=Z.variance(psi,mpo)
        error=abs(E+mpo.const-L*shift-r['E_physical'])
        assert error<1e-8,(p,error)
        data[str(L)]=dict(record=p.name,record_sha256=sha,variance=variance,
            variance_per_site=variance/L,energy_agreement=error,seconds=time.time()-start,
            caveat='Variance tests eigenstate quality; alone it does not identify the global ground state or a CFT.')
        atomic_json(output,data)
        print(L,data[str(L)],flush=True)


if __name__=='__main__':main()
