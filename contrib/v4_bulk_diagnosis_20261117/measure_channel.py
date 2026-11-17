#!/usr/bin/env python3
"""Secondary XX+YY valley-order channel, explicitly including q=0."""
import argparse,time,json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from phase_common import HERE,ROOT,Z,write_json
from observables import environments,operator_matrix,spectral_metrics


def run(checkpoint,tag,outroot=HERE):
    start=time.time();psi=Z.MPS.load(str(checkpoint));L=psi.L;gp=Z.gauge_phases(L,'stag')
    P=np.kron(Z.PX,Z.PX)+np.kron(Z.PY,Z.PY)
    ops=[Z.rotate_bond(P,gp[j],gp[j+1]) for j in range(L-1)]
    lo,hi=L//4,3*L//4;ids=list(range(lo,hi-1))
    d=operator_matrix(psi,environments(psi),ops,2,ids)
    d['half']=spectral_metrics(d['profile'],d['connected'],ids)
    qlo,qhi=3*L//8,5*L//8;ix=[i for i,j in enumerate(ids) if qlo<=j and j+2<=qhi]
    d['quarter']=spectral_metrics(d['profile'],np.array(d['connected'])[np.ix_(ix,ix)],[ids[i] for i in ix])
    for window in ('half','quarter'):
        met=d[window];k=int(np.argmax(met['S_connected']))
        met['connected_peak_including_zero']=dict(q=met['q_grid'][k],S=met['S_connected'][k],S_over_M=met['S_connected'][k]/met['M'])
        met['S_zero_over_M']=met['S_connected'][0]/met['M']
        met['uniform_order_weight']=met['S_zero_over_M']+met['mean']**2
    d.update(tag=tag,L=L,P=psi.P,definition='Physical XX+YY; secondary, unregistered valley-polarization/bond-order diagnostic.',
        checkpoint=str(checkpoint),runtime_seconds=time.time()-start)
    write_json(outroot/'channels'/(tag+'.json'),d)
    print(tag,'S(0)/M',d['half']['S_zero_over_M'],'seconds',d['runtime_seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('checkpoint',type=Path);ap.add_argument('--tag',required=True)
    ap.add_argument('--run-root',type=Path,default=HERE);args=ap.parse_args()
    with threadpool_limits(limits=1):run(args.checkpoint,args.tag,args.run_root)
