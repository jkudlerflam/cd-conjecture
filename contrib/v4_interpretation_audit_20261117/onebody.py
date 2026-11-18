#!/usr/bin/env python3
"""Bounded read-only measurements; never calls DMRG or modifies a state."""
import json,time,sys,signal,resource
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,ROOT,PRIOR,OLD,read,write,digest
sys.path.insert(0,str(ROOT/'codex/work_C190'))
import z2dmrg as Z


def timeout(*_):raise TimeoutError('60-second per-state measurement cap')


def measure(cp,tag,L,P,chi):
    started=time.monotonic();signal.alarm(60)
    try:
        psi=Z.MPS.load(str(cp));assert psi.L==L and psi.P==P and psi.check()
        m=Z.measure(psi,None,gph=Z.gauge_phases(L,'stag'),gwin=L//4,cp_rmax=2)
    finally:signal.alarm(0)
    for label in ('b','f'):
        G=np.array(m['G'+label+'_re'])+1j*np.array(m['G'+label+'_im'])
        assert np.max(abs(G-G.conj().T))<1e-10
        assert np.linalg.eigvalsh(G).min()>-1e-8
        if label=='f':assert np.linalg.eigvalsh(G).max()<1+1e-8
    row=dict(tag=tag,L=L,P=P,chi=chi,checkpoint=str(cp.relative_to(ROOT)),
        checkpoint_sha256=digest(cp),elapsed_seconds=time.monotonic()-started,measurement=m,
        peak_RSS_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9,
        scope='Read-only C190 contractions; original unrotated physical frame; central-quarter window. No optimization.')
    write(HERE/'onebody'/(tag+'.json'),row)
    print(tag,round(row['elapsed_seconds'],3),'seconds',round(row['peak_RSS_GB'],3),'GB',flush=True)
    return row


def main():
    signal.signal(signal.SIGALRM,timeout);start=time.monotonic()
    with threadpool_limits(limits=1):
        tag='reused_L96_P0_chi256';r=measure(OLD/'runtime/legacy_V4_a1_L96_P0.npz',tag,96,0,256)
        old=read(OLD/'outputs/legacy_V4_a1_L96_P0.json')['meas']
        err={k:float(np.max(abs(np.array(r['measurement'][k])-np.array(old[k]))))
            for k in ('dens','Gb_re','Gb_im','Gf_re','Gf_im')}
        assert max(err.values())<1e-8
        projection=r['elapsed_seconds']*(320/256)**3*(36/24)**2*(144/96)
        bench=dict(errors_against_archive=err,benchmark_seconds=r['elapsed_seconds'],
            projected_largest_seconds=projection,peak_RSS_GB=r['peak_RSS_GB'],
            permitted=projection<60 and r['peak_RSS_GB']<2,total_cap_seconds=180)
        write(HERE/'measurement_benchmark.json',bench)
        # Reuse the other archived central-quarter matrix without remeasurement.
        d=read(OLD/'outputs/legacy_V4_a1_L96_P1.json')
        write(HERE/'onebody/reused_L96_P1_chi256.json',dict(tag='reused_L96_P1_chi256',L=96,P=1,chi=256,
            measurement=d['meas'],elapsed_seconds=0,scope='Reused previously validated archived central-quarter matrix.'))
        if not bench['permitted']:return
        cases=[(48,P,192,f'v4_L48_P{P}_s1_chi192') for P in (0,1)]
        cases += [(144,P,256,f'v4_L144_P{P}_s1_chi256') for P in (0,1)]
        cases += [(144,1,320,'v4_L144_P1_s1_c320')]
        for L,P,chi,tag in cases:
            if time.monotonic()-start+projection>180:raise TimeoutError('Total measurement budget exhausted')
            r=measure(PRIOR/'runtime'/(tag+'.npz'),tag,L,P,chi)
            if r['peak_RSS_GB']>=2:raise MemoryError('Measurement RSS cap reached; stop further work')
    write(HERE/'measurement_runtime.json',dict(total_seconds=time.monotonic()-start,
        workers=1,BLAS_threads=1,optimization_sweeps=0))


if __name__=='__main__':main()
