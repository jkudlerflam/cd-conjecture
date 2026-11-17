#!/usr/bin/env python3
"""Audit and measure a preserved state without altering its checkpoint."""
import argparse,hashlib,json,time,resource,sys
from pathlib import Path
from threadpoolctl import threadpool_limits
from phase_common import HERE,OLD,Z,build,write_json
from observables import observe


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parity',type=int,required=True);args=ap.parse_args()
    P=args.parity;tag=f'legacy_V4_a1_L96_P{P}';cp=OLD/'runtime'/(tag+'.npz')
    archived=json.loads((OLD/'outputs'/(tag+'.json')).read_text());cfg=archived['cfg']
    assert (cfg['V'],cfg['a'],cfg['kappa'],cfg['L'],cfg['P'])==(4,1,1,96,P)
    with threadpool_limits(limits=1):
        psi=Z.MPS.load(str(cp));assert psi.L==96 and psi.P==P and psi.check()
        mpo,gp,shift,b=build(cfg);energy=Z.energy(psi,mpo)+mpo.const-96*shift
        assert abs(energy-archived['hist'][-1]['E_phys'])<1e-8
        obs=observe(psi,gp,mpo.h2)
        assert max(abs(x-y) for x,y in zip(obs['operators']['density']['profile'],archived['meas']['dens']))<1e-8
        obs.update(tag=f'reused_L96_P{P}_chi256',chi=256,boundary='full',source=str(cp.relative_to(HERE.parents[1])),
            source_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),E=energy,
            archived_output=str((OLD/'outputs'/(tag+'.json')).relative_to(HERE.parents[1])),
            peak_RSS_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9)
        write_json(HERE/'observations'/(obs['tag']+'.json'),obs)
        print(obs['tag'],'seconds',obs['runtime_seconds'],'density',obs['operators']['density']['half'])


if __name__=='__main__':main()
