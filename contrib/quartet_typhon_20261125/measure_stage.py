"""Measure an immutable stage in a separate, bounded process."""
import argparse
import gzip
import json
from pathlib import Path
import time
from common import write_json
from run_stage import load,rss_mb
import meas_quartet,measure_matrix


def dump(path,obj):
    path=Path(path);tmp=path.with_name(path.name+'.tmp')
    with gzip.open(tmp,'wt') as f:json.dump(obj,f)
    tmp.replace(path)


def main(c):
    t0=time.time();out=Path(c['results']);L=c['L']
    stem=f"L{L}_p{c['parity']}_chi{c['chi']}_{c.get('tag','target')}"
    saved=load(c['output_state']);psi=saved['psi'];meta=json.loads(saved['metadata_json'])
    refs=c.get('refs',[L//2-17,L//2-8,L//2])
    path=out/(stem+'_quartet.json.gz')
    if not path.exists():
        m=meas_quartet.measure(psi,L,refs,meas_quartet.DS_DEFAULT)
        m.update(chi=int(max(psi.chi)),parity=c['parity'],h=c.get('h',1.783),E=meta['E'],
                 state_sweep_converged=meta['converged'])
        dump(path,m)
    if c.get('matrix'):
        path=out/(stem+'_matrix.json.gz')
        if not path.exists():
            m=measure_matrix.measure(psi,c.get('matrix_refs',[L//2-17]),rmax=c.get('matrix_rmax',96),
                                     centered=c.get('centered_matrix',False))
            m.update(parity=c['parity'],h=c.get('h',1.783),E=meta['E'],state_sweep_converged=meta['converged'])
            dump(path,m)
    write_json(out/(stem+'_measurement.json'),{'seconds':time.time()-t0,'max_rss_mb':rss_mb(),'complete':True})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args()
    main(json.loads(Path(a.config).read_text()))
