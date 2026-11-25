"""Bounded follow-up on saved states with every Hermitian probe centred at i+1/2."""
import argparse
import gzip
import json
import os
from pathlib import Path
import resource
import time
import h5py
from tenpy.tools import hdf5_io
from common import write_json,source_hashes,sha256
import measure_matrix


def main(config):
    out=Path(os.environ['TYPHON_RESULTS'])
    write_json(out/'provenance.json',{'source_sha256':source_hashes(),'input':config,'job_id':os.environ.get('SLURM_JOB_ID')})
    rows=[]
    for case in config['cases']:
        t=time.time();state=Path(case['state'])
        with h5py.File(state,'r') as f:saved=hdf5_io.load_from_hdf5(f)
        psi=saved['psi'];meta=json.loads(saved['metadata_json'])
        m=measure_matrix.measure(psi,case['refs'],rmax=case.get('rmax',64),centered=True)
        m.update(parity=meta['parity'],h=meta['h'],E=meta['E'],state_sha256=sha256(state))
        path=out/(case['label']+'_matrix.json.gz');tmp=path.with_name(path.name+'.tmp')
        with gzip.open(tmp,'wt') as f:json.dump(m,f)
        tmp.replace(path)
        rows.append({'label':case['label'],'L':psi.L,'chi':int(max(psi.chi)),'h':meta['h'],
                     'seconds':time.time()-t,'max_rss_mb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                     'state_sha256':m['state_sha256']})
        write_json(out/'measurements.json',rows)
        del psi,saved


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args()
    main(json.loads(Path(a.config).read_text()))
