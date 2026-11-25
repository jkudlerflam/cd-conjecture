"""Read-only final pilot diagnostics; never optimize or overwrite a saved MPS."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from common import write_json,source_hashes,sha256
from run_stage import load
from measure_stage import main as measure


def main(config):
    out=Path(os.environ['TYPHON_RESULTS'])
    write_json(out/'provenance.json',{'source_sha256':source_hashes(),'input':config,'job_id':os.environ.get('SLURM_JOB_ID')})
    cases=[]
    for case in config['cases']:
        state=case['state'];saved=load(state);meta=json.loads(saved['metadata_json'])
        assert meta['delta']==3/16 and meta['L']==128 and meta['parity']==0
        assert meta['h']==(1.783 if case['tag']=='target' else 1.772)
        psi=saved['psi']
        parity=float(np.real(psi.expectation_value_term([('Sigmaz',i) for i in range(psi.L)])))
        norm_error=float(np.linalg.norm(psi.norm_test()))
        assert abs(parity-1)<1e-9
        del psi,saved
        c={k:meta[k] for k in ('L','parity','chi','h')}
        c.update(output_state=state,results=str(out),tag=case['tag'],matrix=True,centered_matrix=True,
                 refs=case.get('refs',[47,56,64]),matrix_refs=case.get('matrix_refs',[47,56]),matrix_rmax=64)
        write_json(out/(case['tag']+'_config.json'),c)
        measure(c)
        cases.append({'tag':case['tag'],'state':state,'state_sha256':sha256(state),
                      'stage_metadata':meta,'measurement_complete':True,'parity_expectation':parity,
                      'saved_state_norm_error':norm_error})
        write_json(out/'measured_states.json',cases)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args()
    main(json.loads(Path(a.config).read_text()))
