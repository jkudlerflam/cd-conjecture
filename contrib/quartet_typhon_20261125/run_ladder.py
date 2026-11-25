"""Run isolated chi stages, preserving all states outside fetched results."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from common import HERE,source_hashes,write_json


def main(config):
    root=Path(os.environ['TYPHON_RESULTS']).parent
    results=root/'results';states=root/'states';states.mkdir(exist_ok=True)
    seconds=config.get('wall_seconds',3600)
    os.environ['JOB_END_EPOCH']=str(time.time()+seconds-30)
    write_json(results/'provenance.json',{'source_sha256':source_hashes(),'input':config,'job_id':os.environ.get('SLURM_JOB_ID')})
    if config.get('validate',True):
        subprocess.run([sys.executable,str(HERE/'validate.py'),'--output',str(results/'validation.json')],check=True)
    previous=config.get('input_state')
    outcomes=[]
    for chi in config['chis']:
        if float(os.environ['JOB_END_EPOCH'])-time.time()<config.get('minimum_stage_seconds',300):
            break
        c={k:v for k,v in config.items() if k not in ('chis','validate','wall_seconds')}
        c.update(chi=chi,results=str(results),input_state=previous,
                 output_state=str(states/f"L{c['L']}_p{c['parity']}_chi{chi}.h5"))
        c['measure']=chi in config.get('measure_chis',config['chis'])
        cp=results/f'config_chi{chi}.json';write_json(cp,c)
        subprocess.run([sys.executable,str(HERE/'run_stage.py'),str(cp)],check=True)
        summary=results/f"L{c['L']}_p{c['parity']}_chi{chi}_{c.get('tag','target')}_summary.json"
        row=json.loads(summary.read_text());outcomes.append(row)
        previous=c['output_state']
        if row['shelved'] or (not row['converged'] and config.get('stop_on_nonconvergence',False)):
            break
    write_json(results/'ladder.json',{'stages':outcomes,'requested_chis':config['chis'],
               'all_stages_converged':len(outcomes)==len(config['chis']) and all(x['converged'] for x in outcomes),
               'next_input_state':previous})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args()
    main(json.loads(Path(a.config).read_text()))
