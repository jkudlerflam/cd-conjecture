"""Bounded L24/chi4 I/O smoke test, one sweep, explicitly unconverged."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from common import HERE,write_json,source_hashes


def main():
    t=time.time()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d);(p/'results').mkdir()
        c={'L':24,'parity':0,'h':1.783,'chis':[4],'wall_seconds':600,'reserve_seconds':60,
           'max_sweeps':0,'validate':False,'refs':[0],'matrix_refs':[0],'matrix_rmax':12,
           'matrix':True,'tag':'workflow_test','seed':1}
        write_json(p/'config.json',c)
        env={**os.environ,'TYPHON_RESULTS':str(p/'results')}
        try:
            log=subprocess.check_output([sys.executable,str(HERE/'run_ladder.py'),str(p/'config.json')],env=env,text=True,stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as error:
            print(error.output)
            raise
        result=json.loads((p/'results/ladder.json').read_text());row=result['stages'][0]
        assert not result['all_stages_converged'] and not row['converged'] and row['measured']
        assert row['dE_last'] is None and row['sweeps']==1
        assert len(list((p/'states').glob('*.h5')))==1
        assert len(list((p/'results').glob('*quartet.json.gz')))==1
        assert len(list((p/'results').glob('*matrix.json.gz')))==1
        write_json(HERE/'workflow_validation.json',{'passed':True,'seconds':time.time()-t,
            'scope':'L24 chi4 one-sweep I/O test; not convergence evidence','converged':False,
            'separate_measurement_completed':True,'initial_delta_E_is_null':True,'source_sha256':source_hashes()})
        print('workflow passed',time.time()-t)


if __name__=='__main__':main()
