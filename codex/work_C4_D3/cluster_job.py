#!/usr/bin/env python3
"""One SLURM array task; no network, no cluster-specific Python dependencies."""
import argparse,json,os,time
from pathlib import Path
from d3 import Sector,measure,append
p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=int(os.environ.get('SLURM_ARRAY_TASK_ID','0')));p.add_argument('--manifest',default='cluster_manifest.json');p.add_argument('--dry-run',action='store_true');a=p.parse_args()
job=json.loads(Path(a.manifest).read_text())[a.index]
if a.dry_run:print(json.dumps(job,indent=2));raise SystemExit
out=Path('cluster_results');out.mkdir(exist_ok=True)
sec=Sector(job['family'],job['L'],target=job.get('charge'),cartan=job.get('cartan',False),cache_blocks=False)
result=measure(sec,job['coeff'],job['label'],levels=job.get('levels',12),entropy=job.get('entropy',True))
result['cluster_job']=job;result['slurm_job_id']=os.environ.get('SLURM_JOB_ID')
append(str(out/f'task_{a.index:04d}.jsonl'),result)
