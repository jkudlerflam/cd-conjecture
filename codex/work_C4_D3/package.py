#!/usr/bin/env python3
import hashlib,json,sys,zipfile
from pathlib import Path
from datetime import datetime,timezone
import numpy,scipy
from analyze import FILES
codes=['d3.py','extra_scan.py','analyze.py','refinement.py','refine_sp4.py','global_velocity.py','catalogue.py','plots.py','sp4_followup.py','sp4_towers.py','sp4_tower_analysis.py','cluster_job.py','make_cluster_manifest.py','near_saturation.py','diagnostics.py','u1_occupations.py','crosschecks.py','verify_outputs.py','fix_charge.py','small_sizes.py','package.py']
artifacts=['REPORT.md','METHODS.md','README.md','D3_STATE.md','requirements.txt','typhon.sbatch','cluster_manifest.json','catalogue.csv','calibrations_summary.json','campaign_counts.json','selftest.json','crosschecks.json','verification.json','survey.png','survey.svg','sp4_towers.jsonl','sp4_two_velocity_test.json','sp4_crossings.jsonl','u1_occupations.jsonl','charge_audit.jsonl','u1_global_velocity.jsonl','refined_global_velocity.jsonl','deformation_global_velocity.jsonl','deformation_velocity_input.jsonl','invalid_spin1_L9.jsonl']
paths=sorted(set(codes+artifacts+FILES))
start=datetime.fromisoformat('2026-10-21T20:43:23+00:00');end=datetime.now(timezone.utc)
meta={'start_utc':start.isoformat(),'finish_utc':end.isoformat(),'wall_minutes':(end-start).total_seconds()/60,'local_limit_minutes':75,'cluster_jobs_submitted':0,'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,'closed_book':'No pre-existing project files or external references were read. All files in this manifest were created in this campaign.','retained_measurement_records':3526,'quarantined_wrong_sector_file':'invalid_spin1_L9.jsonl','sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}}
assert meta['wall_minutes']<75
Path('provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
with zipfile.ZipFile('d3_deliverables.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in paths+['provenance.json']:z.write(p,p)
print(json.dumps({'wall_minutes':meta['wall_minutes'],'packaged_files':len(paths)+1,'zip_bytes':Path('d3_deliverables.zip').stat().st_size}))
