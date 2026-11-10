"""Write or verify the final artifact manifest; local states remain Git-ignored."""
import argparse,datetime,hashlib,json,platform,subprocess
from collections import Counter
from pathlib import Path
import numpy,scipy,tenpy
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
target=HERE/'manifest.json'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
if a.verify:
    d=json.loads(target.read_text());errors=[]
    for name,want in d['files'].items():
        f=HERE/name
        if not f.exists() or sha(f)!=want:errors.append(name)
    for r in d['local_checkpoints']:
        f=HERE/r['path']
        if not f.exists() or sha(f)!=r['sha256']:errors.append(r['path'])
    print(json.dumps(dict(verified=not errors,errors=errors)))
    raise SystemExit(bool(errors))
records=[json.loads(f.read_text()) for f in (HERE/'results').glob('*.json')]
assert all(r['status']!='running' for r in records),'Workers still active'
audit=json.loads((HERE/'audit_report.json').read_text());assert not audit['errors']
counts=dict(Counter(r['status'] for r in records));assert counts==audit['counts'],'Run final audit again'
vumps=[json.loads(f.read_text()) for f in (HERE/'vumps_results').glob('*.json')]
assert all(r['status']!='running' for r in vumps),'VUMPS checks still active'
vcounts=dict(Counter(r['status'] for r in vumps));assert vcounts==audit['vumps_counts'],'Run final audit again'
files={str(f.relative_to(HERE)):sha(f) for f in sorted(HERE.rglob('*')) if f.is_file()
       and '__pycache__' not in f.parts and 'checkpoints_local' not in f.parts
       and f!=target and not f.name.endswith('.tmp')}
checkpoints=[dict(path=str(f.relative_to(HERE)),sha256=sha(f),bytes=f.stat().st_size)
             for f in sorted((HERE/'checkpoints_local').glob('*.pkl'))]
d=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
       git_parent_at_manifest=subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
       input_commit_original='7014ef41',input_commit_after_shared_rebase='367f2e69',
       versions=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,tenpy=tenpy.__version__),
       result_status_counts=counts,accepted_records=sum(r.get('converged',False) for r in records),
       vumps_status_counts=vcounts,accepted_vumps_records=sum(r.get('converged',False) for r in vumps),
       files=files,local_checkpoints=checkpoints,
       checkpoint_policy='Workspace-only pickle states; Git-ignored. Checksums and metadata are committed; no remote computation or publication.')
target.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(dict(files=len(files),checkpoints=len(checkpoints),counts=counts)))
