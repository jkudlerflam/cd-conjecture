import json
import fr_lib as F
from falsify import run
from pathlib import Path
p=Path('falsify_results.jsonl')
rows=list(map(json.loads,p.read_text().splitlines()))
invalid=[z for z in rows if z['name']=='XXZ-1.0']
Path('excluded_incomplete_mixture.jsonl').write_text(''.join(json.dumps(z)+'\n' for z in invalid))
p.write_text(''.join(json.dumps(z)+'\n' for z in rows if z['name']!='XXZ-1.0'))
with p.open('a') as f:
    for L in [4,6,8,10,12]:run('XXZ-1.0',F.xxz_bond(-1),L,f,L<=6)
