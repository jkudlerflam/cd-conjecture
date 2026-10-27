import json
from d3 import Sector,measure,append
from pathlib import Path
for family,source,out in [('sp4','refine_sp4.jsonl','sp4_L4.jsonl'),('orbital','scan_orbital_dense.jsonl','orbital_L4.jsonl')]:
 rows=[json.loads(x) for x in open(source)];seen=set();sec=Sector(family,4)
 for r in rows:
  if r['label'] in seen:continue
  seen.add(r['label']);append(out,measure(sec,r['coeff'],r['label'],8))
