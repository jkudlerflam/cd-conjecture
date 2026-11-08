import json, sys
R=[json.loads(l) for f in sys.argv[1:] for l in open(f)]
R=[r for r in R if 'c' in r]
R.sort(key=lambda r:-max(r['c']))
for r in R[:20]: print(r['label'], [round(x,3) for x in r['p']], [round(x,3) for x in r['c']], [round(x,2) for x in r['gapL']])
