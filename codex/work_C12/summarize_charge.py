import json
rows=json.load(open('compressibility_check_results.json'))
out=[]
for r in rows:
 n0=r['L']//3;e={s['N']:s['E'] for s in r['sectors']}
 if n0 not in e:continue
 low=max((e[n0]-en)/(n0-n) for n,en in e.items() if n<n0)
 high=min((en-e[n0])/(n-n0) for n,en in e.items() if n>n0)
 out.append(dict(L=r['L'],delta=r['delta'],N=n0,mu_lower=low,mu_upper=high,stable_against_tested_sectors=low<high))
print(out)
with open('charge_intervals.json','w') as f:json.dump(out,f,indent=2)
