"""Enumerate the non-pair margins of the centered references."""
import json
import numpy as np
from itertools import combinations
rows=json.load(open('centered_scan_results.json'));out=[]
for row in rows:
 n=row['n'];k=np.array(row['p']);mom=np.concatenate([-k,k])
 charges=np.concatenate([-np.eye(n,dtype=int),np.eye(n,dtype=int)])
 pairs=list(combinations(range(2*n),2));vals=[];resonant=0
 for a,b in combinations(pairs,2):
  m=charges[list(a)].sum(axis=0)-charges[list(b)].sum(axis=0)
  if not np.any(m):resonant+=1;continue
  vals.append(abs(np.angle(np.exp(1j*(mom[list(a)].sum()-mom[list(b)].sum())))))
 d2=min(abs(np.angle(np.exp(1j*(a-b)))) for a,b in combinations(mom,2))
 out.append(dict(n=n,delta2=float(d2),delta4_nonpair=float(min(vals)),resonant_pair_comparisons=resonant))
json.dump(out,open('centered_margins.json','w'),indent=2)
print(out)
