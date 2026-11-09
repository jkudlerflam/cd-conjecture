import json,glob,numpy as np
def wavg(S,Lb):
    m=len(S)//2  # index of the middle cut (cut m+1 blocks)
    return float(np.mean(S[m-1:m+2]))
rows=[]
for f in sorted(glob.glob('bd_*.jsonl')):
    for l in open(f):
        r=json.loads(l); S=np.array(r['S'])
        rows.append((r['wa'],r['phi'],r['mu'],r.get('t2',-1.0),r['Lb'],wavg(S,r['Lb']),r.get('n_mid'),r['hist'][-1][2],f))
rows.sort()
for x in rows: print(x)
# two-size c from window averages
from collections import defaultdict
g=defaultdict(dict)
for x in rows: g[x[:4]][x[4]]=x[5]
print("two-size c = 6 [Sbar(L2)-Sbar(L1)]/ln(L2/L1):")
for k,v in g.items():
    Ls=sorted(v)
    for i in range(len(Ls)-1):
        a,b=Ls[i],Ls[i+1]
        print(k, f"{3*a}->{3*b}", round(6*(v[b]-v[a])/np.log(b/a),3))
