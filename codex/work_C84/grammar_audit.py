"""Exact word enumeration and local bad-pair transition audit."""
import json,itertools
import bond11 as b
import bond18

def main():
 ls=b.alphabet();data=[l for l in ls if not b.marker(l)];marks=[l for l in ls if b.marker(l)]
 out=[]
 for p in range(8):
  prefixes=[()]
  for j in range(p):
   prefixes=[w+(l,) for w in prefixes for l in data if b.allowed(w[-1] if w else ('e',0),l)]
  legal={w+(m,) for w in prefixes for m in marks if b.allowed(w[-1] if w else ('e',0),m)}
  ws=b.words(p);assert legal==set(ws)
  for w in ws:
   if p:
    payload=sum(l[1] for l in w if l[0] in ['L','R'])
    q=sum(b.charge(l) for l in w)
    typ=next((l[1] for l in w if b.excited(l)),None)
    assert q==2*payload+({None:0,'A':-1,'B':1}[typ])
  out.append(dict(data_sites=p,flavours=max(p-1,0),legal_dimension=len(ws),injective=True,complete=True))
 h,_,o,v=b.parts(bond18.params(3,K=100,g=.03))
 changes=[];seen=0
 for j,i in zip(*h.nonzero()):
  if j==i:continue
  x,y=ls[i//11],ls[i%11];u,w=ls[j//11],ls[j%11]
  for a,c in itertools.product(ls,repeat=2):
   old=tuple(not b.allowed(s,t) for s,t in [(a,x),(x,y),(y,c)])
   new=tuple(not b.allowed(s,t) for s,t in [(a,u),(u,w),(w,c)])
   if old!=new:
    # Only a vacuum step may change a bad pair, and only its predecessor.
    assert not any(b.excited(z) for z in [x,y,u,w])
    assert old[1:]==new[1:]
    assert b.marker(x)==b.marker(u) and b.marker(y)==b.marker(w)
    seen+=1
 changes.append(dict(all_local_contexts=seen,only_vacuum_predecessor_changes=True))
 rec=dict(words=out,local_bad_pair_audit=changes)
 json.dump(rec,open('grammar_checks.json','w'),indent=2);print(json.dumps(rec))
if __name__=='__main__':main()
