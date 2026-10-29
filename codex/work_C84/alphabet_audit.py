"""Bounded quotient audit and exact additive-charge linear algebra."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,itertools
from fractions import Fraction as Q
import numpy as np
import bond11 as b
import bond18

def rank(rows):
 A=[[Q(x) for x in row] for row in rows];r=0
 for j in range(len(A[0])):
  k=next((k for k in range(r,len(A)) if A[k][j]),None)
  if k is None:continue
  A[k],A[r]=A[r],A[k];v=A[r][j];A[r]=[x/v for x in A[r]]
  for k in range(r+1,len(A)):
   v=A[k][j]
   if v:A[k]=[x-v*y for x,y in zip(A[k],A[r])]
  r+=1
  if r==len(A):break
 return r

def main():
 p=bond18.params(3,D=4,Delta=.5,K=100,g=.04,td=1e-4)
 ls=b.alphabet();h,_,o,v=b.parts(p);rows=[]
 for mat in [h,v]:
  for j,i in zip(*mat.nonzero()):
   if j==i:continue
   row=[0]*11
   for k in [j//11,j%11]:row[k]+=1
   for k in [i//11,i%11]:row[k]-=1
   rows.append(row)
 r=rank(rows);assert r==8
 basis=[[1]*11,[int(b.marker(l)) for l in ls],[b.charge(l) for l in ls]]
 assert rank(basis)==3 and all(sum(x*y for x,y in zip(row,z))==0 for row in rows for z in basis)
 out=dict(additive_charge_rank=r,additive_charge_dimension=11-r,charge_basis=['identity','marker','payload_charge'],single_letter_merges=[])
 ws=b.words(3)
 for a,c in itertools.combinations(ls,2):
  if b.charge(a)!=b.charge(c):continue
  f=lambda l:a if l==c else l
  new=list(dict.fromkeys(map(f,ls)));ix={l:i for i,l in enumerate(new)};d=len(new)
  images=[tuple(map(f,w)) for w in ws]
  rec=dict(letters=[a,c],parent_dimension=len(ws),image_dimension=len(set(images)),marker_data_merge=b.marker(a)!=b.marker(c))
  adj={l:set() for l in new}
  for mat in [h,v]:
   for j,i in zip(*mat.nonzero()):
    for oldidx,newidx in [(i//11,j//11),(i%11,j%11)]:
     x,y=f(ls[oldidx]),f(ls[newidx]);adj[x].add(y);adj[y].add(x)
  unseen=set(new);components=0
  while unseen:
   todo=[unseen.pop()];components+=1
   while todo:
    for y in adj[todo.pop()] & unseen:unseen.remove(y);todo.append(y)
  rec['single_site_transition_components']=components
  if rec['marker_data_merge']:assert components==1
  if a[0] in ['L','R'] and c[0] in ['L','R']:
   qh=np.zeros((d*d,d*d),complex);seen={}
   for j,i in zip(*h.nonzero()):
    key=(ix[f(ls[j//11])]*d+ix[f(ls[j%11])],ix[f(ls[i//11])]*d+ix[f(ls[i%11])])
    assert key not in seen or abs(seen[key]-h[j,i])<1e-12
    seen[key]=h[j,i]
   for key,val in seen.items():qh[key]=val
   def allowed(x,y):return any(b.allowed(l,r) for l in ls for r in ls if f(l)==x and f(r)==y)
   data=[l for l in new if l[0] in ['L','R']]
   words=[w+(('e',0),) for w in itertools.product(data,repeat=3)];idx={w:i for i,w in enumerate(words)}
   H=np.zeros((len(words),len(words)),complex);fixed=('e',0)
   for col,w in enumerate(words):
    full=(fixed,)+w
    for j in range(4):
     H[col,col]+=100*(not allowed(full[j],full[j+1]))
     ci=ix[full[j]]*d+ix[full[j+1]]
     for ri in np.flatnonzero(qh[:,ci]):
      z=list(full);z[j]=new[ri//d];z[j+1]=new[ri%d]
      if z[0]==fixed and tuple(z[1:]) in idx:H[idx[tuple(z[1:])],col]+=qh[ri,ci]
   rec['vacuum_dimension']=len(H);rec['lowest_four_energies']=np.linalg.eigvalsh(H)[:4].tolist()
   rec['parent_vacuum_energy']=-2*np.cos(np.pi/5)
   legal=[i for i,w in enumerate(words) if all(allowed(x,y) for x,y in zip((fixed,)+w,w))]
   rec['legal_vacuum_dimension']=len(legal)
   rec['lowest_four_energies_infinite_K']=np.linalg.eigvalsh(H[np.ix_(legal,legal)])[:4].tolist()
   HL=H[np.ix_(legal,legal)];ee,uu=np.linalg.eigh(HL)
   A=np.rint(HL.real).astype(int);assert np.array_equal(A,HL)
   w=np.rint(10000*uu[:,0].real).astype(int)
   numerator=int(w@A@w);denominator=int(w@w)
   assert numerator < -2*denominator
   rec['exact_rayleigh_witness']=dict(vector=w.tolist(),numerator=numerator,denominator=denominator,strictly_below_minus_two=True)
  out['single_letter_merges'].append(rec)
 assert len(out['single_letter_merges'])==10
 json.dump(out,open('alphabet_checks.json','w'),indent=2);print(json.dumps(out))
if __name__=='__main__':main()
