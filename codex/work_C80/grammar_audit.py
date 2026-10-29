"""Exact finite word audit for the shifted encoding, no diagonalization."""
import json
import bond_shift as b
import bond18 as parent

def main():
 out=[]
 for n in range(7):
  ls=b.alphabet(12);data=[l for l in ls if l[0] not in ['e','m']];marks=[l for l in ls if l[0] in ['e','m']]
  prefixes=[()]
  for j in range(n):
   prefixes=[w+(l,) for w in prefixes for l in data if b.allowed(w[-1] if w else ('e',0),l,12)]
  legal={w+(m,) for w in prefixes for m in marks if b.allowed(w[-1] if w else ('e',0),m,12)}
  image=b.words(n,12)
  assert legal==set(image) and len(image)==(5*n+3)*2**n
  for old,new in zip(parent.words(n),image):
   oldbits=sum((l[2] or 0) for l in old)
   newbits=sum(l[1] for l in new if l[0] in ['L','R','e'])
   assert oldbits==newbits
   assert sum(parent.charge(l) for l in old)==sum(b.charge(l) for l in new)
  out.append(dict(n=n,legal_dimension=len(legal),injective=True,complete=True,bit_number_preserved=True,charge_preserved=True))
 json.dump(out,open('grammar_checks.json','w'),indent=2);print(json.dumps(out))
if __name__=='__main__':main()
