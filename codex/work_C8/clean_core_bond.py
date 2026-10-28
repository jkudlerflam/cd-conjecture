"""Sparse column oracle for the complete positive core bulk and boundary matrices.
No sector restriction is imposed on the supplied columns.
"""
from dataclasses import dataclass,replace
from collections import defaultdict
from clean_readout_machine import Machine,Config,ROLES,STATIC_ALLOWED
from sweeping_clock_checks import ALPHABET,PAIRS,LEFT,RIGHT,RULES

@dataclass(frozen=True)
class Site:
    clock:str
    role:str
    cursor:int
    data:int
    ready:int
    record:int
    head:int=0  # 0 ordinary, then (control,p,x,phase) in lexicographic order

class CoreBond:
    def __init__(self,q,J=1.,eta=1.,kclock=1.,krole=1.):
        assert q%4==0
        self.q=q;self.M=Machine();self.J=J;self.eta=eta;self.kclock=kclock;self.krole=krole
        self.overlay=1+8*self.M.Q;self.w=64*(q+1)*self.overlay;self.d=15*self.w
        self.cidx={c:i for i,c in enumerate(ALPHABET)};self.ridx={r:i for i,r in enumerate(ROLES)}
    def encode(self,s):
        z=self.cidx[s.clock]
        for value,radix in ((self.ridx[s.role],8),(s.cursor,2),(s.data,2),(s.ready,2),(s.record,self.q+1),(s.head,self.overlay)):z=z*radix+value
        return z
    def decode(self,z):
        vals=[]
        for radix in (self.overlay,self.q+1,2,2,2,8):z,v=divmod(z,radix);vals.append(v)
        head,record,ready,data,cursor,role=vals
        return Site(ALPHABET[z],ROLES[role],cursor,data,ready,record,head)
    def head_decode(self,h):
        z,phase=divmod(h-1,2);z,x=divmod(z,2);control,p=divmod(z,2)
        return self.M.names[control],p,x,phase
    def head_encode(self,control,p,x,phase):return 1+(((self.M.cindex[control]*2+p)*2+x)*2+phase)
    def config(self,sites):
        pos=next(j for j,s in enumerate(sites) if s.head)
        control,p,x,phase=self.head_decode(sites[pos].head)
        return Config([s.role for s in sites],[[s.cursor,s.data,s.ready,s.record] for s in sites],pos,control,p,x,phase)
    def fromconfig(self,c,sites):
        return tuple(Site(s.clock,s.role,*c.tape[j],self.head_encode(c.control,c.p,c.x,c.phase) if c.pos==j else 0) for j,s in enumerate(sites))
    def G(self,a,b):
        if bool(a.head)+bool(b.head)!=1:return a,b
        c=self.config((a,b));self.M.bond(c,0,self.q);return self.fromconfig(c,(a,b))
    def F(self,a):
        if not a.head:return a
        c=self.config((a,));self.M.onsite(c,0,self.q);return self.fromconfig(c,(a,))[0]
    def charge(self,a):
        return a.data+int(self.q//2<=a.record<self.q)+(self.head_decode(a.head)[2] if a.head else 0)
    def good_input(self,a):
        return (a.cursor==int(a.role=='S') and a.ready==0 and
                (a.role in ('S','D') or a.data==0) and
                a.record==(0 if a.role in ('A','B','E','F') else self.q) and
                a.head==(self.head_encode('L0',0,0,0) if a.role=='L' else 0))
    def onsite(self,a):
        out=defaultdict(float)
        if a.clock in ('C0','R1'):
            out[a]+=self.J/2;out[replace(a,clock='R1' if a.clock=='C0' else 'C0')]-=self.J/2
        if a.clock=='I' and not self.good_input(a):out[a]+=self.eta
        return out
    def bulk(self,i,j):
        a,b=self.decode(i),self.decode(j);out=defaultdict(float)
        out[(a,b)]+=self.kclock*int((a.clock,b.clock) not in PAIRS)+self.krole*int((a.role,b.role) not in STATIC_ALLOWED)
        for source,target,label in RULES:
            if (a.clock,b.clock) not in (source,target):continue
            dest=target if (a.clock,b.clock)==source else source
            aa,bb=a,b
            if label=='out':aa,bb=self.G(a,b)
            elif label=='back':bb=self.F(b)
            out[(a,b)]+=self.J/2
            out[(replace(aa,clock=dest[0]),replace(bb,clock=dest[1]))]-=self.J/2
        for aa,v in self.onsite(a).items():out[(aa,b)]+=v/2
        for bb,v in self.onsite(b).items():out[(a,bb)]+=v/2
        return {(self.encode(aa),self.encode(bb)):v for (aa,bb),v in out.items() if v}
    def boundary(self,i,side):
        a=self.decode(i);out=defaultdict(float)
        out[a]+=self.kclock*int(a.clock not in (LEFT if side=='L' else RIGHT))+self.krole*int(a.role!=side)
        for aa,v in self.onsite(a).items():out[aa]+=v/2
        edges=[]
        if side=='L':
            for bit in ('0','1'):edges.extend([('R'+bit,'D'+bit,False),('E'+bit,'C'+bit,True)])
        else:
            edges.append(('I','R0',False))
            for bit in ('0','1'):edges.append(('D'+bit,'E'+bit,False))
        for src,dst,F in edges:
            if a.clock not in (src,dst):continue
            aa=self.F(a) if F else a
            aa=replace(aa,clock=dst if a.clock==src else src)
            out[a]+=self.J/2;out[aa]-=self.J/2
        return {self.encode(aa):v for aa,v in out.items() if v}

if __name__=='__main__':
    import json
    import numpy as np
    rng=np.random.default_rng(2026);rows=[]
    for q in (4,8):
        H=CoreBond(q);nonzeros=0;err=0.
        for _ in range(1600):
            i,j=map(int,rng.integers(0,H.d,2));a,b=H.decode(i),H.decode(j)
            assert H.encode(a)==i and H.encode(b)==j
            for (k,l),value in H.bulk(i,j).items():
                err=max(err,abs(H.bulk(k,l).get((i,j),0)-value))
                assert H.charge(a)+H.charge(b)==H.charge(H.decode(k))+H.charge(H.decode(l));nonzeros+=1
            for side in ('L','R'):
                for k,value in H.boundary(i,side).items():
                    err=max(err,abs(H.boundary(k,side).get(i,0)-value))
                    assert H.charge(a)==H.charge(H.decode(k));nonzeros+=1
        assert err==0
        rows.append(dict(q=q,core_site_dimension=H.d,sampled_bulk_columns=1600,checked_entries=nonzeros,hermiticity_error=err))
    open('clean_core_bond_checks.json','w').write(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))
