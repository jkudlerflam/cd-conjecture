"""A complete constant-control endpoint/color shuttle, conditional on marked outputs.
No input encoder or blank-to-program compiler is claimed.
All local maps are permutations on the full local Hilbert space.
"""
from dataclasses import dataclass
from itertools import product
import json

MARKS=('L','LI','LP','LC','RC','RP','RI','R','B')
# Each phase: starting mark, target mark, movement, swap quantum records, flag toggle.
FORWARD=(('L','LP',1,0,0),('LP','L',-1,1,0),
         ('L','LC',1,0,0),('LC','LI',-1,1,0),
         ('LI','RP',1,0,0),('RP','R',1,1,0),
         ('R','RC',-1,0,0),('RC','RI',1,1,0))
READOUT=(('RI','R',1,0,1),('R','L',-1,0,1),
         ('L','R',1,0,1),('R','L',-1,0,1),('L','RI',1,0,0))
REVERSE=tuple((t,s,-d,w,0) for s,t,d,w,_ in reversed(FORWARD))
PHASES=FORWARD+READOUT+REVERSE

class Shuttle:
    def __init__(self):
        self.symbols=tuple(product(MARKS,range(2)))
        self.index={s:i for i,s in enumerate(self.symbols)}
        self.B=len(self.symbols);self.Q=len(PHASES)
        self.partial={};used={}
        for k,(start,target,direction,swap,toggle) in enumerate(PHASES):
            assert target==PHASES[(k+1)%self.Q][0]
            for mark,flag in self.symbols:
                if mark==start:continue
                kp=(k+1)%self.Q if mark==target else k
                fp=flag^(toggle if mark==target else 0)
                src=k*self.B+self.index[mark,flag]
                dst=kp*self.B+self.index[mark,fp]
                assert src not in self.partial and dst not in used
                self.partial[src]=dst;used[dst]=src
        self.table=[-1]*(self.B*self.Q)
        for src,dst in self.partial.items():self.table[src]=dst
        for mark in MARKS:
            domain={i for i in range(len(self.table)) if self.symbols[i%self.B][0]==mark and i not in self.partial}
            image={i for i in range(len(self.table)) if self.symbols[i%self.B][0]==mark and i not in used}
            common=domain&image
            for i in common:self.table[i]=i
            for i,j in zip(sorted(domain-common),sorted(image-common)):self.table[i]=j
        assert sorted(self.table)==list(range(len(self.table)))
        self.inverse=[0]*len(self.table)
        for i,j in enumerate(self.table):self.inverse[j]=i

    def collide(self,k,mark,flag,inverse=False):
        v=(self.inverse if inverse else self.table)[k*self.B+self.index[mark,flag]]
        return v//self.B,self.symbols[v%self.B][1]

    def layout(self,n,a,b):
        assert 2<=a and a+1<b and b+1<=n-3
        marks=['B']*n
        for i,m in ((0,'L'),(1,'LI'),(a,'LP'),(a+1,'LC'),(b,'RC'),(b+1,'RP'),(n-2,'RI'),(n-1,'R')):marks[i]=m
        return marks

    def initial(self,n,a,b,records=None):
        return Configuration(self.layout(n,a,b),[0]*n,list(range(n)) if records is None else list(records),0,0,0)

    def step(self,c):
        assert c.phase==0
        start,target,d,swap,toggle=PHASES[c.control]
        j=c.position+d
        assert 0<=j<len(c.marks)
        if swap:c.records[c.position],c.records[j]=c.records[j],c.records[c.position]
        c.position=j
        src=c.control*self.B+self.index[c.marks[j],c.flags[j]]
        assert src in self.partial
        c.control,c.flags[j]=self.collide(c.control,c.marks[j],c.flags[j])

    def bond(self,c,j):
        """The full ordered-sweep G for one head, including inverse blocks.
        Record swaps are ordinary graded-space swaps (no fermionic sign).
        """
        d=PHASES[c.control][2];swap=PHASES[c.control][3]
        forward=c.phase==0 and ((d==1 and c.position==j) or (d==-1 and c.position==j+1))
        backward=c.phase==1 and ((d==1 and c.position==j+1) or (d==-1 and c.position==j))
        if forward or backward:
            if swap:c.records[j],c.records[j+1]=c.records[j+1],c.records[j]
            c.position=j+1 if c.position==j else j
            c.phase^=1

    def onsite(self,c,j):
        if c.position!=j:return
        c.control,c.flags[j]=self.collide(c.control,c.marks[j],c.flags[j],inverse=c.phase==0)
        c.phase^=1

    def manifest(self):
        return dict(scope='Marked-output shuttle only; no encoder or autonomous initializer',
                    static_marks=MARKS,phases=PHASES,ordinary_neutral_states=self.B,
                    direction_fixed_controls=self.Q,prescribed_collisions=len(self.partial),
                    complete_collision_entries=len(self.table),collision_permutation=self.table,
                    head_overlay_factor=(1+2*self.Q)*self.B,
                    outer_clock_factor=15*(1+2*self.Q)*self.B)

@dataclass
class Configuration:
    marks:list
    flags:list
    records:list
    position:int
    control:int
    phase:int=0
    def snapshot(self):return tuple(self.flags),tuple(self.records),self.position,self.control,self.phase

if __name__=='__main__':
    s=Shuttle()
    with open('endpoint_shuttle_table.json','w') as f:json.dump(s.manifest(),f,indent=2)
    print(json.dumps({k:v for k,v in s.manifest().items() if k not in ('collision_permutation','phases')},indent=2))
