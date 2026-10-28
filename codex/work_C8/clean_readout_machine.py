"""Autonomous graded, two-port Clifford/color readout machine.
No program tape. Logical cursor and q-state counters live on the tape.
Whole-space move gates and completed onsite collisions are explicit.
"""
from dataclasses import dataclass
from itertools import product
import json

ROLES=('L','A','B','S','D','E','F','R')
BITS=tuple(product(range(2),repeat=3))  # cursor, data, READY
STATIC_ALLOWED={('L','A'),('A','B'),('B','S'),('S','D'),('S','E'),('D','D'),('D','E'),('E','F'),('F','R')}

@dataclass(frozen=True)
class Control:
    name:str
    direction:int
    gate:str='I'
    inverse:bool=False

@dataclass
class Config:
    roles:list
    tape:list  # [cursor,data,READY,record]
    pos:int
    control:str
    p:int=0
    x:int=0
    phase:int=0
    def snapshot(self):return tuple(tuple(t) for t in self.tape),self.pos,self.control,self.p,self.x,self.phase

# Local ops act on (data,READY,record,p,x); cursor changes are classical rules.
def local_op(op,z,q,inverse=False):
    data,ready,v,p,x=z;s=q//2;sgn=-1 if inverse else 1
    if op=='I':pass
    elif op=='PARITY':p^=data
    elif op=='TAKE':data,x=x,data
    elif op=='LOAD':
        if (x,v)==(0,q):x,v=0,0
        elif (x,v)==(0,0):x,v=0,q
        elif (x,v)==(1,q):x,v=0,s
        elif (x,v)==(0,s):x,v=1,q
    elif op=='ODD_P':
        if s<=v<q:v=s+(v-s+sgn*p)%s
    elif op=='EVEN_P':
        if v<s:v=(v+sgn*p)%s
    elif op=='FLAG':ready^=1
    else:raise ValueError(op)
    return data,ready,v,p,x


def record_gate(gate,a,b,q,inverse=False):
    """a=departure record, b=arrival record. Explicit full-space permutations."""
    s=q//2;sgn=-1 if inverse else 1
    if gate=='I':pass
    elif gate=='SWAP':a,b=b,a
    elif gate=='EVEN_ADD':
        if a<s and b<s:a=(a+sgn*b)%s
    elif gate=='ODD_ADD':
        if s<=a<q and b<s:a=s+(a-s+sgn*b)%s
    elif gate=='SUM':
        if a<s and b<s:b=(b+sgn*a)%s
    elif gate=='INC':
        if a<s and b<s:
            if inverse:
                old=(b-2)%s;a=(a-int(old+2>=s))%s;b=old
            else:a=(a+int(b+2>=s))%s;b=(b+2)%s
    else:raise ValueError(gate)
    return a,b

class Machine:
    def __init__(self):
        self.controls={};self.rules={};self.ranges=set()
        self._build()
        self.names=tuple(self.controls);self.cindex={k:i for i,k in enumerate(self.names)}
        self.Q=len(self.names)
        # C acts on control,static role,cursor, plus the local quantum factors.
        # Each prescribed block has a permutation quantum gate. Complete the
        # neutral block injection within each immutable static role.
        self.neutral=[(c,r,k) for c in self.names for r in ROLES for k in range(2)]
        self.nindex={v:i for i,v in enumerate(self.neutral)}
        self.table=[None]*len(self.neutral)
        for src,(dst,op,inv) in self.rules.items():self.table[self.nindex[src]]=(self.nindex[dst],op,inv)
        for role in ROLES:
            dom={i for i,z in enumerate(self.neutral) if z[1]==role and z not in self.rules}
            ran={i for i,z in enumerate(self.neutral) if z[1]==role and z not in self.ranges}
            common=dom&ran
            for i in common:self.table[i]=(i,'I',False)
            assert len(dom)==len(ran)
            for i,j in zip(sorted(dom-common),sorted(ran-common)):self.table[i]=(j,'I',False)
        assert sorted(z[0] for z in self.table)==list(range(len(self.table)))
        self.inverse=[None]*len(self.table)
        for i,(j,op,inv) in enumerate(self.table):self.inverse[j]=(i,op,not inv)

    def control(self,name,d,gate='I',inverse=False):
        assert name not in self.controls
        self.controls[name]=Control(name,d,gate,inverse)

    def rule(self,c,r,k,cp,kp=None,op='I',inverse=False):
        if kp is None:kp=k
        src=(c,r,k);dst=(cp,r,kp)
        assert src not in self.rules,('domain',src)
        assert dst not in self.ranges,('range',src,dst)
        self.rules[src]=(dst,op,inverse);self.ranges.add(dst)

    @staticmethod
    def is_cursor(role,k):return role in ('S','D') and k==1

    def segment(self,c,start,end,next_c,scan='I',endop='I',reverse=False):
        """start/end predicates on (role,cursor), no semantic rule at start."""
        for role,k in product(ROLES,range(2)):
            if start(role,k):continue
            if end(role,k):self.rule(c,role,k,next_c,op=endop,inverse=reverse)
            else:self.rule(c,role,k,c,op=scan if role in ('S','D') else 'I',inverse=reverse)

    def _build(self):
        cur=self.is_cursor
        for side,d,port,near,far,other in (('L',1,'L','A','B','R'),('R',-1,'R','F','E','L')):
            pref=lambda i:f'{side}{i}'
            # Each forward phase moves from start to end; operation on arrival.
            segs=[(port,'CUR',d,'I','PARITY','TAKE'),
                  ('CUR',port,-d,'I','I','LOAD'),
                  (port,near,d,'EVEN_ADD','I','I'),
                  (near,far,d,'SWAP','I','I'),
                  (far,port,-d,'I','I','I'),
                  (port,near,d,'ODD_ADD','I','I'),
                  (near,far,d,'SWAP','I','I'),
                  (far,port,-d,'I','I','ODD_P'),
                  (port,far,d,'I','I','I'),
                  (far,near,-d,'SUM','I','EVEN_P'),
                  (near,port,-d,'I','I','I')]
            pred=lambda x:cur if x=='CUR' else (lambda r,k:r==x)
            for i,(start,end,di,gate,scan,endop) in enumerate(segs):
                self.control(pref(i),di,gate)
                # Reverse move must apply the inverse gate with swapped
                # departure/arrival convention, encoded by REV_ prefix.
                self.control('~'+pref(i),-di,'REV_'+gate,True)
            self.control(side+'OUT',d)
            self.control(side+'BACK',-d)
            # Build forward rules except final arrival, then exact reverse
            # collision blocks with opposite movement controls.
            saved=[]
            for i,(start,end,di,gate,scan,endop) in enumerate(segs):
                c=pref(i)
                for role,k in product(ROLES,range(2)):
                    if pred(start)(role,k):continue
                    if pred(end)(role,k):
                        if i==10:continue
                        cp=pref(i+1);op=endop
                    else:cp=c;op=scan if role in ('S','D') else 'I'
                    self.rule(c,role,k,cp,op=op)
                    saved.append((c,role,k,cp,op))
            for c,role,k,cp,op in saved:self.rule('~'+cp,role,k,'~'+c,op=op,inverse=True)
            # Readout: the last preparatory collision is identity on records.
            for k in range(2):self.rule(pref(10),port,k,side+'OUT',op='FLAG')
            self.segment(side+'OUT',pred(port),pred(other),side+'BACK')
            self.segment(side+'BACK',pred(other),pred(port),'~'+pref(10),endop='FLAG')
            # Initial reverse collision, after the logical input is restored.
            nextc='TO_R' if side=='L' else 'R_NEAR'
            for k in range(2):self.rule('~'+pref(0),port,k,nextc)
        # Travel and counter updates. All data-independent.
        extra=(('TO_R',1,'I'),('R_NEAR',-1,'I'),('R_INC',-1,'INC'),
               ('TO_L',-1,'I'),('L_NEAR',1,'I'),('L_INC',1,'INC'),
               ('FIND',1,'I'),('NEXT',1,'I'),('RETURN',-1,'I'),('RESET',-1,'I'))
        for c,d,g in extra:self.control(c,d,g)
        p=lambda v:lambda r,k:r==v
        self.segment('TO_R',p('L'),p('R'),'R0')
        self.segment('R_NEAR',p('R'),p('F'),'R_INC')
        self.segment('R_INC',p('F'),p('E'),'TO_L')
        self.segment('TO_L',p('E'),p('L'),'L_NEAR')
        self.segment('L_NEAR',p('L'),p('A'),'L_INC')
        self.segment('L_INC',p('A'),p('B'),'FIND')
        for role,k in product(ROLES,range(2)):
            if role=='B':continue
            if cur(role,k):self.rule('FIND',role,k,'NEXT',0)
            else:self.rule('FIND',role,k,'FIND')
        # NEXT is one move, so only the following initialized cases are needed.
        self.rule('NEXT','D',0,'RETURN',1)
        self.rule('NEXT','E',0,'RESET',0)
        for role,k in product(ROLES,range(2)):
            if cur(role,k):continue
            if role=='L':self.rule('RETURN',role,k,'L0')
            else:self.rule('RETURN',role,k,'RETURN')
        for role,k in product(ROLES,range(2)):
            if role=='E':continue
            if role=='S' and k==0:self.rule('RESET',role,k,'RETURN',1)
            elif role=='S':continue
            else:self.rule('RESET',role,k,'RESET')

    def initialize(self,q,bits):
        assert q%4==0 and len(bits)>=1
        roles=['L','A','B','S']+['D']*(len(bits)-1)+['E','F','R']
        tape=[[int(r=='S'),0,0,0 if r in ('A','B','E','F') else q] for r in roles]
        for j,bit in enumerate(bits):tape[j+3][1]=bit
        return Config(roles,tape,0,'L0')

    def collide(self,c,q,inverse=False,require=True):
        src=(c.control,c.roles[c.pos],c.tape[c.pos][0])
        if require and not inverse:assert src in self.rules,('unprescribed',src,c.snapshot())
        table=self.inverse if inverse else self.table
        dst,op,inv=table[self.nindex[src]]
        control,role,cursor=self.neutral[dst]
        t=c.tape[c.pos]
        t[1],t[2],t[3],c.p,c.x=local_op(op,(t[1],t[2],t[3],c.p,c.x),q,inv)
        t[0]=cursor;c.control=control

    def move(self,c,q):
        ctl=self.controls[c.control];j=c.pos+ctl.direction
        assert 0<=j<len(c.roles),('exit',c.control,c.pos)
        a,b=c.tape[c.pos][3],c.tape[j][3]
        if ctl.gate.startswith('REV_'):
            b,a=record_gate(ctl.gate[4:],b,a,q,ctl.inverse)
        else:a,b=record_gate(ctl.gate,a,b,q,ctl.inverse)
        c.tape[c.pos][3],c.tape[j][3]=a,b;c.pos=j

    def step(self,c,q):self.move(c,q);self.collide(c,q)

    def bond(self,c,j,q):
        ctl=self.controls[c.control];d=ctl.direction
        forward=c.phase==0 and ((d==1 and c.pos==j) or (d==-1 and c.pos==j+1))
        backward=c.phase==1 and ((d==1 and c.pos==j+1) or (d==-1 and c.pos==j))
        if forward:self.move(c,q);c.phase=1
        elif backward:
            # Adjoint of the forward move, at the same fixed control.
            src=c.pos-d;a,b=c.tape[src][3],c.tape[c.pos][3]
            if ctl.gate.startswith('REV_'):b,a=record_gate(ctl.gate[4:],b,a,q,not ctl.inverse)
            else:a,b=record_gate(ctl.gate,a,b,q,not ctl.inverse)
            c.tape[src][3],c.tape[c.pos][3]=a,b;c.pos=src;c.phase=0

    def onsite(self,c,j,q):
        if c.pos==j:
            self.collide(c,q,inverse=c.phase==0,require=False);c.phase^=1

    def manifest(self):
        # Ordinary: role,cursor,data,READY,record. Head: p,x and fresh/updated.
        factor=8*2*2*2*(1+8*self.Q)
        return dict(scope='Autonomous clean-readout core controller; no CFT claim',
                    controls=[vars(v) for v in self.controls.values()],Q=self.Q,
                    ordinary_factor=64,head_bit_factor=4,
                    work_factor=factor,clock_factor=15*factor,
                    prescribed_blocks=len(self.rules),complete_blocks=len(self.table),
                    neutral_labels=self.neutral,collision_blocks=self.table,
                    static_roles=ROLES,static_allowed=sorted(STATIC_ALLOWED))

if __name__=='__main__':
    m=Machine();manifest=m.manifest()
    with open('clean_readout_table.json','w') as f:json.dump(manifest,f,indent=2)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('controls','neutral_labels','collision_blocks')},indent=2))
    c=m.initialize(4,[1,0]);initial=c.snapshot()
    for i in range(10000):
        m.step(c,4)
        if c.snapshot()==initial:print('cycle',i+1);break
    else:raise AssertionError('No cycle')
