"""Sparse columns of the complete C8 periodic NN bond.

The on-site basis is selector x (three typed vacancies or three occupied
tags carrying core x filter) x two equalizer qubits. No promised sector
is imposed. Preparation of normalizations is intended for small q only.
For general q the same coefficients are finite spectral definitions.
"""
from collections import defaultdict
from dataclasses import dataclass,replace
import numpy as np
from clean_core_bond import CoreBond
from clean_readout_counts import counts
from length_filter_checks import filter_data,plus
from vacancy_parent import charge_hamiltonian
from rational_target import target


@dataclass(frozen=True)
class PhysicalSite:
    selector:int              # #, L, X, R = 0,1,2,3
    wrapper:int               # H0,H1,H2, then (tag,core,filter)
    A:int
    B:int


def qubit_bond(phi):
    v=np.array([0,1,-np.exp(-1j*phi),0])/np.sqrt(2)
    return np.outer(v,v.conj())


def compensator(states,values):
    u,v=states
    chi=1-abs(np.vdot(u,v))**2
    return (values[0]*(np.eye(2)-np.outer(v,v.conj()))+
            values[1]*(np.eye(2)-np.outer(u,u.conj())))/chi


class PeriodicBond:
    selector_arrows=((1,0),(0,0),(0,1),(1,1))
    wrapper_arrows=((0,0),(1,1),(2,2),(0,1),(1,1),(1,2))

    def __init__(self,q=4,U=1.,vacancy_hop=None,J0=1.):
        self.core=CoreBond(q)
        self.q=q;self.m=q*q//8;self.k=self.m+6;self.b=self.k+1
        self.Q=self.core.d;self.s=6*self.Q+3;self.d=16*self.s
        self.U=U;self.w=.9*U/(self.k+1) if vacancy_hop is None else vacancy_hop
        assert 0<self.w<U/(self.k+1)
        self.J0=J0
        phase,hf,af,bf,eps=filter_data(self.k)
        self.filter_phase=phase;self.filter_scale=2*U/eps
        self.hf=hf*self.filter_scale;self.af=af*self.filter_scale;self.bf=bf*self.filter_scale
        self.phi=1/self.b;self.he=qubit_bond(self.phi)
        self.pin=np.eye(2)-np.outer(plus(0),plus(0).conj())
        self.holes=[];self.chis=[];self.charge_energies=[];self.r=[];self.s2=[]
        for p in (1,2):
            hs,hc,chi,ev=charge_hamiltonian(self.k,p,U,self.w)
            self.holes.append(hs);self.chis.append(chi)
            self.charge_energies.append(float(chi@hc@chi))
            self.r.append(sum(float(z*z) for h,z in zip(hs,chi) if 0 not in h))
            self.s2.append(sum(float(z*z) for h,z in zip(hs,chi) if 0 not in h and 1 not in h))
        a=lambda l:4*np.sin(np.pi/(2*l))**2
        b=self.b;ec=a(b+1)-a(b);muc=2-((b+1)*a(b)-b*a(b+1))
        eta=min((a(b-1)-2*a(b)+a(b+1))/2,(a(b)-2*a(b+1)+a(b+2))/3)
        self.t=2*U*(1/eta+b/muc+1/abs(ec))
        e1,e2=self.charge_energies
        self.marker_mu=self.t*muc+4*U-(b+1)*e1+b*e2
        self.base_density=self.t*ec+e2-e1
        self.Ksel=4*self.t+9*U
        self.states=[plus(-(b-1)*self.phi),plus(-b*self.phi)]
        self.ML=compensator(self.states,[1/x for x in self.r]);self.MR=self.ML.conj()
        self.MD=compensator(self.states,[1/x for x in self.s2])
        lookup={h:i for i,h in enumerate(self.holes[1])}
        self.beta_wrap=sum(z*self.chis[1][lookup[(h[0],self.k+1)]]
                           for h,z in zip(self.holes[0],self.chis[0]))**2
        j=np.arange(1,b)
        self.head_overlap=float(np.sum(2*np.sin(np.pi*j/b)*np.sin(np.pi*j/(b+1))/np.sqrt(b*(b+1))))
        self.beta_helix=np.cos(self.phi/2)**(2*b)
        self.marker_norm=self.head_overlap**2*self.beta_wrap*self.beta_helix
        self.marker_u=J0/self.marker_norm
        c=counts(q);self.counts=c
        self.alpha=np.asarray(c['port_counts'],float)/c['clock_vertices']
        self.beta=np.asarray(c['density_counts'],float)/c['clock_vertices']
        tg=target(self.m);self.theta=tg['theta'];self.mu=tg['mu']
        self.target=tg
        self.J=np.array([-.5*np.array([[np.sin(t),np.cos(t)],[-np.cos(t),np.sin(t)]])
                         for t in self.theta])
        self.payload_energy=float(np.sum(-2*np.sin(tg['width'])/np.pi+
                                         self.mu*(tg['width']/np.pi-.5)))
        # Minimizes the first-order physical energy density at nu=1/2.
        self.kappa=2*J0/np.pi-self.payload_energy

    def encode(self,s):return ((s.selector*self.s+s.wrapper)*2+s.A)*2+s.B
    def decode(self,z):
        z,B=divmod(int(z),2);z,A=divmod(z,2);selector,wrapper=divmod(z,self.s)
        assert 0<=selector<4
        return PhysicalSite(selector,wrapper,A,B)
    def occupied(self,tag,core,filt):return 3+(tag*self.Q+core)*2+filt
    def parts(self,w):
        if w<3:return None
        z,filt=divmod(w-3,2);tag,core=divmod(z,self.Q)
        return tag,core,filt
    def arrow(self,w):return self.wrapper_arrows[w if w<3 else 3+self.parts(w)[0]]
    def legal(self,a,b):return self.selector_arrows[a][1]==self.selector_arrows[b][0]
    def charge(self,s):
        p=self.parts(s.wrapper)
        return 0 if p is None else self.core.charge(self.core.decode(p[1]))

    def wrapper(self,x,y):
        out=defaultdict(complex);px,py=self.parts(x),self.parts(y)
        out[x,y]+=6*self.U*(self.arrow(x)[1]!=self.arrow(y)[0])
        out[x,y]+=self.U*((x<3)!=(y<3))
        # Each vacancy/particle move also changes the vacancy's arrow type.
        for tag,hl,hr in ((0,0,1),(1,1,1),(2,1,2)):
            if x==hl and py is not None and py[0]==tag:out[y,hr]-=self.w
            if y==hr and px is not None and px[0]==tag:out[hl,x]-=self.w
        if px is not None and py is not None and (px[0],py[0]) in ((0,1),(0,2),(1,1),(1,2)):
            tx,cx,fx=px;ty,cy,fy=py
            for (dx,dy),v in self.core.bulk(cx,cy).items():
                out[self.occupied(tx,dx,fx),self.occupied(ty,dy,fy)]+=v
            for idx in range(4):
                v=self.hf[idx,2*fx+fy]
                if v:out[self.occupied(tx,cx,idx//2),self.occupied(ty,cy,idx%2)]+=v
        for side,p in (('L',px),('R',py)):
            if p is None or p[0]!=(0 if side=='L' else 2):continue
            tag,c,f=p
            out[x,y]-=2*self.U
            for d,v in self.core.boundary(c,side).items():
                dest=self.occupied(tag,d,f)
                out[(dest,y) if side=='L' else (x,dest)]+=v
            mat=self.af if side=='L' else self.bf
            for g in range(2):
                dest=self.occupied(tag,c,g)
                out[(dest,y) if side=='L' else (x,dest)]+=mat[g,f]
        return {key:v for key,v in out.items() if v}

    def core_port(self,c,side,f,component):
        s=self.core.decode(c);out=defaultdict(complex)
        if s.role!=side or not s.ready or s.head:return {}
        a,b=divmod(2*f,self.q//2)
        for p in (0,1):
            odd=self.q//2+b+p
            if s.record==a:dest=odd;x=(-1)**p;y=1j*x
            elif s.record==odd:dest=a;x=(-1)**p;y=-1j*x
            else:continue
            v=(x,y)[component] if side=='L' else (y,-x)[component]
            out[self.core.encode(replace(s,record=dest))]+=v
        return dict(out)

    def port(self,s,side,f,component):
        p=self.parts(s.wrapper)
        if p is None or p[0]!=(0 if side=='L' else 2):return {}
        tag,c,fil=p;out=defaultdict(complex)
        mat=self.ML if side=='L' else self.MR
        bit=s.B if side=='L' else s.A
        for dest,v in self.core_port(c,side,f,component).items():
            for g in range(2):
                d=replace(s,wrapper=self.occupied(tag,dest,fil),**({'B':g} if side=='L' else {'A':g}))
                out[d]+=v*mat[g,bit]/self.alpha[0 if side=='L' else 1,f]
        return dict(out)

    def charge_port(self,s,side,f,raising):
        """The raw record creation/annihilation map, with compensators.

        Keeping charge definite avoids roundoff cancellation between XX, YY,
        XY and YX terms. Right-side dressing is supplied by its code map.
        """
        p=self.parts(s.wrapper)
        if p is None or p[0]!=(0 if side=='L' else 2):return {}
        tag,c,fil=p;core=self.core.decode(c)
        if core.role!=side or not core.ready or core.head:return {}
        row,col=divmod(2*f,self.q//2)
        mat=self.ML if side=='L' else self.MR
        bit=s.B if side=='L' else s.A;out=defaultdict(complex)
        for parity in (0,1):
            odd=self.q//2+col+parity
            old,new=(row,odd) if raising else (odd,row)
            if core.record!=old:continue
            dest=self.core.encode(replace(core,record=new))
            for g in range(2):
                d=replace(s,wrapper=self.occupied(tag,dest,fil),**({'B':g} if side=='L' else {'A':g}))
                out[d]+=(-1)**parity*mat[g,bit]/self.alpha[0 if side=='L' else 1,f]
        return dict(out)

    def density_value(self,x,y,f):
        px,py=self.parts(x),self.parts(y)
        if px is None or py is None or (px[0],py[0])!=(0,1):return 0
        a,b=self.core.decode(px[1]),self.core.decode(py[1])
        if a.role!='L' or b.role!='A' or not a.ready or a.head or b.head:return 0
        row,col=divmod(2*f,self.q//2)
        value=0
        for p in (0,1):
            if b.record!=(row+col+p)%(self.q//2):continue
            if a.record==row:value+=1
            elif a.record==self.q//2+col+p:value-=1
        return value

    def column(self,i,j,base=1.,perturbation=0.):
        a,b=self.decode(i),self.decode(j);out=defaultdict(complex)
        legal=self.legal(a.selector,b.selector)
        internal=legal and b.selector!=0
        cut=legal and b.selector==0
        def put(x,y,v):
            if v:out[self.encode(x),self.encode(y)]+=v
        if base:
            put(a,b,base*(self.Ksel*(not legal)+self.marker_mu/2*((a.selector==0)+(b.selector==0))))
            if (a.selector,b.selector) in ((1,2),(2,3)):
                u,v=(2,3) if a.selector==1 else (1,2)
                put(replace(a,selector=u),replace(b,selector=v),-base*self.t)
            if internal:
                for (x,y),v in self.wrapper(a.wrapper,b.wrapper).items():
                    put(replace(a,wrapper=x),replace(b,wrapper=y),base*v)
                for track in ('A','B'):
                    old=2*getattr(a,track)+getattr(b,track)
                    for new in range(4):
                        put(replace(a,**{track:new//2}),replace(b,**{track:new%2}),base*self.he[new,old])
            if cut:
                for z in range(2):
                    put(replace(a,B=z),b,base*self.pin[z,a.B])
                    put(a,replace(b,A=z),base*self.pin[z,b.A])
        if perturbation:
            eps=perturbation
            put(a,b,eps*self.kappa/2*((a.selector==0)+(b.selector==0)))
            # U_b on both equalizer bits of the transferred (left) site.
            phase=np.exp(1j*self.b*self.phi*(a.A+a.B-1))
            if (a.selector,b.selector)==(0,1) and a.wrapper==0:
                put(replace(a,selector=3,wrapper=2),replace(b,selector=0),-eps*self.marker_u*phase)
            if (a.selector,b.selector)==(3,0) and a.wrapper==2:
                put(replace(a,selector=0,wrapper=0),replace(b,selector=1),-eps*self.marker_u*phase.conjugate())
            if cut:
                for f in range(self.m):
                    for raising in (True,False):
                        hopping=-eps*np.exp((1j if raising else -1j)*self.theta[f])
                        for aa,va in self.charge_port(a,'R',f,raising).items():
                            for bb,vb in self.charge_port(b,'L',f,not raising).items():
                                put(aa,bb,hopping*va*vb)
            if a.selector==0 and internal:
                value=sum(-self.mu[f]/2*self.density_value(a.wrapper,b.wrapper,f)/self.beta[0,f]
                          for f in range(self.m))
                for z in range(2):put(replace(a,B=z),b,eps*value*self.MD[z,a.B])
        return {key:v for key,v in out.items() if v}

    def report(self):
        return dict(q=self.q,m=self.m,core_length=self.k,cell_lengths=[self.b,self.b+1],
                    core_dimension=self.Q,site_dimension=self.d,
                    selector_t=self.t,selector_mu=self.marker_mu,selector_penalty=self.Ksel,
                    base_energy_density=self.base_density,filter_scale=self.filter_scale,
                    vacancy_hop=self.w,endpoint_occupations=self.r,
                    two_site_occupations=self.s2,marker_visibility=self.marker_norm,
                    marker_coefficient=self.marker_u,marker_field=self.kappa,
                    payload_energy_per_cell=self.payload_energy)
