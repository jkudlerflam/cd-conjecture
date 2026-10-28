"""Adversarial finite checks of the full periodic sparse bond oracle."""
import copy
import json
from dataclasses import replace
from itertools import product
import numpy as np
from periodic_clean_bond import PeriodicBond,PhysicalSite
from clean_core_bond import Site
from clean_readout_checks import clock_stream
from arrow_wrapper_checks import local_bond


def small_wrapper_comparison(H):
    # A different two-state core tests every wrapper column, including illegal
    # arrows and oriented endpoint pins, against an independent dense builder.
    h=np.zeros((4,4),complex);h[1,1]=h[2,2]=.7;h[1,2]=h[2,1]=-.7
    a=np.diag([0.,.31]);b=np.diag([0.,.47])
    class Tiny:
        d=2
        def bulk(self,i,j):return {divmod(z,2):h[z,2*i+j] for z in range(4) if h[z,2*i+j]}
        def boundary(self,i,side):
            mat=a if side=='L' else b
            return {z:mat[z,i] for z in range(2) if mat[z,i]}
    T=copy.copy(H);T.core=Tiny();T.Q=2;T.s=15;T.d=240
    hp=np.kron(h,np.eye(4))+np.kron(np.eye(4),T.hf)
    # Reorder (core1,core2,filter1,filter2) into occupied physical sites.
    perm=[8*x+4*y+2*f+g for x in range(2) for f in range(2) for y in range(2) for g in range(2)]
    hp=hp[np.ix_(perm,perm)]
    ap=np.kron(a,np.eye(2))+np.kron(np.eye(2),T.af)
    bp=np.kron(b,np.eye(2))+np.kron(np.eye(2),T.bf)
    dense=local_bond(4,hp,ap,bp,T.U,T.w,2*T.U,6*T.U)[0]
    measured=np.zeros_like(dense)
    for x,y in product(range(15),repeat=2):
        for (u,v),z in T.wrapper(x,y).items():measured[15*u+v,15*x+y]=z
    err=float(np.max(abs(measured-dense)))
    assert err<1e-11
    return dict(core_dimension=2,wrapper_dimension=15,columns=225,max_error=err)


def actual_columns(H):
    rng=np.random.default_rng(9028);pairs=[]
    for _ in range(1600):pairs.append(tuple(map(int,rng.integers(H.d,size=2))))
    # Drive each cut channel and every color input, both marker orientations,
    # arbitrary clock letters, and both equalizer components explicitly.
    def occupied(role,record,tag,sel,A,B,ready=1):
        c=H.core.encode(Site('D0',role,0,0,ready,record,0))
        return PhysicalSite(sel,H.occupied(tag,c,0),A,B)
    for x,y,A,B,C,D in product(range(H.q+1),range(H.q+1),range(2),range(2),range(2),range(2)):
        for sel in (2,3):
            pairs.append((H.encode(occupied('R',x,2,sel,A,B)),H.encode(occupied('L',y,0,0,C,D))))
        for sel in (1,2):
            pairs.append((H.encode(occupied('L',x,0,0,A,B)),H.encode(occupied('A',y,1,sel,C,D,0))))
    for sel,wrap,other in ((0,0,1),(3,2,0)):
        for A,B,C,D in product(range(2),repeat=4):
            pairs.append((H.encode(PhysicalSite(sel,wrap,A,B)),H.encode(PhysicalSite(other,1,C,D))))
    errors=[0.,0.];entries=[0,0];changed=0
    for part,(base,eps) in enumerate(((1.,0.),(0.,1.))):
        cache={}
        def column(i,j):
            key=(i,j)
            if key not in cache:cache[key]=H.column(i,j,base=base,perturbation=eps)
            return cache[key]
        for i,j in pairs:
            a,b=H.decode(i),H.decode(j)
            assert (H.encode(a),H.encode(b))==(i,j)
            for (u,v),value in column(i,j).items():
                aa,bb=H.decode(u),H.decode(v)
                rev=column(u,v).get((i,j),0).conjugate()
                errors[part]=max(errors[part],abs(value-rev)/(1+abs(value)))
                assert H.charge(a)+H.charge(b)==H.charge(aa)+H.charge(bb),(part,a,b,aa,bb,value)
                assert (a.selector==0)+(b.selector==0)==(aa.selector==0)+(bb.selector==0)
                if (a.selector,b.selector)!=(aa.selector,bb.selector):
                    changed+=1
                    for left,right in product(range(4),repeat=2):
                        old=(left,a.selector,b.selector,right);new=(left,aa.selector,bb.selector,right)
                        assert [H.legal(*old[k:k+2]) for k in range(3)]==[H.legal(*new[k:k+2]) for k in range(3)]
                entries[part]+=1
        assert errors[part]<3e-13
    return dict(input_columns_per_part=len(pairs),entries=entries,
                max_relative_hermiticity_errors=errors,selector_changing_entries=changed,
                charge_marker_and_legality_failures=0)


def history_compression(H):
    m=H.m;D=2**m;assert H.q==4
    bits=list(product(range(2),repeat=m));cs=[H.core.M.initialize(H.q,x) for x in bits]
    ports=np.zeros((2,m,2,D,D),complex);density=np.zeros((m,D,D),complex)
    def core_site(c,j,clock):
        head=H.core.head_encode(c.control,c.p,c.x,c.phase) if c.pos==j else 0
        return H.core.encode(Site(clock,c.roles[j],*c.tape[j],head))
    C=0
    for clocks,lab,j in clock_stream(H.k):
        for c in cs:
            if lab=='out':H.core.M.bond(c,j,H.q)
            elif lab=='back':H.core.M.onsite(c,j+1,H.q)
            elif lab=='finish':H.core.M.onsite(c,0,H.q)
        keys={c.snapshot():i for i,c in enumerate(cs)}
        for x,c in enumerate(cs):
            for side,pos in (('L',0),('R',H.k-1)):
                ci=core_site(c,pos,clocks[pos]);old=c.tape[pos][3]
                for f,u in product(range(m),range(2)):
                    for dest,v in H.core_port(ci,side,f,u).items():
                        c.tape[pos][3]=H.core.decode(dest).record
                        y=keys.get(c.snapshot())
                        if y is not None:ports[0 if side=='L' else 1,f,u,y,x]+=v
                        c.tape[pos][3]=old
            ci,cj=core_site(c,0,clocks[0]),core_site(c,1,clocks[1])
            for f in range(m):density[f,x,x]+=H.density_value(H.occupied(0,ci,0),H.occupied(1,cj,0),f)
        C+=1
    for side,f in product(range(2),range(m)):ports[side,f]/=C*H.alpha[side,f]
    for f in range(m):density[f]/=C*H.beta[0,f]
    # The physical helix compensators are contracted separately in their
    # actual selected states. This avoids constructing 4^(2b+1) vectors.
    interpolation=0.
    for p,state in enumerate(H.states):
        interpolation=max(interpolation,abs(H.r[p]*np.vdot(state,H.ML@state)-1),
                          abs(H.r[p]*np.vdot(state.conj(),H.MR@state.conj())-1),
                          abs(H.s2[p]*np.vdot(state,H.MD@state)-1))
    actual=np.zeros((D*D,D*D),complex)
    for f,u,v in product(range(m),range(2),range(2)):
        actual+=H.J[f,u,v]*np.kron(ports[1,f,u],ports[0,f,v])
    field=-sum(H.mu[f]/2*density[f] for f in range(m))
    actual+=np.kron(field,np.eye(D))+np.kron(np.eye(D),field)
    # Independent Fock-space construction with one common fermion ordering.
    fermions=[]
    for j in range(2*m):
        c=np.zeros_like(actual)
        for x in range(D*D):
            occupation=[(x>>(2*m-1-i))&1 for i in range(2*m)]
            if occupation[j]:c[x^(1<<(2*m-1-j)),x]=(-1)**sum(occupation[:j])
        fermions.append(c)
    reference=np.zeros_like(actual)
    for f in range(m):
        term=np.exp(1j*H.theta[f])*fermions[f].conj().T@fermions[m+f]
        reference-=term+term.conj().T
        for j in (f,m+f):reference+=H.mu[f]*(fermions[j].conj().T@fermions[j]-.5*np.eye(D*D))
    error=float(np.max(abs(actual-reference)))
    assert error<1e-12 and interpolation<1e-10
    # Direct minimization of the projected thermodynamic energy checks kappa.
    nu=np.linspace(.05,.95,901)
    e=(H.kappa+H.payload_energy-2*H.J0*np.sin(np.pi*nu)/np.pi)/(H.b+nu)
    optimum=float(nu[np.argmin(e)])
    assert abs(optimum-.5)<1e-12
    return dict(clock_vertices=C,logical_dimension=D,two_cell_dimension=D*D,
                max_two_cell_Fock_matrix_error=error,
                max_compensation_error=float(interpolation),
                selected_reference_length_density=optimum)


if __name__=='__main__':
    H=PeriodicBond()
    out=dict(parameters=H.report(),wrapper=small_wrapper_comparison(H),
             columns=actual_columns(H),compression=history_compression(H),
             scope='Finite local matrix and projection checks, no CFT estimate.')
    with open('periodic_bond_checks.json','w') as f:
        json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))
