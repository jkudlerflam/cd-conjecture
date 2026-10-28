#!/usr/bin/env python3
"""Literal lambda=0 quantum selector; independent finite-ring construction.

Only writes beside this file. No existing numerical implementation imported
except an explicit, read-only 25-entry bond comparison in validate().
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import importlib.util
import itertools as it
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
# occupied alphabet L0,R0,Xd_v0,M_empty,M_v; vacancy=-1.
LFT, RGT, TOK, ME, MV = range(5)

def parameters(b):
    a=lambda l: 2*(1-math.cos(math.pi/(l+1)))
    delta=lambda l:(l+1)*a(l)-l*a(l+1)
    mu=2-(delta(b)+delta(b-1))/2
    return mu, (mu-2*math.cos(math.pi/(b+1)))/b

def frozen(pattern,mu):
    return sum(mu-2*math.cos(math.pi/(l+1)) for l in pattern)

def compositions(n,k):
    for cuts in it.combinations(range(n+k-1),k-1):
        ends=(-1,)+cuts+(n+k-1,)
        yield tuple(ends[j+1]-ends[j]-1 for j in range(k))

def charge(N,holes):
    gaps=list(compositions(holes,N)); ix={g:i for i,g in enumerate(gaps)}
    rows=[];cols=[]
    for col,g in enumerate(gaps):
        for j,v in enumerate(g):
            if v:
                for s in (-1,1):
                    out=list(g);out[j]-=1;out[(j+s)%N]+=1
                    rows.append(ix[tuple(out)]);cols.append(col)
    H=sp.coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(len(gaps),)*2).tocsr()
    length=N+holes
    E=-2*math.sin(math.pi*N/length)/math.sin(math.pi/length)
    alpha=N/length-math.sin(math.pi*N/length)**2/(length*N*math.sin(math.pi/length)**2)
    return gaps,ix,H,E,alpha

def word(pattern,walkers):
    out=[]
    for length,r in zip(pattern,walkers):
        if r==length-1:out += [LFT]*(length-1)+[MV]
        else:out += [LFT]*r+[TOK]+[RGT]*(length-2-r)+[ME]
    return tuple(out)

def primitive_cells(pattern):
    k=len(pattern)
    return next(p for p in range(1,k+1) if k%p==0 and pattern==pattern[:p]*(k//p))

def orbit_projection(pattern,gaps,gix,walkers,wix):
    """Physical momentum 0: quotient repeated material roots, normalized orbits.

    A marker-rooted covering is repeated r times iff the cell-length word
    repeats r times. Nonperiodic words have one distinguished root. Physical
    translation has already been projected to momentum 0; its phases are 1.
    """
    W=len(walkers);d=len(gaps)*W;p=primitive_cells(pattern);r=len(pattern)//p
    if r==1:return sp.identity(d,format='csr'),r
    shift=sum(pattern[:p])
    gp=np.array([gix[g[shift:]+g[:shift]] for g in gaps])
    wp=np.array([wix[w[p:]+w[:p]] for w in walkers])
    perm=(gp[:,None]*W+wp[None,:]).ravel()
    canon=np.arange(d,dtype=np.int64);cur=canon.copy()
    for _ in range(1,r):
        cur=perm[cur];canon=np.minimum(canon,cur)
    _,labels,counts=np.unique(canon,return_inverse=True,return_counts=True)
    U=sp.coo_matrix((1/np.sqrt(counts[labels]),(np.arange(d),labels)),shape=(d,len(counts))).tocsr()
    assert np.max(abs((U.T@U-sp.identity(U.shape[1])).data),initial=0)<1e-12
    return U,r

def model(pattern,holes,mu):
    N=sum(pattern);gaps,gix,Hg,Ec,alpha=charge(N,holes)
    walkers=list(it.product(*(range(l) for l in pattern)));wix={w:i for i,w in enumerate(walkers)}
    W=len(walkers);G=len(gaps);gaparray=np.array(gaps,dtype=np.int16);A=(gaparray==0)
    diag=np.zeros(G);starts=np.cumsum([0]+list(pattern[:-1]))
    Hc=sp.kron(Hg,sp.identity(W,format='csr'),format='csr')
    V=sp.csr_matrix(Hc.shape)
    for ci,(start,length) in enumerate(zip(starts,pattern)):
        marker=start+length-1
        # Literal h23 distributes mu/2 to EACH occupied-pair endpoint.
        diag += mu/2*(A[:,(marker-1)%N].astype(float)+A[:,marker].astype(float))
        for pos in range(length-1):
            rows=[];cols=[]
            for col,w in enumerate(walkers):
                if w[ci]==pos:
                    v=list(w);v[ci]+=1;row=wix[tuple(v)]
                    rows.extend((row,col));cols.extend((col,row))
            F=sp.coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(W,W)).tocsr()
            V += sp.kron(sp.diags(A[:,start+pos].astype(float)),F,format='csr')
    V += sp.diags(np.repeat(diag,W))
    U,repeats=orbit_projection(tuple(pattern),gaps,gix,walkers,wix)
    if repeats>1:Hc=(U.T@Hc@U).tocsr();V=(U.T@V@U).tocsr()
    asym=max(np.max(abs((X-X.T).data),initial=0) for X in (Hc,V))
    assert asym<1e-12
    # Exact cell sine vectors; charge PF vector is calculated in the full bath.
    ec,vc=eigsh(Hg,k=1,which='SA',tol=1e-12,v0=np.ones(G))
    assert abs(ec[0]-Ec)<1e-9
    vwalk=np.ones(W)
    for ci,length in enumerate(pattern):
        vwalk *= np.sqrt(2/(length+1))*np.sin(math.pi*(np.array(walkers)[:,ci]+1)/(length+1))
    v0=np.asarray(U.T@np.kron(abs(vc[:,0]),vwalk)).ravel();v0/=np.linalg.norm(v0)
    first=float(v0@(V@v0))
    assert abs(first-alpha*frozen(pattern,mu))<2e-9
    return Hc,V,v0,dict(N=N,L=N+holes,pattern=list(pattern),cover_dimension=G*W,dimension=Hc.shape[0],root_repeats=repeats,charge_energy=Ec,alpha=alpha,frozen_energy=frozen(pattern,mu),trial_first_order=first,hermiticity_error=asym)

def solve_pattern(pattern,holes,mu,Js):
    tic=time.monotonic();Hc,V,v0,rec=model(pattern,holes,mu);rec['build_seconds']=time.monotonic()-tic
    rows=[]
    for J in sorted(Js,reverse=True):
        tic=time.monotonic();H=Hc+J*V
        es,vs=eigsh(H,k=1,which='SA',v0=v0,tol=2e-12,maxiter=10000,ncv=30)
        v0=vs[:,0];energy=float(es[0]);res=float(np.linalg.norm(H@v0-energy*v0))
        row=dict(J=J,energy=energy,residual=res,correction=energy-rec['charge_energy']-J*rec['alpha']*rec['frozen_energy'],seconds=time.monotonic()-tic)
        row['correction_over_J2']=row['correction']/J**2;rows.append(row)
        print(json.dumps(dict(pattern=list(pattern),L=rec['L'],dimension=rec['dimension'],**row)),flush=True)
    rec['sweep']=rows
    return rec

def reduced_bond(mu,K=10):
    h=np.zeros((25,25));markers={ME,MV};arrows=((0,0),(1,1),(0,1),(1,0),(0,0))
    for a,b in it.product(range(5),repeat=2):
        h[5*a+b,5*a+b]=mu/2*((a in markers)+(b in markers))+K*(arrows[a][1]!=arrows[b][0])
    for (a,b),(c,d) in (((TOK,RGT),(LFT,TOK)),((TOK,ME),(LFT,MV))):
        h[5*a+b,5*c+d]=h[5*c+d,5*a+b]=-1
    return h

def direct_physical(pattern,holes,mu,J):
    N=sum(pattern);length=N+holes
    strings=set()
    for w in it.product(*(range(l) for l in pattern)):
        s=word(pattern,w)
        strings.update(s[k:]+s[:k] for k in range(N))
    states=[]
    for positions in it.combinations(range(length),N):
        for s in sorted(strings):
            out=[-1]*length
            for x,a in zip(positions,s):out[x]=a
            states.append(tuple(out))
    ix={s:i for i,s in enumerate(states)};h=reduced_bond(mu)
    moves={i:[(j,h[j,i]) for j in np.nonzero(h[:,i])[0]] for i in range(25)}
    rows=[];cols=[];vals=[]
    for col,s in enumerate(states):
        for x in range(length):
            y=(x+1)%length;a,b=s[x],s[y]
            if (a<0)!=(b<0):
                out=list(s);out[x],out[y]=b,a
                rows.append(ix[tuple(out)]);cols.append(col);vals.append(-1)
            elif a>=0 and b>=0:
                for j,amp in moves[5*a+b]:
                    out=list(s);out[x],out[y]=divmod(j,5)
                    rows.append(ix[tuple(out)]);cols.append(col);vals.append(J*amp)
    H=sp.coo_matrix((vals,(rows,cols)),shape=(len(states),)*2).tocsr()
    e,v=eigsh(H,k=1,which='SA',tol=1e-12,v0=np.ones(len(states)))
    return float(e[0]),len(states),float(np.linalg.norm(H@v[:,0]-e[0]*v[:,0]))

def validate():
    # Read-only comparison to literal source, not used to build matrices above.
    path=ROOT/'contrib/explicit_candidate/bond23.py'
    spec=importlib.util.spec_from_file_location('literal_bond23',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    mu=parameters(3)[0];h=m.bond(m.params(lam=0,musel=mu))
    inds=[m.IX[x] for x in [('L',None,0),('R',None,0),('Xd','v',0),('Me',None,None),('M','v',None)]]
    ids=np.array([23*a+b for a,b in it.product(inds,repeat=2)])
    others=np.setdiff1d(np.arange(529),ids)
    diff=float(np.max(abs(h[np.ix_(ids,ids)]-reduced_bond(mu))))
    leak=float(np.max(abs(h[np.ix_(others,ids)])))
    assert max(diff,leak)<1e-12
    checks=[]
    for pattern,holes in [((2,2),4),((1,3),4),((4,),4),((1,1,2),4),((3,3),3)]:
        Hc,V,v0,rec=model(pattern,holes,mu);J=.073
        eq,v=eigsh(Hc+J*V,k=1,which='SA',v0=v0,tol=1e-12)
        ed,dim,res=direct_physical(pattern,holes,mu,J)
        check=dict(pattern=list(pattern),L=sum(pattern)+holes,physical_dimension=dim,quotient_dimension=Hc.shape[0],quotient_energy=float(eq[0]),direct_energy=ed,error=float(eq[0]-ed),direct_residual=res)
        assert abs(check['error'])<1e-10;checks.append(check)
    # Vacuum transition form factors at b=3: sum bonds has no excited/ground entry.
    l=3;cell=np.diag(-np.ones(l-1),1)+np.diag(-np.ones(l-1),-1)+mu*np.eye(l)
    e,u=np.linalg.eigh(cell);terms=[]
    for j in range(l):
        f=np.zeros((l,l))
        if j<l-1:f[j,j+1]=f[j+1,j]=-1
        if j==l-2 or j==l-1:f+=mu/2*np.eye(l)
        terms.append(f)
    ff=np.array([[u[:,a]@f@u[:,0] for f in terms] for a in range(1,l)])
    return dict(literal_bond_max_error=diff,literal_bond_leakage=leak,ring_checks=checks,cell_offdiagonal_formfactor_sum_max=float(np.max(abs(ff.sum(axis=1)))))

def necklaces(total):
    allpatterns=set()
    for k in range(1,total+1):
        for g in compositions(total-k,k):
            p=tuple(x+1 for x in g)
            canonical=min(p[j:]+p[:j] for j in range(k))
            allpatterns.add(canonical)
    return sorted(allpatterns,key=lambda p:(len(p),p))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--suite',choices=['validation','small','large','extended','giant'],default='small');ap.add_argument('--output');args=ap.parse_args()
    tic=time.monotonic();out=dict(status='numerical diagnostic, not an n=23 proof',construction='independent contact-weighted lambda=0 vacuum selector, physical momentum-zero root quotient',validation=validate())
    jobs=[]
    if args.suite=='small':
        for b,N in [(2,4),(2,6),(3,6)]:
            selected=(b,)*(N//b)
            for p in necklaces(N):jobs.append((b,p,N))
    elif args.suite=='large':
        for b,N,patterns in [(2,8,[(2,2,2,2),(1,2,2,3),(1,2,3,2),(4,2,2),(8,)]),(3,9,[(3,3,3),(2,3,4),(6,3),(9,)])]:
            for p in patterns:jobs.append((b,p,N))
    elif args.suite=='giant':
        for b,N in [(2,4),(2,6),(2,8),(3,6),(3,9)]:
            for p in [(b,)*(N//b),(N,)]:jobs.append((b,p,N))
    elif args.suite=='extended':
        for p in [(2,2,2,2,2),(1,3,2,2,2),(1,2,3,2,2),(10,)]:jobs.append((2,p,10))
    Js=[.1,.05,.02,.01,.005,.002]
    results=[]
    for b,p,holes in jobs:
        mu,eb=parameters(b);rec=solve_pattern(tuple(p),holes,mu,Js);rec.update(selected_b=b,mu_selector=mu,selected_frozen_density=eb);results.append(rec)
        out.update(results=results,seconds=time.monotonic()-tic)
        (HERE/(args.output or f'{args.suite}_results.json')).write_text(json.dumps(out,indent=2)+'\n')
    comparisons=[]
    for r in results:
        ref=next(x for x in results if x['selected_b']==r['selected_b'] and x['N']==r['N'] and x['L']==r['L'] and x['pattern']==[r['selected_b']]*(r['N']//r['selected_b']))
        k=sum(x!=r['selected_b'] for x in r['pattern']);bare=r['frozen_energy']-ref['frozen_energy']
        for x,y in zip(r['sweep'],ref['sweep']):
            delta=x['energy']-y['energy'];remainder=delta-x['J']*r['alpha']*bare
            comparisons.append(dict(b=r['selected_b'],N=r['N'],L=r['L'],pattern=r['pattern'],wrong_cells=k,J=x['J'],bare_frozen_excess=bare,exact_excess=delta,remainder=remainder,remainder_over_J2=remainder/x['J']**2,remainder_over_Jk=remainder/(x['J']*k) if k else None))
    out.update(results=results,comparisons=comparisons,seconds=time.monotonic()-tic)
    (HERE/(args.output or f'{args.suite}_results.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(suite=args.suite,seconds=out['seconds'],comparisons=len(comparisons),validation=out['validation'])),flush=True)
if __name__=='__main__':main()
