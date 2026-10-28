"""Independent C18 checks of the two-port adaptation of C8's sweep clock.
No files outside this directory are imported or written.
"""
import itertools
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal

BULK = ('L','0','1','C0','C1','R0','R1','I','Z',
        'A0','A1','D0','D1','E0','E1')
HEAD = set(('C0','C1','R0','R1','I','D0','D1','E0','E1'))
LEFT, RIGHT = 'PL', 'PR'
ALL = BULK + (LEFT, RIGHT)
PAIRS = {('L','L'), ('I','Z'), ('Z','Z'), (RIGHT,LEFT)}
OLDHEAD = ('C0','C1','R0','R1','I')
for h in OLDHEAD: PAIRS.add(('L',h))
for h in OLDHEAD[:-1]:
    for x in '01': PAIRS.add((h,x))
for x,y in itertools.product('01',repeat=2):
    PAIRS.add((x,y))
    for h in 'ADE': PAIRS.add(('A'+x,h+y))
    for h in 'DE': PAIRS.add((h+x,y))
START = ('L',)+OLDHEAD+('A0','A1','D0','D1','E0','E1')
END = ('0','1','Z')+OLDHEAD+('D0','D1','E0','E1')
PAIRS.update((LEFT,x) for x in START)
PAIRS.update((x,RIGHT) for x in END)
RULES=[]
def rule(a,b): RULES.append((tuple(a.split()),tuple(b.split())))
for x in '01':
    rule('C1 '+x,'L C'+x)
    rule('L R'+x,'R0 '+x)
    rule('PL R'+x,'PL D'+x)
    rule('D'+x+' PR','E'+x+' PR')
    rule('PL E'+x,'PL C'+x)
    for y in '01':
        rule('D'+x+' '+y,'A'+x+' D'+y)
        rule('A'+x+' E'+y,'E'+x+' '+y)
rule('I Z','L I')
rule('I PR','R0 PR')
rule('C0','R1')


def defects(word):
    return tuple(int(p not in PAIRS) for p in zip(word[:-1],word[1:]))


def neighbors(word, periodic=False):
    N=len(word)
    for a,b in RULES:
        width=len(a)
        for i in range(N if periodic else N-width+1):
            inds=tuple((i+j)%N for j in range(width))
            if len(set(inds)) != width: continue
            old=tuple(word[j] for j in inds)
            target=b if old==a else a if old==b else None
            if target is not None:
                nxt=list(word)
                for j,x in zip(inds,target): nxt[j]=x
                yield tuple(nxt)


def C(n): return (2*n+4)*2**n+n-4

def a(n): return 4*np.sin(np.pi/(2*(C(n)+1)))**2

def selection(n):
    b=n+2
    dn=(b+1)*a(n)-b*a(n+1)
    dp=b*a(n-1)-(b-1)*a(n)
    delta=(dn+dp)/2
    gamma=(a(n-1)-2*a(n)+a(n+1))/2
    energy=(a(n)-delta)/b
    gap=4*np.sin(3*np.pi/(2*(C(n)+1)))*np.sin(np.pi/(2*(C(n)+1)))
    return delta,gamma,energy,gap


def local_defect_checks():
    # Equality of defect patterns for all contexts of each rewrite.
    count=0
    for aa,bb in RULES:
        for l,r in itertools.product(ALL,repeat=2):
            assert defects((l,)+aa+(r,))==defects((l,)+bb+(r,))
            count+=1
    return count


def open_graph(n):
    words=list(itertools.product(BULK,repeat=n))
    index={w:i for i,w in enumerate(words)}
    adj=[[] for _ in words]
    ds=[]
    for i,w in enumerate(words):
        full=(LEFT,)+w+(RIGHT,)
        sig=defects(full); ds.append(sum(sig))
        for v in neighbors(full):
            assert v[0]==LEFT and v[-1]==RIGHT
            assert defects(v)==sig
            j=index[v[1:-1]]; adj[i].append(j)
        assert len(set(adj[i]))==len(adj[i])
        h=sum(x in HEAD for x in w)
        assert len(adj[i])<=2*h
        assert h<=sum(sig)+1
    seen=set(); sizes=[]; legal=[]; illegal_min=None; cycles=0
    for i in range(len(words)):
        if i in seen: continue
        stack=[i]; seen.add(i); comp=[]
        while stack:
            j=stack.pop(); comp.append(j)
            for k in adj[j]:
                if k not in seen: seen.add(k); stack.append(k)
        sizes.append(len(comp))
        edgecount=sum(len(adj[j]) for j in comp)//2
        if edgecount>=len(comp): cycles+=1
        if ds[i]==0:
            assert edgecount==len(comp)-1
            assert max(len(adj[j]) for j in comp)<=2
            legal.append(len(comp))
        elif n<=3:
            ids={j:k for k,j in enumerate(comp)}
            H=np.eye(len(comp))*8*ds[i]
            for j in comp:
                for k in adj[j]: H[ids[j],ids[k]]-=1
            ev=float(np.linalg.eigvalsh(H)[0])
            assert ev>=-2+(8-2)*ds[i]-1e-10
            illegal_min=ev if illegal_min is None else min(illegal_min,ev)
    assert legal==[C(n)]
    return dict(n=n,words=len(words),components=len(sizes),legal_vertices=legal[0],
                illegal_min_energy=illegal_min,cyclic_components=cycles)


def ring_checks(L):
    legal=0; edges=0; maxheads=0
    for w in itertools.product(ALL,repeat=L):
        sig=tuple(int((w[i],w[(i+1)%L]) not in PAIRS) for i in range(L))
        h=sum(x in HEAD for x in w); v=sum(sig)
        nb=list(neighbors(w,True)); edges+=len(nb)
        assert len(nb)<=2*h
        # After cuts, each legal fragment has <= markers+1 heads.
        assert h<=w.count(LEFT)+v
        for z in nb:
            assert tuple(int((z[i],z[(i+1)%L]) not in PAIRS) for i in range(L))==sig
        if v==0:
            legal+=1
            assert w.count(LEFT)==w.count(RIGHT)
            if LEFT not in w: assert h==0
        maxheads=max(maxheads,h)
    return dict(L=L,words=len(ALL)**L,legal_words=legal,directed_moves=edges)


def numeric_selection():
    results=[]
    for n in range(2,41):
        b=n+2; delta,gamma,en,gap=selection(n)
        assert 0<delta<2
        vals=[((a(k)-delta)/(k+2)-en,k) for k in range(1,101) if k!=n]
        best=min(vals)
        assert best[1]==n+1
        assert abs(best[0]/(gamma/(b+1))-1)<1e-10
        excess=[(k+2)*v for v,k in vals]
        assert abs(min(excess)/gamma-1)<1e-10
        if n in (2,3,4,6,8,12,20,40):
            results.append(dict(n=n,b=b,C=C(n),delta=delta,gamma=gamma,
                                density_margin=gamma/(b+1),clock_gap=gap,
                                shifted_energy=a(n)))
    return results


def endpoint_pins():
    ans=[]
    for n in range(1,9):
        c=C(n); d=np.full(c,2.); d[0]+=1
        ev=eigh_tridiagonal(d,-np.ones(c-1),select='i',select_range=(0,0),
                           lapack_driver='stebz',tol=1e-15)[0][0]
        ref=4*np.sin(np.pi/(2*(c+.5)))**2
        shift=ev-a(n)
        analytic=ref-a(n)
        assert abs(ev-ref)<3e-15
        # A first-slice pin expectation is an upper bound on the lift.
        weight=2*np.sin(np.pi/(c+1))**2/(c+1)
        assert 0<shift<=weight+3e-15
        gap=4*np.sin(3*np.pi/(2*(c+1)))*np.sin(np.pi/(2*(c+1)))
        bound=gap*weight/(gap+1)
        assert shift>=bound-3e-15
        ans.append(dict(n=n,C=c,lift=shift,exact_lift=analytic,
                        proved_lower_bound=bound,variational_upper_bound=weight,
                        C3_lift=shift*c**3))
    return ans


def plateau_counterexample():
    # A monotone but nonconvex clock length skips lengths under selection.
    cc=np.array([9,9,29],float)
    aa=4*np.sin(np.pi/(2*(cc+1)))**2
    # b=3,4,5. Middle selects iff d_next<d_prev, which fails.
    dprev=4*aa[0]-3*aa[1]; dnext=5*aa[1]-4*aa[2]
    assert dnext>dprev
    return dict(vertices=cc.tolist(),left_crossing=float(dprev),right_crossing=float(dnext))

if __name__=='__main__':
    start=time.monotonic()
    out=dict(local_contexts=local_defect_checks(),
             open=[open_graph(n) for n in range(1,5)],
             rings=[ring_checks(L) for L in (3,4)],
             selection=numeric_selection(),endpoint_pins=endpoint_pins(),
             monotonicity_counterexample=plateau_counterexample())
    out['seconds']=time.monotonic()-start
    Path('clock_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
