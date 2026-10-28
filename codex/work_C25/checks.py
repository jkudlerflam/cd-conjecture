"""Independent finite checks of the C25 invariant reductions."""
import sys
sys.dont_write_bytecode=True
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,solve
import composition as c

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]])
Y=np.array([[0.,-1j],[1j,0.]])
Z=np.diag([1.,-1.]); plus=np.ones(2)/np.sqrt(2)
def tensor(items):
    a=np.array([[1.]])
    for b in items: a=np.kron(a,b)
    return a
def at(op,j,N): return tensor([op if k==j else I for k in range(N)])

def main():
    out={}
    # Raw full accessible charge-clock matrices, independent of block solver.
    ch=c.charge(2,.16); s=.1*ch['gap']; maxerr=0.; coderes=0.
    for side in ('L','R'):
        H0=np.kron(ch['H'],np.eye(c.C))
        for j in range(3): H0+=s*np.kron(np.diag(ch['D'][:,j]),c.PARTS[j])
        gs=np.kron(ch['chi'],c.PHI)
        coderes=max(coderes,np.linalg.norm(H0@gs))
        for i in (0,1):
            H=c.core_matrix(ch,side,i,s)
            src=np.kron(ch['occ'][side]*ch['chi'],c.PHI)
            e,a=c.measure(ch,side,i,s,(np.array([0.]),np.array([1.]),0.))
            for z in (0.,1.,3.):
                direct=s*np.dot(src,solve(H+s*z*np.eye(len(H)),src,assume_a='pos'))
                spectral=s*np.sum(a/(e+s*z))
                maxerr=max(maxerr,abs(direct-spectral))
    assert coderes<1e-13 and maxerr<2e-7
    out['raw_source_inverse_check']=dict(dimension=450,absolute_scaled_error=maxerr,
                                        code_residual=coderes,s=s)
    # Full 2^N helix parent, without using magnon-sector reduction.
    charges={p:c.charge(p,.16) for p in (1,2)}
    M=c.compensator(charges[1]['r'],charges[2]['r'],.2)
    hx=[]
    singlet=np.array([0,1,-1,0])/np.sqrt(2)
    hh=np.outer(singlet,singlet)
    pin=I-np.outer(plus,plus)
    for N in (5,6):
        H=at(pin,N-1,N)
        for j in range(N-1): H+=tensor([I]*j+[hh]+[I]*(N-j-2))
        theta=-(N-1)*.2; rot=np.diag([1,np.exp(1j*theta)])
        local=rot.conj().T@M@rot
        g=tensor([plus[:,None]]*N).ravel(); src=at(local,0,N)@g
        e,v=eigh(H); a=abs(v.conj().T@src)**2
        er,ar,_=c.helix_measure(N,M,.2)
        errors=[abs(np.sum(a*np.exp(-t*e))-np.sum(ar*np.exp(-t*er))) for t in (0,.1,1,10)]
        assert max(errors)<2e-8
        hx.append(dict(N=N,full_dimension=2**N,reduced_dimension=N+1,
                       max_heat_kernel_error=max(errors)))
    out['helix_checks']=hx
    # Literal proposed G1 local check and complete 32-dimensional cell.
    import importlib.util
    spec=importlib.util.spec_from_file_location('demonstrator',c.ROOT/'../../contrib/routeA_demonstrator/cell.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    H,ports,G,q=mod.cell_ops(1.,1.,.5)
    ev=np.linalg.eigvalsh(H)
    port_error=max(np.linalg.norm(G.T@o@G) for side in ports.values() for o in side.values())
    assert np.linalg.norm(H@G)<1e-12 and ev[2]>0 and port_error<1e-12
    out['G1_positive_check']=dict(local_dimension=8,cell_dimension=32,
                                 nullity=int(np.sum(abs(ev)<1e-12)),gap=float(ev[2]),
                                 code_residual=float(np.linalg.norm(H@G)),port_error=float(port_error))
    # Literal Pauli-to-Bell source labels.
    bell=np.array([1,0,0,1])/np.sqrt(2)
    S=np.column_stack([bell]+[np.kron(f,I)@bell for f in (X,Y,Z)])
    ops=[S.conj().T@np.kron(f,I)@S for f in (X,Y,Z)]
    vp=sum(np.kron(o,o) for o in ops); source=vp[:,0]
    expected=np.zeros(16); expected[[5,10,15]]=1
    assert np.linalg.norm(source-expected)<1e-12
    out['mixed_term_selection_check']=dict(Bell_pair_dimension=16,
                                          source_labels=[5,10,15],ground_label=0,
                                          source_error=float(np.linalg.norm(source-expected)),
                                          mixed_coefficient=0.)
    # Shortest virtual vacancy paths give the small-w tunneling coefficient.
    paths=[]
    for p in (1,2):
        ch=c.charge(p,.16); cfg=ch['configs']; N=ch['N'];
        start=cfg.index(tuple(range(p))); end=cfg.index(tuple(range(N-p,N)))
        H=ch['H']; adj=[np.flatnonzero(H[j]<-1e-10).tolist() for j in range(len(cfg))]
        dist={start:0}; queue=[start]
        for j in queue:
            for v in adj[j]:
                if v not in dist: dist[v]=dist[j]+1; queue.append(v)
        weight={start:1.}
        for j in queue:
            for v in adj[j]:
                if dist[v]!=dist[j]+1: continue
                bits=[x in cfg[v] for x in range(N)]
                cost=sum(bits[x]!=bits[x+1] for x in range(N-1))-1
                fac=1. if v==end else 1./cost
                weight[v]=weight.get(v,0.)+weight[j]*fac
        paths.append(dict(p=p,min_hops=dist[end],gap_leading_coefficient=2*weight[end]))
    out['charge_gap_leading_paths']=paths
    (c.ROOT/'checks_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
