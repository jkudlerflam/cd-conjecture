"""Full microscopic ED and independent squeezing checks, all outputs local."""
from microscopic import *
import json
import time
import scipy.linalg as la
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh


def complete_spectrum(H):
    """Exact invariant graph components, diagonalize every block with multiplicity."""
    adjacency=H.copy();adjacency.data=np.ones(len(adjacency.data))
    count,lab=connected_components(adjacency,directed=False)
    sizes=np.bincount(lab);order=np.argsort(lab,kind='stable')
    starts=np.r_[0,np.cumsum(sizes)[:-1]]
    local=np.empty(len(lab),int);local[order]=np.arange(len(lab))-np.repeat(starts,sizes)
    coo=H.tocoo();assert np.all(lab[coo.row]==lab[coo.col])
    eigen=[]
    for size in np.unique(sizes):
        comps=np.flatnonzero(sizes==size)
        pos=np.full(count,-1,int);pos[comps]=np.arange(len(comps))
        mask=sizes[lab[coo.row]]==size
        blocks=np.zeros((len(comps),size,size),complex)
        blocks[pos[lab[coo.row[mask]]],local[coo.row[mask]],local[coo.col[mask]]]=coo.data[mask]
        eigen.append(np.linalg.eigvalsh(blocks).ravel())
    out=np.sort(np.concatenate(eigen));assert len(out)==H.shape[0]
    return out,dict(components=int(count),largest_block=int(max(sizes)))


def squeezed(h,ow,L,t=1.,J=.02,mu=0.):
    N=len(ow[0]);positions=list(it.combinations(range(L),N));pi={x:i for i,x in enumerate(positions)}
    wi={w:i for i,w in enumerate(ow)};D=len(ow);d=46
    h=h.tocsc();rows=[];cols=[];vals=[];physical=[]
    for ip,x in enumerate(positions):
        for iw,w in enumerate(ow):
            col=ip*D+iw;ww=[0]*L
            for j,s in enumerate(x):ww[s]=w[j]
            physical.append(tuple(ww))
            rows.append(col);cols.append(col);vals.append(-mu*N)
            for j,s in enumerate(x):
                for step in (-1,1):
                    dest=(s+step)%L
                    if dest in x:continue
                    xp=list(x);xp[j]=dest;xp=tuple(sorted(xp));wp=w
                    if s==0 and dest==L-1:wp=w[1:]+w[:1]
                    if s==L-1 and dest==0:wp=w[-1:]+w[:-1]
                    rows.append(pi[xp]*D+wi[wp]);cols.append(col);vals.append(-t)
                k=(j+1)%N
                if (x[k]-s)%L!=1:continue
                c=(w[j]-1)*d+w[k]-1
                for z in range(h.indptr[c],h.indptr[c+1]):
                    rr=h.indices[z];wp=list(w);wp[j]=rr//d+1;wp[k]=rr%d+1;wp=tuple(wp)
                    assert wp in wi
                    rows.append(ip*D+wi[wp]);cols.append(col);vals.append(J*h.data[z])
    mat=sp.coo_matrix((vals,(rows,cols)),shape=(len(physical),len(physical))).tocsr()
    mat.eliminate_zeros();return mat,physical


def sector_check(n):
    letters,charges,h=core(n,readout=0.)
    ow=legal_words(n,letters);N=n+2;L=N+1
    Hsq,words=squeezed(h,ow,L,J=.02)
    Hphys=restricted_physical(lift(h,J=.02,mu=0.),words)
    error=Hphys-Hsq
    err=float(np.max(np.abs(error.data),initial=0))
    assert err<1e-12
    # Every offset orbit contains exactly one word starting with A.
    canonical=[w for w in ow if letters[w[0]-1][0]=='A']
    ci={w:i for i,w in enumerate(canonical)};D=len(canonical)
    def canon(w):
        j=next(j for j,a in enumerate(w) if letters[a-1][0]=='A')
        return w[j:]+w[:j]
    # Q=0 internal shift and zero-twist charge ground state, a uniform one-hole state.
    U=sp.coo_matrix((np.full(len(words),1/np.sqrt(N*L)),
                    (np.arange(len(words)),[ci[canon(tuple(a for a in w if a))] for w in words])),
                   shape=(len(words),D)).tocsr()
    Ui=sp.coo_matrix((np.full(len(ow),1/np.sqrt(N)),
                    (np.arange(len(ow)),[ci[canon(w)] for w in ow])),shape=(len(ow),D)).tocsr()
    # Internal microscopic ring, not a logical chain. Convert one-based occupied labels.
    Hinternal=restricted_physical(h,[tuple(a-1 for a in w) for w in ow])
    H0,_=squeezed(h,ow,L,J=0.)
    V=(Hphys-H0)/.02
    word_index={w:i for i,w in enumerate(words)}
    translated=[]
    for word in words:
        occupied=tuple(a for a in word if a)
        rotated=occupied[-1:]+occupied[:-1]
        ri=iter(rotated)
        translated.append(word_index[tuple(next(ri) if a else 0 for a in word)])
    T=sp.coo_matrix((np.ones(len(words)),(translated,np.arange(len(words)))),
                    shape=Hphys.shape).tocsr()
    commutator=V@T-T@V
    twist_commutator=float(np.sqrt(np.sum(abs(commutator.data)**2)))
    assert twist_commutator>1
    effective=U.getH()@V@U
    target=(N-1)/N*(Ui.getH()@Hinternal@Ui)
    projerr=float(np.max(np.abs((effective-target).data),initial=0))
    residual=float(np.max(np.abs((H0@U+2*U).data),initial=0))
    assert projerr<1e-11 and residual<1e-12
    # One vacancy has L allowed momenta in each internal cyclic sector.
    predicted=[]
    for r in range(N):
        Q=2*np.pi*r/N
        beta=Q+np.pi*(N-1)
        ks=(2*np.pi*np.arange(L)+beta)/L
        # Fill all but one of the L fermion orbitals.
        es=-2*np.cos(ks)
        predicted.extend(np.repeat(es.sum()-es,D))
    ev0,parts=complete_spectrum(H0)
    twist_error=float(np.max(np.abs(ev0-np.sort(predicted))))
    assert twist_error<1e-12
    rng=np.random.default_rng(431+n)
    numerical=[]
    E1=float(la.eigvalsh(target.toarray())[0])
    for J in (.02,.01):
        H=H0+J*V
        values,vectors=eigsh(H,k=2**n+3,which='SA',v0=rng.normal(size=H.shape[0]),tol=2e-11)
        index=np.argsort(values);values=values[index];vectors=vectors[:,index]
        res=float(np.max(np.linalg.norm(H@vectors-vectors*values,axis=0)))
        numerical.append(dict(J=J,energies=values.tolist(),residual=res,
                              first_order_remainder_over_J2=float((values[0]+2-J*E1)/J**2)))
    return dict(n=n,L=L,N=N,dimension=len(words),initialized_cell_dimension=D,
                squeeze_entry_error=err,first_order_entry_error=projerr,
                interaction_internal_shift_commutator_frobenius=twist_commutator,
                charge_ground_residual=residual,J0_full_spectrum_error=twist_error,
                J0_blocks=parts,finite_J=numerical)


def main():
    start=time.monotonic();out={'full_rings':[],'initialized_sectors':[]}
    for n in (1,2):
        letters,q,h=core(n)
        ev,parts=complete_spectrum(h)
        b=h.tocoo();pairq=np.repeat(q,46)+np.tile(q,46)
        comm=float(np.max(np.abs((pairq[b.row]-pairq[b.col])*b.data),initial=0))
        h47=lift(h);H=ring_full3(h47)
        full,blocks=complete_spectrum(H)
        multiplicity=int(np.count_nonzero(abs(full-full[0])<1e-9))
        assert ev[0]>-1e-11 and comm==0
        assert abs(full[0]+2.34)<1e-11 and multiplicity==6
        out['full_rings'].append(dict(n=n,L=3,dimension=47**3,
          bond_min_eigenvalue=float(ev[0]),charge_commutator=comm,
          ground_energy=float(full[0]),ground_multiplicity=multiplicity,
          next_energy=float(full[multiplicity]),**blocks))
        print(out['full_rings'][-1],flush=True)
        out['initialized_sectors'].append(sector_check(n))
        print(out['initialized_sectors'][-1],flush=True)
    out['elapsed_seconds']=time.monotonic()-start
    with open('checks.json','w') as f:json.dump(out,f,indent=2)


if __name__=='__main__':main()
