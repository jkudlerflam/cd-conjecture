"""Connected density/bond matrices and pair contractions in a fixed bulk window."""
import time
import numpy as np
from phase_common import Z


def block_tensor(psi,i,width):
    t=psi.M[i]
    for j in range(i+1,i+width):t=np.tensordot(t,psi.M[j],axes=([-1],[0]))
    return t.reshape(t.shape[0],2**width,t.shape[-1])


def insert(E,T,O):
    u=np.tensordot(E,T,axes=([1],[0]))
    u=np.tensordot(u,O,axes=([1],[1])).transpose(0,2,1)
    return np.tensordot(T.conj(),u,axes=([0,1],[0,1]))


def right_effect(T,O):
    u=np.tensordot(T,O,axes=([1],[1])).transpose(0,2,1)
    return T.conj().reshape(T.shape[0],-1)@u.reshape(T.shape[0],-1).T


def close(E,R):return np.sum(E*R)


def environments(psi):
    left=[np.ones((1,1))]
    for A in psi.M:left.append(Z._apply_site(left[-1],A,None))
    assert abs(left[-1][0,0]-1)<1e-8
    return left


def operator_matrix(psi,left,ops,width,indices):
    """All profiles plus symmetrized raw covariance on selected supports."""
    effects=[];means=[]
    for j,O in enumerate(ops):
        R=right_effect(block_tensor(psi,j,width),O)
        effects.append(R);means.append(float(np.real(close(left[j],R))))
    raw=np.zeros((len(indices),len(indices)))
    for ii,i in enumerate(indices):
        T=block_tensor(psi,i,width)
        raw[ii,ii]=np.real(close(left[i],right_effect(T,ops[i]@ops[i])))
        E=insert(left[i],T,ops[i]);cursor=i+width
        for jj in range(ii+1,len(indices)):
            j=indices[jj]
            if j<i+width:
                total=width+j-i
                O=np.kron(ops[i],np.eye(2**(j-i)))@np.kron(np.eye(2**(j-i)),ops[j])
                val=close(left[i],right_effect(block_tensor(psi,i,total),O))
            else:
                while cursor<j:
                    E=Z._apply_site(E,psi.M[cursor],None);cursor+=1
                val=close(E,effects[j])
            raw[ii,jj]=raw[jj,ii]=np.real(val)
    m=np.array(means)[indices];conn=raw-np.outer(m,m)
    assert np.min(np.linalg.eigvalsh(conn))>-1e-8
    return dict(profile=means,indices=list(indices),raw=raw.tolist(),connected=conn.tolist())


def spectral_metrics(profile,connected,indices):
    x=np.array(indices,float);m=np.array(profile)[indices];M=len(m)
    delta=m-m.mean();C=np.array(connected)
    qs=np.linspace(0,np.pi,4097);phase=np.exp(1j*np.outer(qs,x-x[0]))
    amp=phase@delta
    sd=np.abs(amp)**2/M
    sc=np.real(np.einsum('qi,ij,qj->q',phase.conj(),C,phase,optimize=True))/M
    assert sc.min()>-1e-8
    allowed=np.flatnonzero(qs>=2*np.pi/M)
    peakd=allowed[np.argmax(sd[allowed])];peakc=allowed[np.argmax(sc[allowed])]
    def peak(k):return dict(q=float(qs[k]),q_over_pi=float(qs[k]/np.pi),
        amplitude=float(2*abs(amp[k])/M),S_disconnected=float(sd[k]),S_connected=float(sc[k]),
        disconnected_over_M=float(sd[k]/M),connected_over_M=float(sc[k]/M),
        total_over_M=float((sd[k]+sc[k])/M))
    distances=[dict(r=r,connected_mean=float(np.diag(C,r).mean()),
        connected_RMS=float(np.sqrt(np.mean(np.diag(C,r)**2)))) for r in range(M)]
    return dict(M=M,mean=float(m.mean()),range=float(np.ptp(m)),RMS=float(np.sqrt(np.mean(delta**2))),
        profile_peak=peak(peakd),connected_peak=peak(peakc),q_grid=qs[::16].tolist(),
        S_connected=sc[::16].tolist(),S_disconnected=sd[::16].tolist(),distance_average=distances)


def pair_correlations(psi,left,gph,lo,hi):
    O=np.kron(np.kron(Z.SP,Z.I2),Z.SP);effects={};means={};phases={}
    for j in range(lo,hi-2):
        effects[j]=right_effect(block_tensor(psi,j,3),O)
        phases[j]=gph[j]*gph[j+2]
        means[j]=phases[j]*close(left[j],effects[j])
    result=[]
    for i in range(lo,min(lo+3,hi-6)):
        E=insert(left[i],block_tensor(psi,i,3),O.T);cursor=i+3;row=[]
        for j in range(i+4,hi-2):
            while cursor<j:E=Z._apply_site(E,psi.M[cursor],None);cursor+=1
            raw=np.conj(phases[i])*phases[j]*close(E,effects[j])
            disconnected=np.conj(means[i])*means[j];conn=raw-disconnected
            row.append(dict(r=j-i,j=j,raw=[float(raw.real),float(raw.imag)],
                disconnected=[float(disconnected.real),float(disconnected.imag)],
                connected=[float(conn.real),float(conn.imag)]))
        result.append(dict(i=i,mean_pair=[float(means[i].real),float(means[i].imag)],correlations=row))
    return result


def observe(psi,gph,bonds):
    start=time.time();L=psi.L;lo=L//4;hi=3*L//4
    left=environments(psi);n=np.diag([0.,1.]);zz=np.kron(Z.ZZ_,Z.ZZ_)
    out=dict(L=L,P=psi.P,bulk_half=[lo,hi],bulk_quarter=[3*L//8,5*L//8],operators={})
    for name,ops,width in [('density',[n]*L,1),('bond',bonds,2),('kinetic',[b-zz for b in bonds],2)]:
        indices=list(range(lo,hi-width+1));d=operator_matrix(psi,left,ops,width,indices)
        d['definition']={'density':'n=(1-Z)/2','bond':'physical interaction bond, excluding onsite fields',
                         'kinetic':'transverse bond B-ZZ'}[name]
        d['half']=spectral_metrics(d['profile'],d['connected'],indices)
        qlo,qhi=out['bulk_quarter'];ix=[k for k,j in enumerate(indices) if qlo<=j and j+width<=qhi]
        qindices=[indices[k] for k in ix];C=np.array(d['connected'])[np.ix_(ix,ix)]
        d['quarter']=spectral_metrics(d['profile'],C,qindices)
        out['operators'][name]=d
    out['pair']=pair_correlations(psi,left,gph,lo,hi)
    out['runtime_seconds']=time.time()-start
    out['conventions']='Bond matrices are symmetrized; raw, connected and profile-disconnected Fourier terms are distinct. Pair A_j=b_j b_(j+2).'
    return out
