"""C22: exact one-bus resolvent of a local, empty-bus ground code.

Only NumPy/SciPy. The one-bus restriction is exact for H2, not a
finite-coupling truncation claim. Physical many-bus ED is separate.
"""
import json
import numpy as np
from scipy.linalg import eigh


def kron_all(xs):
    ans = np.ones((1, 1), complex)
    for x in xs:
        ans = np.kron(ans, x)
    return ans


def data_ops(n):
    I = np.eye(2)
    Z = np.diag([1., -1.])
    sm = np.array([[0., 1.], [0., 0.]])
    cs, zs, spins = [], [], []
    for j in range(n):
        cs.append(kron_all([Z if i < j else sm if i == j else I for i in range(n)]))
        spins.append(kron_all([sm if i == j else I for i in range(n)]))
        zs.append(kron_all([Z if i == j else I for i in range(n)]))
    return cs, zs, spins, kron_all([Z]*n)


def one_bus(n, delta=4., t=1., g=.1, phases=(0., 0., 0.), gauged=True, minus_ratio=1.):
    """Species 0,+,-. Site n is a marker without a data qubit."""
    N, m = 2**n, n+1
    cs, zs, spins, P = data_ops(n)
    H = delta*np.eye(3*m*N, dtype=complex)
    def sl(j, a):
        return slice((3*j+a)*N, (3*j+a+1)*N)
    for j in range(n):
        for a in range(3):
            hop = np.eye(N) if gauged or a != 0 else zs[j]
            val = -t*np.exp(1j*phases[a])*hop
            H[sl(j+1,a), sl(j,a)] += val
            H[sl(j,a), sl(j+1,a)] += val.conj().T
        c = cs[j] if gauged else spins[j]
        H[sl(j,1),sl(j,0)] += g*c
        H[sl(j,0),sl(j,1)] += g*c.conj().T
        H[sl(j,2),sl(j,0)] += g*minus_ratio*c.conj().T
        H[sl(j,0),sl(j,2)] += g*minus_ratio*c
    sources = []
    for j in (0,n):
        row = []
        for a in range(3):
            S = np.zeros((len(H), N), complex)
            S[sl(j,a)] = P if gauged and j == n and a == 0 else np.eye(N)
            row.append(S)
        sources.append(row)
    return H, sources


def h2(n, link_phase=np.pi/2, **kw):
    H, (L, R) = one_bus(n, **kw)
    e, U = eigh(H)
    A = [U.conj().T @ s for s in R]
    B = [U.conj().T @ s for s in L]
    W = (np.kron(A[0], B[0])+np.exp(1j*link_phase)*np.kron(A[1], B[2])
         +np.exp(-1j*link_phase)*np.kron(A[2], B[1]))
    inv = 1/(e[:,None]+e[None,:]).reshape(-1)
    K = -W.conj().T @ (inv[:,None]*W)
    return K, e


def h2_g2(n, link_phase=np.pi/2, **kw):
    """Exact coefficient of g^2 in H2, before endpoint visibility."""
    kw.pop('g', None)
    A, (L,R) = one_bus(n,g=0.,**kw)
    H1,_ = one_bus(n,g=1.,**kw)
    B = H1-A
    e,U=eigh(A)
    B=U.conj().T@B@U
    L=[U.conj().T@s for s in L]
    R=[U.conj().T@s for s in R]
    W=(np.kron(R[0],L[0])+np.exp(1j*link_phase)*np.kron(R[1],L[2])
       +np.exp(-1j*link_phase)*np.kron(R[2],L[1]))
    den=e[:,None]+e[None,:]
    X=(W/den.reshape(-1,1)).reshape(len(A),len(A),-1)
    Y=np.einsum('ab,bck->ack',B,X)+np.einsum('ab,cbk->cak',B,X)
    Y=Y.reshape(len(A)**2,-1)
    return -Y.conj().T@(Y/den.reshape(-1,1))


def walker_cell(n, delta=4., t=1., g=.1, phases=(0.,0.,0.), t0=1.):
    """Complete legal cell. Vacuum walker sector plus three excited sectors.

    No truncation: a legal cell has one token by the local arrow grammar.
    The vacuum ground energy is subtracted. `delta` is excited-sector
    diagonal energy after that subtraction.
    """
    N,m=2**n,n+1
    A=-t0*(np.diag(np.ones(m-1),1)+np.diag(np.ones(m-1),-1))
    ev,uv=eigh(A)
    Hv=np.kron(A-ev[0]*np.eye(m),np.eye(N))
    He,sources=one_bus(n,delta,t,g,phases,gauged=True)
    from scipy.linalg import block_diag
    H=block_diag(Hv,He)
    G=np.zeros((len(H),N),complex)
    G[:m*N]=np.kron(uv[:,0,None],np.eye(N))
    if uv[0,0]<0:
        G=-G
    ports=[]
    for end,j in enumerate((0,n)):
        S=sources[end]
        V=np.zeros((m*N,N),complex)
        V[j*N:(j+1)*N]=np.eye(N)
        forward=[]
        for a in range(3):
            F=np.zeros_like(H)
            F[m*N:,:m*N]=S[a]@V.conj().T
            forward.append(F)
        F0=forward[0]+forward[0].conj().T
        Fp=forward[1]+forward[2].conj().T
        ports.append((F0,Fp,Fp.conj().T))
    return H,G,ports,float(abs(uv[0,0]))


def quadratic_fit(K, modes):
    cs, _, _, _ = data_ops(modes)
    D = len(K)
    majoranas = [x for c in cs for x in (c+c.conj().T, -1j*(c-c.conj().T))]
    scalar = np.trace(K).real/D
    fit = scalar*np.eye(D, dtype=complex)
    coeff = np.zeros((2*modes,2*modes))
    for a in range(2*modes):
        for b in range(a+1,2*modes):
            op = 1j*majoranas[a]@majoranas[b]
            v = np.trace(op@K).real/D
            coeff[a,b] = v
            fit += v*op
    # number-conserving one-body matrix, after removing the scalar vacuum
    vac = 0
    inds = [1 << (modes-1-j) for j in range(modes)]
    one = K[np.ix_(inds,inds)]-K[vac,vac]*np.eye(modes)
    return dict(scalar=scalar, residual=float(np.linalg.norm(K-fit,2)),
                quadratic_norm=float(np.linalg.norm(fit-scalar*np.eye(D),2))), one, coeff


def diagnostic(n, g=.1, phases=(0.,0.,0.)):
    K,e = h2(n,g=g,phases=phases)
    fit, one, _ = quadratic_fit(K,2*n)
    cs,_,_,P = data_ops(2*n)
    charge = sum(c.conj().T@c for c in cs)
    fit.update(n=n,g=g,phases=phases,gap=float(e[0]),
               parity_error=float(np.linalg.norm(K@P-P@K)),
               charge_error=float(np.linalg.norm(K@charge-charge@K)),
               hopping_real=one[:n,n:].real.tolist(), hopping_imag=one[:n,n:].imag.tolist(),
               field_left=one[:n,:n].real.tolist(), field_right=one[n:,n:].real.tolist())
    M = one[:n,:n]+one[n:,n:]
    T = one[:n,n:]
    bands=np.array([np.linalg.eigvalsh(M+np.exp(1j*k)*T+np.exp(-1j*k)*T.conj().T)
                    for k in np.linspace(-np.pi,np.pi,4001)])
    fit.update(band_min=bands.min(axis=0).tolist(),band_max=bands.max(axis=0).tolist(),
               hopping_singular_values=np.linalg.svd(T,compute_uv=False).tolist())
    return fit


if __name__ == '__main__':
    out=[]
    for n in (1,2,3):
        for g in (.2,.1,.05):
            d=diagnostic(n,g)
            out.append(d)
            print(json.dumps(d),flush=True)
    for n in (1,2):
        d=diagnostic(n,.1,(.13,.31,-.07))
        out.append(d)
        print(json.dumps(d),flush=True)
    with open('bus_diagnostics.json','w') as f:
        json.dump(out,f,indent=2)
