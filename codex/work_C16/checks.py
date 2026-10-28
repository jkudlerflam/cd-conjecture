"""Small C16 checks. Only NumPy/SciPy and the standard library; no network."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import itertools
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

P = [np.array([[0, 1], [1, 0]], complex),
     np.array([[0, -1j], [1j, 0]], complex),
     np.diag([1., -1.]).astype(complex)]
I = np.eye(2)

def local_bond(J, B):
    return sum(J[a,b]*np.kron(P[a], P[b]) for a in range(3) for b in range(3)) + sum(B[a]*(np.kron(P[a], I)+np.kron(I, P[a]))/2 for a in range(3))

def ring_matrix(h, L):
    rows, cols, vals = [], [], []
    for s in range(1 << L):
        for i in range(L):
            j = (i+1) % L
            x = 2*((s >> i)&1)+((s >> j)&1)
            for y in range(4):
                val = h[y,x]
                if abs(val) > 1e-14:
                    target = (s & ~(1 << i) & ~(1 << j)) | ((y >> 1) << i) | ((y & 1) << j)
                    rows.append(target); cols.append(s); vals.append(val)
    return coo_matrix((vals, (rows, cols)), shape=(1 << L,)*2).tocsr()

def allowed(s, L):
    return (s & ((s << 1) | (s >> (L-1)))) == 0

def blockade_basis(L, N=None):
    if N is None:
        return [s for s in range(1 << L) if allowed(s,L)]
    return [s for positions in itertools.combinations(range(L),N)
            if allowed(s := sum(1 << i for i in positions), L)]

def constrained_operators(L, basis):
    index = {s: i for i,s in enumerate(basis)}
    ops = {key: np.zeros((len(basis),)*2, complex) for key in ('N','T','C','X','Y')}
    for col,s in enumerate(basis):
        ops['N'][col,col] = s.bit_count()
        for i in range(L):
            target = s ^ (1 << i)
            if target in index:
                row = index[target]
                ops['X'][row,col] += 1
                ops['Y'][row,col] += 1j*(1-2*((s >> i)&1))
            j = (i+1)%L
            ni,nj = (s >> i)&1,(s >> j)&1
            target = s ^ (1 << i) ^ (1 << j)
            if ni != nj and target in index:
                row = index[target]
                ops['T'][row,col] += 1
                ops['C'][row,col] += 1j*(ni-nj)
    return ops

def hopping(L,N):
    basis = blockade_basis(L,N)
    index = {s: i for i,s in enumerate(basis)}
    rows,cols,vals = [],[],[]
    for col,s in enumerate(basis):
        for i in range(L):
            j=(i+1)%L
            if ((s >> i)&1) != ((s >> j)&1):
                target = s ^ (1 << i) ^ (1 << j)
                if target in index:
                    rows.append(index[target]); cols.append(col); vals.append(-1.)
    H=coo_matrix((vals,(rows,cols)),shape=(len(basis),)*2).tocsr()
    E,vec=eigsh(H,k=1,which='SA',tol=1e-11,v0=np.ones(len(basis)))
    return basis, E[0], vec[:,0], np.linalg.norm(H@vec[:,0]-E[0]*vec[:,0])

def entropy(basis, psi, ell):
    mask=(1 << ell)-1
    left=sorted(set(s & mask for s in basis)); right=sorted(set(s >> ell for s in basis))
    li={s:i for i,s in enumerate(left)}; ri={s:i for i,s in enumerate(right)}
    mat=np.zeros((len(left),len(right)))
    for s,a in zip(basis,psi):
        mat[li[s & mask],ri[s >> ell]]=a
    probs=np.linalg.svd(mat,compute_uv=False)**2
    probs=probs[probs>1e-15]
    return float(-np.dot(probs,np.log(probs)))

out={}
L=7; a,b,D,Bz=.7,-.2,.31,-2.3
J=np.array([[a,D,0],[-D,a,0],[0,0,b]])
B=np.array([0.,0.,Bz])
H=ring_matrix(local_bond(J,B),L)
ones=[1 << i for i in range(L)]
E0=float(H[0,0].real)
exact=np.sort(-4*b-2*Bz+4*a*np.cos(2*np.pi*np.arange(L)/L)+4*D*np.sin(2*np.pi*np.arange(L)/L))
out['magnon']={'L':L,'dispersion_error':float(np.max(np.abs(np.linalg.eigvalsh(H[ones,:][:,ones].toarray())-E0-exact)))}
J[0,2]=.13; J[2,0]=-.13
Hmix=ring_matrix(local_bond(J,B),L)
other=[s for s in range(1 << L) if s.bit_count()!=1]
out['product_not_magnon_invariant']={'vacuum_residual':float(np.linalg.norm((Hmix[:,0]).toarray()[1:])), 'one_magnon_leakage_norm':float(np.linalg.norm(Hmix[other,:][:,ones].toarray()))}
ks=np.array([.1,.05,.025,.0125])
xi=-2*np.sin(ks/2)**2
slow=xi**2/(np.sqrt(xi**2+np.sin(ks)**2)+np.sin(ks))
out['cubic_edge']={'k':ks.tolist(),'positive_slow_energy':slow.tolist(),'effective_z':(np.log(slow[:-1]/slow[1:])/np.log(2)).tolist(),'energy_over_k_cubed':(slow/ks**3).tolist()}

out['blockade_dimensions']=[{'L':L,'dimension':len(blockade_basis(L))} for L in range(4,13)]
L=8; rng=np.random.default_rng(16027)
J=rng.normal(size=(3,3)); B=rng.normal(size=3)
H=ring_matrix(local_bond(J,B),L)
basis=blockade_basis(L); o=constrained_operators(L,basis)
pred=L*(J[2,2]+B[2])*np.eye(len(basis))+(-4*J[2,2]-2*B[2])*o['N']+(J[0,0]+J[1,1])*o['T']+(J[0,1]-J[1,0])*o['C']+(B[0]+J[0,2]+J[2,0])*o['X']+(B[1]+J[1,2]+J[2,1])*o['Y']
out['blockade_projection']={'L':L,'dimension':len(basis),'max_entry_error':float(np.max(np.abs(H[basis,:][:,basis].toarray()-pred)))}

rows=[]
for L in (12,15,18,21,24):
    N=L//3
    basis,E,psi,residual=hopping(L,N)
    half=entropy(basis,psi,L//2)
    cuts=np.arange(3,L//2+1)
    es=np.array([entropy(basis,psi,int(cut)) for cut in cuts])
    logchord=np.log(L/np.pi*np.sin(np.pi*cuts/L))
    cprofile=3*np.polyfit(logchord,es,1)[0]
    row={'L':L,'N':N,'dimension':len(basis),'energy':float(E),'energy_formula':float(-2/np.sin(np.pi/(L-N))), 'eigen_residual':float(residual),'entropy_half':half,'c_entropy_profile':float(cprofile),'c_casimir_exact_bulk_and_velocity':float(-6*L*(E+4*L/(3*np.pi))/(np.pi*3))}
    rows.append(row)
    print(json.dumps(row),flush=True)
Lvals=np.array([r['L'] for r in rows],float)
S=np.array([r['entropy_half'] for r in rows])
x=np.log(Lvals/np.pi*np.sin(np.pi*np.floor(Lvals/2)/Lvals))
# Pure leading fit and a fit allowing a 1/L^2 correction, reported separately.
out['blockade_hopping']={'rho':1/3,'velocity':3.,'K':4/9,'rows':rows,'c_half_size_fit':float(3*np.polyfit(x,S,1)[0]),'c_half_size_fit_with_Lminus2':float(3*np.linalg.lstsq(np.array([x,np.ones(len(x)),Lvals**-2]).T,S,rcond=None)[0][0])}
with open('checks_results.json','w') as f:
    json.dump(out,f,indent=2)
print(json.dumps({k:v for k,v in out.items() if k!='blockade_hopping'},indent=2))
print('Half-size entropy fits:',out['blockade_hopping']['c_half_size_fit'],out['blockade_hopping']['c_half_size_fit_with_Lminus2'])
