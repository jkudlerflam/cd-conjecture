"""Extract g^2 directly by three resolvent solves, without finite differences."""
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import mobile as m
from band_checks import gram

out=[]
for n in (1,2):
    Gm,_=gram(n,1.); k=n+2
    for D in (10.,100.,1000.):
        c=m.core(n,g=0,mass=1.,express=D,loop=0,phase=np.pi/(2*n))
        cg=m.core(n,g=1,mass=1.,express=D,loop=0,phase=np.pi/(2*n))
        ex=c['ex']; H=sp.csc_matrix(c['H'][np.ix_(ex,ex)])
        C=sp.csc_matrix((cg['H']-c['H'])[np.ix_(ex,ex)])
        I=sp.eye(H.shape[0],format='csc')
        HH=sp.kron(H,I,format='csc')+sp.kron(I,H,format='csc')
        CC=sp.kron(C,I,format='csc')+sp.kron(I,C,format='csc')
        Z=sum(v*np.kron((c['F']['R'][fr]@c['G'])[ex],(c['F']['L'][fl]@c['G'])[ex]) for fr,fl,v in m.TERMS)
        lu=splu(HH); X=lu.solve(Z); X=lu.solve(CC@X); X=lu.solve(CC@X)
        K=-Z.conj().T@X
        ids=[2**(2*n-1-j) for j in range(2*n)]
        T=K[np.ix_(ids[:n],ids[n:])]
        normalized=-T*D**(2*n+4)*k*k/2
        out.append({'n':n,'express':D,'pair_source_dimension':HH.shape[0],
                    'relative_Gram_error':float(np.linalg.norm(normalized-Gm,2)/np.linalg.norm(Gm,2)),
                    'eigenvalues_normalized_hermitian_part':np.linalg.eigvalsh((normalized+normalized.conj().T)/2).tolist()})
print(json.dumps(out,indent=2))
open('asymptotic_checks.json','w').write(json.dumps(out,indent=2)+'\n')
