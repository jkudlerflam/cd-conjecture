"""Independent raw source solves and direct cut-overlap checks."""
import json
import numpy as np
from scipy.linalg import solve
import mobile as m

out={'status':'numerical','cases':[]}
for n in (1,2):
    c=m.core(n,g=.45,mass=.5,express=.2,loop=.08)
    chs={p:m.charge(c['k'],p,.16) for p in (1,2)}
    a,b=chs[1],chs[2]; N=a['N']; k=c['k']
    ix1={x:i for i,x in enumerate(a['cfg'])}; ix2={x:i for i,x in enumerate(b['cfg'])}
    overlap=0.
    for ia,va in enumerate(a['cfg']):
        for ib,vb in enumerate(b['cfg']):
            holes=set(va)|{N+j for j in vb}
            if N not in holes: continue
            left=tuple(sorted(j for j in holes if j<N+1))
            right=tuple(sorted(j-N-1 for j in holes if j>=N+1))
            overlap+=a['chi'][ia]*b['chi'][ib]*b['chi'][ix2[left]]*a['chi'][ix1[right]]
    phi=1/N
    plus=lambda x:np.array([1,np.exp(1j*x)])/np.sqrt(2)
    def phases(lengths):
        return [(j*phi,(j-l+1)*phi) for l in lengths for j in range(l)]
    initial=phases((N,N+1)); final=phases((N+1,N))
    hov=np.exp(-1j*N*phi)
    for j,(old,new) in enumerate(zip(initial,final)):
        for track in (0,1):
            v=plus(old[track])
            if j==N: v=v*np.array([1,np.exp(1j*N*phi)])
            hov*=np.vdot(plus(new[track]),v)
    q=c['Q']
    cr={}
    for side in ('L','R'):
        for f,F in c['F'][side].items():
            dq={'0':0,'+':1,'-':-1}[f]
            cr[side+f]=float(np.linalg.norm(q[:,None]*F-F*q[None,:]-dq*F))
    rec={'n':n,'core_dimension':c['d'],
         'minimum_positive_piece_eigenvalue':float(min(np.linalg.eigvalsh(B)[0] for B in list(c['parts'])+[c['ends']])),
         'charge_source_residuals':cr,
         'occupation_transport_error':float(abs(overlap-m.transport(a,b))),
         'helix_transport_error':float(abs(hov-np.cos(phi/2)**(2*N)))}
    if n==1:
        ch=chs[2]; s=.01*ch['gap']; ex=c['ex']; de=int(ex.sum())
        raw=np.kron((ch['H']-np.linalg.eigvalsh(ch['H'])[0]*np.eye(ch['d']))/s,np.eye(de))
        raw=raw.astype(complex)+np.kron(np.eye(ch['d']),c['ends'][np.ix_(ex,ex)])
        for j,B in enumerate(c['parts']):raw+=np.kron(np.diag(ch['D'][:,j]),B[np.ix_(ex,ex)])
        ev,A,res=m.measure(c,ch,s); errs=[]
        for side,idx in (('L',0),('R',1)):
            for f,F in c['F'][side].items():
                Z=np.kron((ch['occ'][idx]*ch['chi'])[:,None],(F@c['G'])[ex])
                direct=Z.conj().T@solve(raw+np.eye(len(raw)),Z,assume_a='her')
                spectral=A[side][f].conj().T@(A[side][f]/(ev[:,None]+1))
                errs.append(float(np.linalg.norm(direct-spectral,2)))
        rec['raw_source_dimension']=len(raw)
        rec['raw_vs_stabilized_source_absolute_error']=max(errs)
    out['cases'].append(rec)
print(json.dumps(out,indent=2))
open('audit_checks.json','w').write(json.dumps(out,indent=2)+'\n')
