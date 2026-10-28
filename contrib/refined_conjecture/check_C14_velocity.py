# Independent check by the coordinating Claude session of Astra C14 (codex/work_C14): spectrum of the SU(N) rank-2 bond and its velocity.
# Independent check (Claude): velocity of the SU(N) rank-2 symmetric chain h = -X^2/8 + 3X/4 - 1.
import itertools, sys, numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
N=int(sys.argv[1]); Ls=[int(x) for x in sys.argv[2:]]
loc=[s for s in itertools.product(range(3),repeat=N) if sum(s)==2]; d=len(loc); idx={s:i for i,s in enumerate(loc)}
E=np.zeros((N,N,d,d))
for a in range(N):
  for b in range(N):
    for j,s in enumerate(loc):
      if s[b]==0: continue
      t=list(s); amp=np.sqrt(t[b]); t[b]-=1; amp*=np.sqrt(t[a]+1); t[a]+=1
      if tuple(t) in idx: E[a,b,idx[tuple(t)],j]+=amp
X=sum(np.kron(E[a,b],E[b,a]) for a in range(N) for b in range(N))
h=-X@X/8+3*X/4-np.eye(d*d); ev=np.linalg.eigvalsh(h); print("N",N,"d",d,"spec h",np.unique(np.round(ev,8)))
col=np.array(loc)
for L in Ls:
  tot=2*L
  if tot%N: print("skip L",L); continue
  target=[tot//N]*N
  confs=[c for c in itertools.product(range(d),repeat=L) if list(col[list(c)].sum(0))==target]
  cidx={c:i for i,c in enumerate(confs)}; D=len(confs)
  hs=sparse.coo_matrix(h); rows=[];cols=[];vals=[]
  for i,c in enumerate(confs):
    for j in range(L):
      k=(j+1)%L; p=c[j]*d+c[k]
      sel=hs.col==p
      for r,v in zip(hs.row[sel],hs.data[sel]):
        nc=list(c); nc[j],nc[k]=divmod(r,d); t=cidx.get(tuple(nc))
        if t is not None: rows.append(t); cols.append(i); vals.append(v)
  H=sparse.csr_matrix((vals,(rows,cols)),shape=(D,D))
  nev=min(30,D-2); w,V=eigsh(H,k=nev,which='SA'); o=np.argsort(w); w=w[o]; V=V[:,o]
  perm=np.array([cidx[tuple(c[1:]+c[:1])] for c in confs])
  # momentum resolution within near-degenerate blocks
  res=[]; i=0
  while i<len(w):
    j=i
    while j+1<len(w) and abs(w[j+1]-w[i])<1e-8: j+=1
    B=V[:,i:j+1]; T=B.conj().T@B[perm] if False else B.conj().T@B[np.argsort(perm)]
    lam=np.linalg.eigvals(T); ks=np.angle(lam)*L/(2*np.pi)
    for kk in ks: res.append((w[i],kk))
    i=j+1
  E0=res[0][0]
  k1=[e for e,kk in res if abs(abs(kk)-1)<0.3]
  print(f"L={L} D={D} E0/L={E0/L:.6f} E0={E0:.6f}", "lowest |k|=2pi/L gap", (min(k1)-E0) if k1 else None, "v_est", (L*(min(k1)-E0)/(2*np.pi)) if k1 else None, "target v", 2*np.pi/N)
