"""Sparse periodic ED in a conserved charge sector, with two estimators.
Calibrates the bond normalization against ULS and spin-1 TB.
Writes ed_calibration.json. Largest Hilbert space 3**12, sector 73789.
"""
import json, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from fusion_checks import bond
from fusion_tba import ground_energy

def solve(h,d,L,charges,target):
    t0=time.time(); power=d**np.arange(L,dtype=np.int64)
    full=np.arange(d**L,dtype=np.int64)
    digits=(full[:,None]//power[None,:])%d
    keep=np.sum(np.asarray(charges)[digits],axis=1)==target
    states=full[keep]; ds=digits[keep]; lookup=np.full(d**L,-1,dtype=np.int64); lookup[states]=np.arange(len(states))
    rows=[]; cols=[]; vals=[]; transitions=np.argwhere(abs(h)>1e-12)
    for site in range(L):
      j=(site+1)%L
      for out,inp in transitions:
        a,b=divmod(int(inp),d); aa,bb=divmod(int(out),d)
        ids=np.flatnonzero((ds[:,site]==a)&(ds[:,j]==b))
        dst=states[ids]+(aa-a)*power[site]+(bb-b)*power[j]
        rr=lookup[dst]; assert np.all(rr>=0)
        rows.append(rr); cols.append(ids); vals.append(np.full(len(ids),h[out,inp]))
    H=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(states),len(states))).tocsr()
    ev,U=eigsh(H,k=1,which='SA',tol=2e-11,v0=np.random.default_rng(71).normal(size=len(states)))
    psi=np.zeros(d**L); psi[states]=U[:,0]
    cuts=range(2,L//2+1); ent=[]
    for cut in cuts:
      sv=np.linalg.svd(psi.reshape(d**cut,d**(L-cut)),compute_uv=False); p=sv**2; p=p[p>1e-14]
      ent.append(float(-p@np.log(p)))
    x=np.log(L/np.pi*np.sin(np.pi*np.array(list(cuts))/L))
    c_ent=3*np.polyfit(x,ent,1)[0] if len(x)>1 else None
    return dict(L=L,dim_sector=len(states),energy=float(ev[0]),cuts=list(cuts),entropy=ent,
                c_spatial_entropy=c_ent,residual=float(np.linalg.norm(H@U[:,0]-ev[0]*U[:,0])),seconds=time.time()-t0)

def main():
    allout=[]
    for N,k,Ls in [(2,1,[8,10,12,14,16]),(2,2,[6,8,10,12]),(3,1,[6,9])]:
      _,h,*_=bond(N,k); d=int(np.sqrt(len(h))); rows=[]
      # For SU(3), base-N charges distinguish all populations at these lengths.
      for L in Ls:
        charges=list(range(d)) if N==2 else [0,1,L+1]
        target=L*(d-1)//2 if N==2 else (L//3)*(L+2)
        row=solve(h,d,L,charges,target); rows.append(row); print(N,k,row,flush=True)
        e0=ground_energy(N,k); v=2*np.pi/N
        row['c_casimir_known_bulk']=-6*L*(row['energy']-e0*L)/(np.pi*v)
      ll=np.array([r['L'] for r in rows],float); ee=np.array([r['energy'] for r in rows])/ll
      X=np.array([np.ones_like(ll),-np.pi/(6*ll**2),ll**-4]).T
      fit=np.linalg.lstsq(X,ee,rcond=None)[0] if len(rows)>=3 else None
      out=dict(N=N,k=k,d=d,c_exact=k*(N*N-1)/(N+k),v_exact=2*np.pi/N,rows=rows)
      if fit is not None: out['fit']=dict(e0=float(fit[0]),A=float(fit[1]),c=float(fit[1]/(2*np.pi/N)))
      allout.append(out); open('ed_calibration.json','w').write(json.dumps(allout,indent=2))
if __name__=='__main__':main()
