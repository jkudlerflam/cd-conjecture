"""Open-chain DMRG of the full parity-class qubit chain in the Fibonacci cusp (C152 bond convention).
Maps the effective-model point (|w|, phi, mu_eff) at given V/Delta to bond couplings:
  Delta = Jx - Jy = 1, V = 4 Jz, |t2| = Delta^2/(2V), w = |w||t2| e^{i phi} = Jx + Jy + 2iD,
  mu = 4Jz + 2hz = |t2| (mu_eff + 6).
Edge sites get the missing half field as an on-site term."""
import numpy as np, json, sys, time
from dmrg_d import dmrg
SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]]); SZ=np.diag([1.,-1.]).astype(complex); I2=np.eye(2)
def bond(Jx,Jy,D,Jz,hz):
    return (Jx*np.kron(SX,SX)+Jy*np.kron(SY,SY)+D*(np.kron(SX,SY)-np.kron(SY,SX))+Jz*np.kron(SZ,SZ)
            +0.5*hz*(np.kron(SZ,I2)+np.kron(I2,SZ)))
def couplings(V,wa,phi,mu_eff,Dl=1.0):
    t2=Dl**2/(2*V); w=wa*t2*np.exp(1j*phi)
    Jz=V/4; JxpJy=w.real; D=w.imag/2; Jx=(JxpJy+Dl)/2; Jy=(JxpJy-Dl)/2
    mu=t2*(mu_eff+6); hz=(mu-4*Jz)/2
    return Jx,Jy,D,Jz,hz
if __name__=="__main__":
    a=sys.argv[2]; cfgs=json.load(open(a)) if a.endswith('.json') else json.loads(a); f=open(f"fd_{sys.argv[1]}.jsonl","a")
    for c in cfgs:
        t0=time.time()
        Jx,Jy,D,Jz,hz=couplings(c['V'],c['wa'],c['phi'],c['mu'])
        h=bond(Jx,Jy,D,Jz,hz); L=c['L']
        ons=[None]*L; ons[0]=0.5*hz*SZ; ons[-1]=0.5*hz*SZ
        out=dmrg(h,L,chis=tuple(c['chis']),cutoff=1e-12,onsite=ons,lanczos_tol=1e-12)
        S=out['S']
        rec=dict(c,couplings=[Jx,Jy,D,Jz,hz],E=out['E'],S=[float(x) for x in S],Smid=float(S[L//2-1]),hist=out['hist'],secs=time.time()-t0)
        f.write(json.dumps(rec)+"\n"); f.flush(); print(json.dumps({k:rec[k] for k in rec if k!='S'}),flush=True)
