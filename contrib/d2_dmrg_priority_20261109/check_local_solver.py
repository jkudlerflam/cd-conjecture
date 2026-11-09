"""Diagnose Davidson convergence on a fixed central two-site effective problem."""
import argparse,json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parents[1]/'codex/work_C183'))
import numpy as np
import dmrg2 as D
import band as B

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--output',required=True);a=ap.parse_args();d=json.loads(Path(a.input).read_text());L=d['L'];cp=d['couplings'];j=L//2-1
    z=np.load(d['mps_checkpoint']);M=[z[f'M{i}'].copy() for i in range(L)]
    if np.max(abs(z['gauge_angles']))>1e-12:raise ValueError('diagnostic currently handles the physical complex gauge only')
    for i in range(L-1,0,-1):
        A=M[i].reshape(M[i].shape[0],-1);Q,R=np.linalg.qr(A.T.conj());M[i]=Q.T.conj().reshape(-1,2,M[i].shape[2]);M[i-1]=np.tensordot(M[i-1],R.T.conj(),axes=([2],[0]))
    M[0]/=np.linalg.norm(M[0])
    for i in range(j):
        A=M[i].reshape(-1,M[i].shape[2]);Q,R=np.linalg.qr(A);M[i]=Q.reshape(M[i].shape[0],2,-1);M[i+1]=np.tensordot(R,M[i+1],axes=([1],[0]))
    h=B.bond(**cp)-(cp['Jz']+cp['hz'])*np.eye(4);Ws,_=D.bond_mpo_list([h]*(L-1));edge=cp['hz']*B.SZ/2
    Ws[0][0,-1]+=edge;Ws[-1][0,0]+=edge;Ws=D.add_parity_string(Ws,-1 if d['parity']=='even' else 1)
    le=np.ones((1,1,1),complex);re=np.ones_like(le)
    for i in range(j):le=D._lupdate(le,M[i],Ws[i])
    for i in range(L-1,j+1,-1):re=D._rupdate(re,M[i],Ws[i])
    th=np.tensordot(M[j],M[j+1],axes=([2],[0]));mv=D._heff_apply(le,Ws[j],Ws[j+1],re,th.shape);diag=D._heff_diag(le,Ws[j],Ws[j+1],re)
    out=dict(input=a.input,L=L,parity=d['parity'],bond=j,dimension=th.size,trials=[])
    for mmax in [24,48,80]:
        t=time.time();e,v,n,res=D.davidson(mv,th.reshape(-1),diag,res_tol=1e-8,maxit=400,mmax=mmax)
        direct=np.linalg.norm(mv(v)-e*v)
        r=dict(mmax=mmax,E=float(e),matvecs=n,reported_residual=float(res),direct_residual=float(direct),seconds=time.time()-t);out['trials'].append(r);print(r,flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
