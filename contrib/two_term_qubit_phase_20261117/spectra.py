"""One-thread, bounded-size diagnostics. Run with OPENBLAS_NUM_THREADS=1."""
import argparse,json,time,os,fcntl,tempfile
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from model import Block


def solve(L,g,k=16):
    start=time.perf_counter();blocks=[];ground=None
    for p,q in ((1,1),(1,-1),(-1,1)):
        b=Block(L,p,q);A,D=b.operators();H=A+g*D
        rng=np.random.default_rng(20261117+L)
        E,V=eigsh(H,k=min(k,H.shape[0]-2),which='SA',tol=2e-12,ncv=min(max(4*k,64),H.shape[0]),v0=rng.normal(size=H.shape[0]))
        order=np.argsort(E);E=E[order];V=V[:,order]
        residual=np.max(np.linalg.norm(H@V-V*E,axis=0))
        T=b.translation(1 if p==1 else 2);groups=[];i=0
        while i<len(E):
            j=i+1
            while j<len(E) and abs(E[j]-E[i])<1e-8:j+=1
            W=V[:,i:j];small=W.conj().T@(T@W)
            leak=float(np.linalg.norm(T@W-W@small))
            groups.append(dict(energy=float(E[i]),multiplicity=j-i,
                translation_angle_over_pi=sorted((np.angle(np.linalg.eigvals(small))/np.pi).tolist()),translation_leakage=leak))
            i=j
        blocks.append(dict(parity=p,flip=q,physical_multiplicity=2 if p==-1 else 1,
            translation_steps=1 if p==1 else 2,dimension=H.shape[0],max_residual=float(residual),levels=groups))
        if ground is None or E[0]<ground[0]:ground=(float(E[0]),b,V[:,0].copy(),H,E.copy(),V.copy(),groups)
    E0,b,v,H,E,V,groups=ground;obs={};poles={}
    for name,sign in (('density_0',np.ones(L)),('density_pi',(-1.)**np.arange(L))):
        # Z switches flip parity; same representative index with no extra factor.
        Ov=(b.z@sign)*v
        obs[name]=float(np.vdot(Ov,Ov).real/L**2)
    for name,ops in (('XY',b.A),('YX',b.B),('ZZ',b.D),('energy',[2*a+g*d for a,d in zip(b.A,b.D)])):
        mean=sum(np.vdot(v,o@v).real for o in ops)/L
        obs[name+'_mean']=float(mean)
        for angle,label in ((np.pi,'pi'),(2*np.pi/L,'2pi_over_L')):
            Ov=sum(np.exp(1j*angle*j)*(o@v) for j,o in enumerate(ops))
            norm=float(np.vdot(Ov,Ov).real-abs(np.vdot(v,Ov))**2)
            obs[name+'_'+label+'_connected']=norm/L**2
            # A weighted mean excitation energy, not an inferred gap/velocity.
            w=Ov-v*np.vdot(v,Ov)
            obs[name+'_'+label+'_mean_excitation']=float(np.vdot(w,H@w).real/norm-E0) if norm>1e-15 else None
            overlaps=abs(V.conj().T@w)**2;i=0;visible=[]
            for group in groups:
                j=i+group['multiplicity'];weight=float(sum(overlaps[i:j]))
                if group['energy']-E0>1e-8 and weight>1e-10:
                    visible.append(dict(gap=group['energy']-E0,weight_over_L2=weight/L**2,
                        fraction_of_total=weight/norm,translation_leakage=group['translation_leakage']))
                i=j
            poles[name+'_'+label]=visible
    return dict(L=L,g=g,seconds=time.perf_counter()-start,ground_energy=E0,
        ground_sector=[b.parity,b.flip],observables=obs,retained_neutral_poles=poles,sectors=blocks)


if __name__=='__main__':
    # Refuse an accidental overlap with another invocation of this diagnostic.
    lock=open(Path(tempfile.gettempdir())/'cd_two_term_phase.lock','w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    parser=argparse.ArgumentParser();parser.add_argument('--lengths',nargs='+',type=int,default=[8,12,16,20]);parser.add_argument('--couplings',nargs='+',type=float,default=[.125,.5,1,4]);parser.add_argument('--levels',type=int,default=16);parser.add_argument('--output',default='spectra.json');args=parser.parse_args()
    assert os.environ.get('OPENBLAS_NUM_THREADS')=='1','Set OPENBLAS_NUM_THREADS=1'
    output=Path(__file__).parent/args.output
    payload=dict(method=f'scipy eigsh, {args.levels} levels in each of three symmetry blocks; fourth exactly paired',blas_threads=1,records=[])
    for L in args.lengths:
        for g in args.couplings:
            r=solve(L,g,k=args.levels);payload['records'].append(r)
            output.write_text(json.dumps(payload,indent=2)+'\n')
            print(L,g,round(r['seconds'],3),r['ground_energy'],flush=True)
