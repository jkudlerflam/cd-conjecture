"""Read-only exact-environment energy/tangent audit of accepted infinite MPS."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import hashlib,json,pickle,time,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import tenpy
from tenpy.networks.uniform_mps import UniformMPS
from tenpy.networks.mpo import MPOTransferMatrix
from tenpy.algorithms.vumps import VUMPSEngine
import run as R
HERE=Path(__file__).resolve().parent
names=[f'candidate_infinite_L2_P0_s17_chi{c}_g_one' for c in [64,96,128]]
names+=['candidate_infinite_L2_P0_s17_chi192_g','ising_critical_infinite_L2_P0_s23_chi16_one','ising_critical_infinite_L2_P0_s23_chi24_a1']
rows=[]
for name in names:
    start=time.time();f=HERE/'results'/f'{name}.json';r=json.loads(f.read_text())
    assert r['status']=='done' and r['converged']
    ck=Path(r['checkpoint']);before=hashlib.sha256(ck.read_bytes()).hexdigest();assert before==r['checkpoint_sha256']
    row=dict(source=str(f.relative_to(HERE)),source_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),checkpoint_sha256=before)
    try:
        with ck.open('rb') as h:s=pickle.load(h)
        psi=s['psi'];u=UniformMPS.from_MPS(psi)
        m=(R.GaugeChain if r.get('gauge',False) else R.PauliChain)(dict(L=r['L'],bc_MPS='infinite',bc_x='periodic',couplings=r['couplings'],conserve=None))
        env,Es,_=MPOTransferMatrix.find_init_LP_RP(m.H_MPO,u,calc_E=True,tol=1e-12,maxiter=100000)
        left,right=VUMPSEngine.tangent_projector_test(SimpleNamespace(psi=u,model=m),env)
        Es=np.asarray(Es);assert max(abs(Es.imag))<1e-10
        row.update(status='done',environment_energies_shifted=Es.real.tolist(),
            energy_disagreement=float(np.max(np.abs(Es.real-r['E_direct_shifted']))),
            tangent_left=[float(x) for x in left],tangent_right=[float(x) for x in right],
            max_tangent=float(max(left+right)),one_site=r['config'].get('one_site',False),
            norm=float(np.linalg.norm(u.norm_test())))
    except Exception as e:
        row.update(status='failed',error=repr(e),traceback=traceback.format_exc())
    assert before==hashlib.sha256(ck.read_bytes()).hexdigest()
    row['seconds']=time.time()-start;rows.append(row)
    out=dict(tenpy_version=tenpy.__version__,complete=len(rows)==len(names),rows=rows)
    tmp=HERE/'tangent_audit.tmp';tmp.write_text(json.dumps(out,indent=2)+'\n');tmp.replace(HERE/'tangent_audit.json')
    print(json.dumps({k:v for k,v in row.items() if k!='traceback'}),flush=True)
