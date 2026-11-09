"""Full-spin MPS number fluctuations and explicitly defined density smoothness.
Number is an observable, not an imposed symmetry; these measurements do not
establish quasiparticle density or criticality.
"""
import argparse,json
from pathlib import Path
import numpy as np

def number_moments(M):
    n=np.diag([0.,1.]);I=np.eye(2)
    e0=np.ones((1,1),complex);e1=np.zeros_like(e0);e2=np.zeros_like(e0)
    def tr(e,A,o):return np.einsum('ab,asc,st,btd->cd',e,A.conj(),o,A,optimize=True)
    for A in M:
        e0,e1,e2=tr(e0,A,I),tr(e1,A,I)+tr(e0,A,n),tr(e2,A,I)+2*tr(e1,A,n)+tr(e0,A,n)
    norm=e0.item().real;mean=e1.item().real/norm;second=e2.item().real/norm
    return dict(norm=norm,N=mean,N2=second,Nvariance=second-mean**2)

def main():
    a=argparse.ArgumentParser();a.add_argument('inputs',nargs='+');a.add_argument('--output',required=True);args=a.parse_args();rows=[]
    for filename in args.inputs:
        d=json.loads(Path(filename).read_text());L=d['L'];z=np.load(d['mps_checkpoint']);M=[z[f'M{j}'] for j in range(L)]
        r=dict(file=filename,L=L,parity=d['parity'],**number_moments(M));rho=(1-np.array(d['Z']))/2
        # Smooth microscopic two-site Friedel oscillations with a four-site box,
        # then inspect only windows wholly inside the central half of the chain.
        central=rho[L//4:3*L//4];smooth=np.convolve(central,np.ones(4)/4,mode='valid')
        r['density_central_half_mean']=float(central.mean());r['density_four_site_smoothed_range']=float(np.ptp(smooth));r['density_four_site_relative_range']=float(np.ptp(smooth)/smooth.mean())
        r['N_from_profiles_difference']=r['N']-d['N']
        if L<=12:
            psi=M[0]
            for A in M[1:]:psi=np.tensordot(psi,A,axes=([-1],[0]))
            p=abs(psi.reshape(-1))**2;ns=np.array([k.bit_count() for k in range(2**L)]);p/=p.sum();v=float(p@ns**2-(p@ns)**2)
            r['dense_number_variance_error']=r['Nvariance']-v
        rows.append(r)
    Path(args.output).write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
