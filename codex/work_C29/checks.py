"""Independent algebra, heat kernel, finite-g, and Bloch checks."""
import json
from itertools import combinations
from fractions import Fraction as Q
import numpy as np
from scipy.linalg import eigh
from model import *


def rational_weights(d):
    z=Q(2*d);a2=Q(2);I20=2/z**3
    def I2c(b2):return 2*z*(z*z+3*b2)/(z*z-b2)**3
    I2cc=(I20+I2c(4*a2))/2
    I1s_a=2*z/(z*z-a2)**2;I1cs_a=2*z/(z*z-4*a2)**2;I0ss_a2=2/(z*(z*z-4*a2))
    F00=(4*I20+4*I2c(a2)+I2cc+4*I1s_a+2*I1cs_a-15*I0ss_a2)/64
    F11=(I2cc-2*I1cs_a+I0ss_a2)/16
    return F00,F11


def full_source(n,g):
    H,G,ports,alpha,E,q=cell(n,g);e,U=eigh(H)
    L,R=ports;W=np.zeros((len(H)**2,4**n),complex)
    for r,(_,_,phi,w) in enumerate(parameters(n)):
        aa=[U.conj().T@F@G for F in R[r]];bb=[U.conj().T@F@G for F in L[r]]
        W+=w*(np.kron(aa[0],bb[0])+np.exp(1j*phi)*np.kron(aa[1],bb[2])+np.exp(-1j*phi)*np.kron(aa[2],bb[1]))
    den=(e[:,None]+e).reshape(-1);iv=np.zeros_like(den);iv[den>1e-9]=1/den[den>1e-9]
    K=-W.conj().T@(iv[:,None]*W)
    return float(np.linalg.norm(K-alpha**4*combined_bond(n,g),2))


def main():
    out={};rng=np.random.default_rng(2901);rows=[]
    for n in (1,2,3):
        ph=rng.uniform(-1,1,size=(n,3));phi=.27
        K=bond(n,4.,ph,phi,coefficient=True);one,res=quadratic_fit(K,2*n)
        T=formula(n,4.,ph,phi)
        rows.append(dict(n=n,excited_cell_dimension=3*(n+1)*2**n,formula_error=float(np.linalg.norm(one[:n,n:]-T,2)),quadratic_residual=res,field_norm=float(np.linalg.norm(one[:n,:n])+np.linalg.norm(one[n:,n:]))))
    out['nonuniform_flux_checks']=rows
    rats=[]
    for d in (3,5):
        a,b=rational_weights(d);F=heat_weights(2,d)
        assert abs(float(a)-F[0,0])<1e-16 and abs(float(b)-F[1,1])<1e-16
        rats.append(dict(d=d,F00=str(a),F11=str(b),ratio=str(a/b)))
    a,b=rational_weights(3);c,d=rational_weights(5)
    out['rational_weights']=rats;out['phase_determinant']=str(a*d-b*c)
    out['complete_source_errors']=[dict(n=n,physical_legal_dimension=(7*(n+1)*2**n)**2,error=full_source(n,.2)) for n in (1,2)]
    out['finite_g']=[]
    for g in (.2,.1,.05):
        K=combined_bond(2,g);o,res=quadratic_fit(K,4);T=o[:2,2:]
        out['finite_g'].append(dict(g=g,nonquadratic_norm=res,normal_relative=float(np.linalg.norm(T@T.conj().T-T.conj().T@T)/np.linalg.norm(T)**2),field_norm=float(np.linalg.norm(o[:2,:2]+o[2:,2:]))))
    T,mu,W=target_data(2);theta=np.angle(np.diag(T));r=np.abs(np.diag(T));kf=np.arccos(-mu/(2*r));points=np.stack([-theta-kf,-theta+kf],axis=1)
    points=(points+np.pi)%(2*np.pi)-np.pi
    out['n2_target']=dict(T_real=T.real.tolist(),T_imag=T.imag.tolist(),mu=mu.tolist(),theta=theta.tolist(),kf=kf.tolist(),fermi_points=points.tolist(),velocities=(2*r*np.sin(kf)).tolist(),normal_error=float(np.linalg.norm(T@T.conj().T-T.conj().T@T)),phase_difference_mod_pi=float((theta[1]-theta[0]+np.pi/2)%np.pi-np.pi/2),pair_transfer_mismatch=float((2*(theta[1]-theta[0])+np.pi)%(2*np.pi)-np.pi))
    T3=sum(w*w*formula(3,d,ph,phi) for d,ph,phi,w in parameters(3))
    flat=points.ravel(); pairs=list(combinations(range(4),2)); sums=[flat[i]+flat[j] for i,j in pairs]
    out['n2_target']['bilinear_mismatch']=float(min(abs((flat[i]-flat[j]+np.pi)%(2*np.pi)-np.pi) for i,j in combinations(range(4),2)))
    out['n2_target']['distinct_pair_sum_mismatch']=float(min(abs((sums[i]-sums[j]+np.pi)%(2*np.pi)-np.pi) for i,j in combinations(range(6),2)))
    out['n3_failure']=dict(T_real=T3.real.tolist(),T_imag=T3.imag.tolist(),normal_norm=float(np.linalg.norm(T3@T3.conj().T-T3.conj().T@T3)),normal_relative=float(np.linalg.norm(T3@T3.conj().T-T3.conj().T@T3)/np.linalg.norm(T3)**2))
    print(json.dumps(out,indent=2))
    with open('checks_results.json','w') as f:json.dump(out,f,indent=2)

if __name__=='__main__':main()
