"""Small matrix checks of comparisons with no imposed bond symmetry."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import logsumexp
from scipy.integrate import quad
from check_rings import logz_xx

ROOT = Path(__file__).resolve().parent

def ring(h, length):
    dim = 2**length
    result = np.zeros((dim, dim), dtype=complex)
    for state in range(dim):
        for i in range(length):
            j = (i+1) % length
            ai, aj = (state >> i) & 1, (state >> j) & 1
            source = 2*ai+aj
            for bi in range(2):
                for bj in range(2):
                    target = state + (bi-ai)*(2**i)+(bj-aj)*(2**j)
                    result[target,state] += h[2*bi+bj,source]
    return result

def main():
    eye = np.eye(2)
    x = np.array([[0,1],[1,0]])
    y = np.array([[0,-1j],[1j,0]])
    z = np.diag([1,-1])
    h = (0.31*np.kron(x,x)+0.23*np.kron(y,z)+0.17*np.kron(z,eye)
         -0.11*np.kron(eye,y)+0.07*np.kron(x,y))
    evh = eigvalsh(h)
    h -= (evh[0]+evh[-1])/2*np.eye(4)
    m = (evh[-1]-evh[0])/2
    o = 0.4*x+0.3*y-0.2*z
    hg = h+np.kron(o,eye)-np.kron(eye,o)+0.37*np.eye(4)
    e3,e6 = eigvalsh(ring(h,3)),eigvalsh(ring(h,6))
    eg3,eg6 = eigvalsh(ring(hg,3)),eigvalsh(ring(hg,6))
    assert np.max(np.abs(eg3-e3-3*0.37)) < 1e-12
    assert np.max(np.abs(eg6-e6-6*0.37)) < 1e-12
    rows=[]
    for beta in [0.0,0.1,0.5,1.0,4.0,16.0]:
        r=2*logsumexp(-beta*e3)-logsumexp(-beta*e6)
        rg=2*logsumexp(-beta*eg3)-logsumexp(-beta*eg6)
        assert abs(r-rg) < 1e-10
        assert abs(r)<=4*beta*m+1e-10
        gt=3*logsumexp(-2*beta*eigvalsh(h))
        assert logsumexp(-beta*e6)<=gt+1e-10
        rows.append({'beta':beta,'R_L3':float(r),'gluing_ceiling':float(4*beta*m),
                     'GT_slack_L6':float(gt-logsumexp(-beta*e6))})
    xx=[]
    for beta in [4.0,8.0,16.0,32.0]:
        pressure,err=quad(lambda k:np.logaddexp(0,-beta*np.cos(k)),0,2*np.pi,
                          epsabs=1e-11,epsrel=1e-11)
        pressure/=2*np.pi
        xi=logz_xx(4,beta)-4*pressure
        lower=beta*(np.sqrt(2)-4/np.pi)-4*np.log(2)
        assert xi>=lower-1e-10
        xx.append({'beta':beta,'Xi_L4':xi,'log_d2':float(np.log(4)),
                   'rigorous_lower_bound':float(lower),'quadrature_error_estimate':err/(2*np.pi)})
    output={'generic_bond_raw_half_spread':float(m),'generic_checks':rows,
             'XX_pressure_correction':xx,'XX_delta4':float(np.sqrt(2)-4/np.pi)}
    twisted = np.zeros((4,4),dtype=complex)
    twisted[1,2] = np.exp(1j*np.pi/4)/2
    twisted[2,1] = np.conj(twisted[1,2])
    et4,et8 = eigvalsh(ring(twisted,4)),eigvalsh(ring(twisted,8))
    assert abs(et4[0]+1) < 1e-12
    assert abs(et8[0]+1/np.sin(np.pi/8)) < 1e-12
    rt = 2*logsumexp(-16*et4)-logsumexp(-16*et8)
    assert rt < 0
    output['twisted_XX'] = {'theta':float(np.pi/4),'E4':float(et4[0]),
                            'E8':float(et8[0]),'R_L4_beta16':float(rt),
                            'doubling_energy_defect':float(et8[0]-2*et4[0])}
    slow = []
    length, tau = 2048, 1.0
    for epsilon in [0.5,0.1,0.02,0.005]:
        v1,v2 = 1/(1+epsilon),epsilon/(1+epsilon)
        def stack_lz(n,b):
            return logz_xx(n,b*v1)+logz_xx(n,b*v2)
        b=tau*length
        r=2*stack_lz(length,b)-stack_lz(2*length,b)
        q=4*stack_lz(length,b)-stack_lz(2*length,2*b)
        slow.append({'epsilon':epsilon,'L':length,'tau':tau,'d':4,'M':0.5,
                     'A':1.0,'B':1/v1+1/v2,'R':r,'Q':q})
    output['two_velocity_XX'] = slow
    (ROOT/'comparisons.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__': main()
