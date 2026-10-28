"""Finite checks of analytic normalization bounds for the C8 candidate."""
import json
import numpy as np
from periodic_clean_bond import PeriodicBond


def case(q):
    H=PeriodicBond(q)
    m,n,b=H.m,H.k,H.b;T=11*m*m+105*m
    lower=np.array([(n-3)*(2*n+1)/(T*(2*n+5)),
                    (n-4)*(2*n+1)/(T*(2*n+5))])
    assert min(H.alpha.flat)>=lower[0] and min(H.beta.flat)>=lower[1]
    assert min(lower)>=1/(80*m)
    droplet=[]
    for p in (1,2):
        hs,chi=H.holes[p-1],H.chis[p-1]
        expected=(1-2*p*H.w/H.U)/2
        left=hs.index(tuple(range(p)))
        right=hs.index(tuple(range(n,n+p)))
        actual=min(chi[left]**2,chi[right]**2)
        assert actual>=expected-1e-12 and expected>1/4
        droplet.append(dict(vacancies=p,measured_endpoint_droplet_weight=float(actual),
                            lower_bound=expected))
    assert min(H.r+H.s2)>1/4
    norms=[float(np.linalg.norm(mat,2)) for mat in (H.ML,H.MR,H.MD)]
    bound=8*np.pi**2*b*b
    assert max(norms)<=bound and H.marker_norm>1/64
    actual_B=(abs(H.kappa)+2*H.marker_u+
              4*norms[0]*norms[1]*float(np.sum(1/(H.alpha[0]*H.alpha[1])))+
              norms[2]*float(np.max(abs(H.mu)/(2*H.beta[0]))))
    polynomial_B=129+2*m+1638400*np.pi**4*b**4*m**3+640*np.pi**2*b*b*m
    assert actual_B<=polynomial_B
    return dict(q=q,m=m,port_visibility_min=float(min(H.alpha.flat)),
                density_visibility_min=float(min(H.beta.flat)),
                visibility_lower_bound=1/(80*m),droplets=droplet,
                compensator_norms=norms,compensator_upper_bound=bound,
                marker_visibility=H.marker_norm,marker_lower_bound=1/64,
                norm_triangle_bound=actual_B,polynomial_norm_bound=float(polynomial_B))


if __name__=='__main__':
    out=[case(q) for q in (4,8)]
    with open('normalization_bounds_checks.json','w') as f:
        json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))
