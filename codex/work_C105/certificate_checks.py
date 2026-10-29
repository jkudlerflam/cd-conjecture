"""Lightweight checks of the XX congruence certificate and SDP obstruction."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import gammaln, polygamma
from xx_halfline import profile


def site_factors(n):
    j = np.arange(n, dtype=float)
    r = 2*np.exp(2*(gammaln((j+2)/2)-gammaln((j+1)/2)))
    return np.sqrt(r*r[::-1]/n)


def run():
    certificate = []
    for n in [8, 16, 32, 64, 128, 256]:
        a = site_factors(n)
        x = np.arange(1, n, dtype=float)
        w = x*(n-x)/n
        eig = eigh_tridiagonal(np.zeros(n), w/2, eigvals_only=True)
        deficit = -w.sum()/np.pi - eig[eig<0].sum()
        dirichlet = np.sum(np.diff(np.r_[0, a, 0])**2)
        allowance = dirichlet/(2*np.pi)
        assert deficit <= allowance+1e-10
        assert np.max(abs(a[:-1]*a[1:]-w)) < 1e-10
        certificate.append(dict(n=n, defect=float(deficit),
                                allowance=float(allowance),
                                dirichlet=float(dirichlet),
                                factor_residual=float(np.max(abs(a[:-1]*a[1:]-w)))))
    sdp = []
    for r in range(2, 9):
        j = np.arange(r)
        k = j[:, None]-j[None, :]
        # Sign of hopping does not change these eigenvalues.
        corr = .5*np.sinc(k/2)
        nu = np.linalg.eigvalsh(corr)
        minimum = float(np.prod(np.minimum(nu, 1-nu)))
        t = .5*2**r*minimum
        sdp.append(dict(r=r, minimum_local_eigenvalue=minimum,
                        affine_step=t, feasible_energy=-float((1+t)/np.pi),
                        strict_improvement=float(t/np.pi)))
    catalan = (polygamma(1, .25)-polygamma(1, .75))/16
    expected = 4*catalan/np.pi**2-1/np.pi
    first = profile(1)[0]
    assert abs(first-expected)<1e-13
    old = json.loads(Path('small_checks_results.json').read_text())['rows']
    extension = []
    for model in ['XX', 'Heisenberg', 'XXZ_half', 'SU3']:
        v = [z for z in old if z['model']==model and z['shape']=='quadratic'][-1]
        a = site_factors(v['n'])
        D = np.sum(np.diff(np.r_[0, a, 0])**2)
        shift = dict(XX=0, Heisenberg=.25, XXZ_half=.125, SU3=0)[model]
        trial_allowance = -(v['e0']+shift)*D/2
        extension.append(dict(model=model, n=v['n'], defect=v['defect'],
                              centered_e0=v['e0']+shift, allowance=float(trial_allowance),
                              passes=bool(v['defect']<=trial_allowance+1e-10)))
    out = dict(certificate=certificate, fixed_range_sdp=sdp,
               first_bond_exact=float(expected), first_bond_quadrature=float(first),
               naive_interacting_extension=extension)
    Path('certificate_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    run()
