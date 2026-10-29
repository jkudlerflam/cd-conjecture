"""Small checks of C65's regulator results. No all-depth RG calculation."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit


def bubble_checks():
    A, cutoff, T = 1.3, 2.0, 0.7
    n = lambda k: expit(-A*k/T)
    base = (1-2*n(cutoff))/(2*math.pi*A)
    cusp = -n(cutoff)*(1-n(cutoff))/(2*math.pi*T)
    rows = []
    errors = []
    for p in (0.08, 0.02, 0.005, 0.00125):
        direct = quad(lambda k: (n(k)-n(k+p))/(2*math.pi*A*p),
                      -cutoff, cutoff-p, epsabs=1e-12)[0]
        edge = (1-2*quad(lambda s: n(cutoff-p*s), 0, 1,
                        epsabs=1e-13)[0])/(2*math.pi*A)
        negative = quad(lambda k: (n(k)-n(k-p))/(-2*math.pi*A*p),
                        -cutoff+p, cutoff, epsabs=1e-12)[0]
        slope = (edge-base)/p
        assert abs(direct-edge) < 3e-13
        assert abs(negative-edge) < 3e-13
        errors.append(abs(slope-cusp))
        rows.append(dict(p=p, direct=direct, boundary=edge,
                         right_slope=slope, cusp_coefficient=cusp))
    assert errors[-1] < errors[0]/30
    assert abs(rows[-1]['right_slope']/cusp-1) < 0.002

    def r(k):
        s = abs(k)/cutoff
        if s <= 1:
            return 1.0
        if s >= 2:
            return 0.0
        t = s-1
        a, b = math.exp(-1/t), math.exp(-1/(1-t))
        return b/(a+b)

    # Integrate the divided difference as an average of n' to avoid cancellation.
    nodes, weights = np.polynomial.legendre.leggauss(12)
    nodes, weights = (nodes+1)/2, weights/2
    def smooth(p):
        def integrand(k):
            occupied = expit(-A*(k+nodes*p)/T)
            diff = np.dot(weights, occupied*(1-occupied))/T
            return r(k)**2*r(k+p)**2*diff/(2*math.pi)
        return quad(integrand, -2*cutoff, 2*cutoff,
                    points=[-cutoff, cutoff, -cutoff-p, cutoff-p],
                    epsabs=2e-12, epsrel=2e-12)[0]
    sm0 = smooth(0)
    smrows = []
    for p in (0.04, 0.02, 0.01, 0.005):
        sp, sn = smooth(p), smooth(-p)
        assert abs(sp-sn) < 1e-11
        smrows.append(dict(p=p, one_sided_slope=(sp-sm0)/p,
                           symmetric_derivative=(sp-sn)/(2*p)))
    assert abs(smrows[-1]['one_sided_slope']) < abs(smrows[0]['one_sided_slope'])/6
    assert max(abs(x['symmetric_derivative']) for x in smrows) < 1e-9
    return dict(sharp=rows, sharp_slope_errors=errors, smooth=smrows)


def cusp_moment_check():
    # Fourier inverse of |p| exp(-|p|), with convention dp/(2*pi).
    kernel = lambda x: (1-x*x)/(math.pi*(1+x*x)**2)
    def exact(R):
        return (math.log(1+R*R)+2/(1+R*R)-math.log(2)-1)/math.pi
    rows = []
    previous = None
    increments = []
    for R in (4., 8., 16., 32., 64., 128.):
        moment = 2*quad(lambda x: x*abs(kernel(x)), 1, R, epsabs=1e-12)[0]
        assert abs(moment-exact(R)) < 1e-12
        rows.append(dict(R=R, truncated_first_moment=moment))
        if previous is not None:
            increments.append(moment-previous)
        previous = moment
    expected_increment = 2*math.log(2)/math.pi
    assert abs(increments[-1]-expected_increment) < 0.001
    return dict(rows=rows, doubling_increment_limit=expected_increment,
                last_doubling_increment=increments[-1])


def annihilator(modes, j):
    a = np.zeros((1 << modes, 1 << modes))
    for state in range(1 << modes):
        if (state >> j) & 1:
            sign = (-1)**((state & ((1 << j)-1)).bit_count())
            a[state ^ (1 << j), state] = sign
    return a


def attenuation_checks():
    modes = 4
    dim = 1 << modes
    d = [annihilator(modes, j) for j in range(modes)]
    occ = [a.T@a for a in d]
    A, B = 1.3, 0.2
    H0 = A/2*sum(occ)
    Qp, Qm = occ[1]-occ[0], occ[3]-occ[2]
    pair = d[0]@d[1]@d[2]@d[3]
    V = B*(pair+pair.T+Qp@Qm)
    r = np.array([0.85, 0.65, 0.75, 0.50])
    # Construct the second-quantized isometry directly in the doubled CAR basis.
    big_d = [annihilator(2*modes, j) for j in range(2*modes)]
    creators = [r[j]*big_d[j].T+math.sqrt(1-r[j]**2)*big_d[j+modes].T
                for j in range(modes)]
    J = np.zeros((1 << (2*modes), dim))
    for state in range(dim):
        vec = np.zeros(1 << (2*modes))
        vec[0] = 1
        for j in reversed(range(modes)):
            if (state >> j) & 1:
                vec = creators[j]@vec
        J[:, state] = vec
    assert np.max(abs(J.T@J-np.eye(dim))) < 2e-15
    # System modes occupy the low bits; even system operators are I_aux tensor O.
    phiV = J.T@np.kron(np.eye(dim), V)@J
    phiH0 = J.T@np.kron(np.eye(dim), H0)@J
    targetV = B*(np.prod(r)*(pair+pair.T)
                  +(r[1]**2*occ[1]-r[0]**2*occ[0])
                   @(r[3]**2*occ[3]-r[2]**2*occ[2]))
    targetH0 = A/2*sum(r[j]**2*occ[j] for j in range(modes))
    residual = max(np.max(abs(phiV-targetV)), np.max(abs(phiH0-targetH0)))
    assert residual < 2e-15
    bound = (abs(B)/A)*H0+abs(B)*np.eye(dim)
    lower_margins = [float(np.linalg.eigvalsh(bound+sign*W)[0])
                     for W in (V, phiV) for sign in (-1, 1)]
    assert min(lower_margins) > -2e-14
    assert np.linalg.eigvalsh(H0-phiH0)[0] > -2e-14
    return dict(system_dimension=dim, doubled_dimension=dim*dim,
                attenuation_identity_residual=float(residual),
                two_sided_bound_minimum_eigenvalues=lower_margins)


if __name__ == '__main__':
    results = dict(scope='Regulator coefficient, moment and finite CAR checks. TP remains unproved.',
                   bubble=bubble_checks(), cusp_moment=cusp_moment_check(),
                   attenuation=attenuation_checks())
    Path(__file__).with_name('checks.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(dict(status='passed',
                         sharp_cusp=results['bubble']['sharp'][-1]['cusp_coefficient'],
                         sharp_slope_error=results['bubble']['sharp_slope_errors'][-1],
                         smooth_last_slope=results['bubble']['smooth'][-1]['one_sided_slope'],
                         attenuation=results['attenuation']), indent=2))
