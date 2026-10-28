"""Canonical pressure and marker-selection checks on scalar physical rings."""
import json
from pathlib import Path

import numpy as np

from canonical_length_checks import atom, low, physical_matrix


def pressure(alpha, J0=1.):
    rho = np.arccos(np.clip(-alpha/(2*J0), -1., 1.))/np.pi
    return -2*J0*np.sin(np.pi*rho)/np.pi-alpha*rho


def cost(rho, J0=1.):
    rho = np.asarray(rho)
    return np.where(rho < 0, 2*J0/np.pi-2*J0*rho,
                    np.where(rho > 1, 2*J0/np.pi+2*J0*(rho-1),
                             2*J0/np.pi*(1-np.sin(np.pi*rho))))


def main():
    b = 3; J0 = 1.; kappa = 2/np.pi; epsilon = .0005
    Delta = 2.; gamma = 1.; C = 1.; J = kappa+4*J0; beta = J0
    delta = .1; eta = .2; a = b+.5
    gd = 2*J0/np.pi*(1-np.cos(np.pi*delta))
    sd = 2*J0*np.sin(np.pi*delta)
    D0 = 2*J0*C*Delta+2*J+4*J0
    assert epsilon*D0 <= Delta/2
    assert epsilon <= gd*Delta/(4*J**2)
    assert epsilon*J <= min(eta/(4*C), eta*gamma/16)
    rho = np.linspace(-4, 6, 10001)
    alpha = -2*J0*np.cos(np.pi*np.clip(rho, 0, 1))
    legendre_error = float(np.max(np.abs(kappa+pressure(alpha)+alpha*rho-cost(rho))))
    assert legendre_error < 4e-15
    exterior = np.abs(rho-.5) >= delta
    tangent_margin = cost(rho)-(gd+sd*(np.abs(rho-.5)-delta))
    assert tangent_margin[exterior].min() > -2e-15
    cases = 0; minimum_margin = np.inf
    for L in range(3, 2001):
        markers = np.arange(1, L//2+1)
        densities = L/markers-b
        # Exact integer exterior test for b=3 and delta=1/10, including equality.
        use = np.abs(10*L-35*markers) >= markers
        margins = markers[use]*(cost(densities[use])-gd/2)-gd*L/(2*(a+delta))
        if len(margins): minimum_margin = min(minimum_margin, float(margins.min()))
        assert np.all(margins >= -3e-12)
        cases += len(margins)
    finite = []
    selected = []
    for L in (10, 14, 17, 18, 21):
        energies = [(atom(L, b)+epsilon*(kappa-2*J0), 1)]
        for M in range(2, L//2+1):
            basis, gaps, H = physical_matrix(M, L, epsilon*J0, b)
            H.setdiag(H.diagonal()+epsilon*kappa*M)
            e, U, residual = low(H)
            E = float(e[0]); N = L-b*M; density = N/M
            lower = epsilon*M*float(cost(density))-2*epsilon*beta-2*M*epsilon**2*J**2/Delta
            assert E >= lower-3e-11
            tilt_checks = []
            for al in (-2., -1., 0., 1., 2.):
                mu = epsilon*al
                gmu = Delta*(1-C*abs(mu))-2*(epsilon*J+abs(mu))
                r = np.arange(M+1)
                eref = epsilon*kappa*M+np.min(-2*epsilon*J0*np.sin(np.pi*r/M)/np.sin(np.pi/M)-mu*r)
                bound = eref-M*epsilon**2*J**2/gmu
                assert E-mu*N >= bound-3e-11
                bare_A = np.array([sum(atom(ell, b) for ell in row) for row in gaps])
                p1 = (gaps == b+1).sum(axis=1)
                Amu = bare_A-mu*(N-p1)
                assert np.all(Amu >= (1-C*abs(mu))*bare_A-1e-12)
                tilt_checks.append(dict(alpha=al, actual=E-mu*N,
                                        atomic_bound=float(bound), gap=gmu))
            finite.append(dict(length=L, markers=M, dimension=len(basis),
                               density=density, ground=E, bulk_bound=lower,
                               residual=residual, tilt=tilt_checks))
            energies.append((E, M))
        E, M = min(energies)
        Mstar = round(L/a); rstar = L-b*Mstar
        BL = 2*beta+np.pi*J0*a*a/(4*Mstar)
        if 0 <= rstar <= Mstar:
            assert E <= epsilon*BL+3e-11
        selected.append(dict(length=L, markers=M, density=L/M-b, ground=E,
                             trial_upper=epsilon*BL,
                             certified_size_condition=bool(gd*L/(2*(a+delta)) > 2*beta+BL)))
    out = dict(status='passed, scalar physical controls only; no central-charge inference',
               parameters=dict(b=b, J0=J0, kappa=kappa, epsilon=epsilon, Delta=Delta,
                               gamma=gamma, C=C, full_norm=J, beta=beta, delta=delta,
                               density_selection_coupling_bound=gd*Delta/(4*J**2)),
               algebra=dict(legendre_error=legendre_error, marker_cases=cases,
                            minimum_marker_margin=minimum_margin),
               sectors=finite, selected=selected)
    Path('marker_density_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(status='passed', sectors=len(finite),
                         largest_dimension=max(row['dimension'] for row in finite),
                         largest_residual=max(row['residual'] for row in finite),
                         algebra=out['algebra'], selected=selected), indent=2))


if __name__ == '__main__':
    main()
