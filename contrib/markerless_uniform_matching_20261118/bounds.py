"""Explicit volume-independent local-normal-form majorants.

These conservative constants certify a finite-order statement, not convergence
of an infinite Schrieffer-Wolff expansion or an interacting CFT.
"""
from math import exp, factorial, sqrt, log


def normal_form_majorant(n, D, chi, order=4, kappa=0., gap=.5):
    c = [D+8., 2*sqrt(2), n*abs(chi)]+[0.]*(order-2)
    generators = []
    for r in range(1, order+1):
        s = c[r]/gap
        generators.append(s)
        new = c.copy()
        for k in range(r+1, order+1):
            new[k] = sum(c[k-l*r]*(2*(k+2)*s)**l/factorial(l)
                         for l in range(k//r+1))
        c = new
    R = order+1
    weight = exp(kappa*R)
    radius = min([1.]+[(1/(4*R*weight*order*s))**(1/r)
                       for r, s in enumerate(generators, 1) if s > 0])
    total = sum(s*radius**r for r, s in enumerate(generators, 1))
    assert 2*R*weight*total <= .5*(1+1e-12)
    initial = D+8.+radius*2*sqrt(2)+radius**2*n*abs(chi)
    analytic_bound = 4*weight*initial
    block_ground_radius = min([1.]+[(gap/(order*2**(r+4)*c[r]))**(1/r)
                                    for r in range(1,order+1) if c[r]>0])
    return dict(order=order, kappa=kappa, gap=gap, generator_bounds=generators,
                coefficient_bounds=c, analytic_radius=radius,
                block_ground_radius=block_ground_radius,
                analytic_bound=analytic_bound,
                remainder_constant=2*analytic_bound/radius**(order+1),
                projected_even_remainder_constant=2*analytic_bound/radius**(order+2))


def explicit_fourth_bound(n, zeta, chi, kappa=0., clock_gap=.5, pair_gap=3.):
    two = 4*zeta/(pair_gap**2*clock_gap)+4*zeta**2/pair_gap**3
    three = 2*two
    mixed = 8*n*abs(chi)*zeta/pair_gap**2+8*n*abs(chi)*sqrt(zeta)/(clock_gap*pair_gap)
    field = n*n*chi*chi/clock_gap
    strength = 2*exp(2*kappa)*(two+mixed)+3*exp(3*kappa)*three+exp(kappa)*field
    return dict(two_cell_port_bound=two, three_cell_connected_port_bound=three,
                mixed_readout_bound=mixed, single_cell_field_squared_bound=field,
                total_local_strength_bound=strength)


def car_weight(n):
    return 1+n*log(8*(2*n)**2)
