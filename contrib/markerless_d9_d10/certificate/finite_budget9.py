"""Finite conversion strength and local two-cell port remainder for the markerless 9-letter members,
following C84 PROOF.md (10) and (17) with the finite-K term deleted (the grammar is exactly invariant, so
there is no finite-K error).  Exact Fractions throughout.
Inherited, not re-derived here: the finite-g tail bound (10) and the two-cell remainder 2^50 lambda^4.
Their hypotheses that depend on the bond are checked below or in the README: source norm sqrt(2) alpha^2,
pair conversion norm 2g, free pair gap >= 3, and a vacuum band isolated by the clock gap > 1/32.
Usage: finite_budget9.py n  [directory holding derivative_n{n}.npz and root_certificate_exact_n{n}.json]"""
import json, math, sys, os
from fractions import Fraction as Q
import numpy as np
from rational_bands import integer_ldl

def main(n, d='.'):
    a = np.load(os.path.join(d, f'derivative_n{n}.npz'), allow_pickle=True)
    R, ds = integer_ldl(a['F_numerator'] + a['F_numerator'].T)
    bits = int(a['bits']); powers = [(int(x).bit_length()-1)//2 for x in ds]
    trace = sum((Q((1 << (2*bits))*sum(int(x)**2 for x in R[i]), 1 << (2*powers[i])) for i in range(n)), Q(0))
    det = int(a['detC'])
    # alpha^4 = (2/(n+2))^2 sin^4(pi/(n+2)) >= 64/(n+2)^6, from sin x >= 2x/pi
    apow = 1
    while Q(64, (n+2)**6) <= Q(1, 1 << apow): apow += 1
    assert Q(64, (n+2)**6) > Q(1, 1 << apow)
    root = json.load(open(os.path.join(d, f'root_certificate_exact_n{n}.json')))
    budget = Q(1, 1 << root['perturbation_budget_power'])
    # clock gap of the selected cell must exceed 1/32 for the contour argument
    clock_gap = 2*math.cos(math.pi/(n+2)) - 2*math.cos(2*math.pi/(n+2))
    assert clock_gap > 1/32 + 1e-6
    def errs(gp, lp):
        unit = Q(1, (1 << apow)*det**2*(1 << (2*gp)))            # lower bound on E_unit = alpha^4 beta^2 g^2
        conv = 4*n*det**2*trace/Q(1 << (2*gp))                   # error_g = 4 n beta^-2 g^2 tau
        lam = 16*n*trace*Q(1 << 50, 1 << (2*lp))/unit            # error_lambda = 16 n tau 2^50 lambda^2 / E_unit
        return unit, conv, lam
    gp = 1
    while errs(gp, 10**6)[1] >= budget/16: gp += 1
    lp = gp
    while errs(gp, lp)[2] >= budget/16: lp += 1
    unit, conv, lam = errs(gp, lp)
    assert conv + lam < budget/8
    lg = lambda x: math.log10(x.numerator) - math.log10(x.denominator)
    out = dict(n=n, alphabet=9, cell_length=n+1, g_power=-gp, lambda_power=-lp, alpha4_lower_power=-apow,
               Eunit_lower_log10=lg(unit), metric_trace_log10=lg(trace), conversion_error_log10=lg(conv),
               two_cell_lambda_error_log10=lg(lam), budget_power=-root['perturbation_budget_power'],
               clock_gap=clock_gap, sum_below_one_eighth_budget=True, finite_K_term='absent (exact invariance)')
    json.dump(out, open(os.path.join(d, f'finite_budget_n{n}.json'), 'w'), indent=2)
    print(json.dumps(out), flush=True)
if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else '.')
