import time, numpy as np
from rt_bond import P
from rt_cell import second_order, one_particle_T, hs_quadratic, herm_part_analysis
for n, td in ((1, 1e-3), (2, 0.0), (2, 1e-3), (3, 0.0), (3, 1e-5)):
    t0 = time.time()
    r = second_order(n, P(td=td))
    T, muA, muB = one_particle_T(r['H2'], n)
    e0, h, res, ncomm = hs_quadratic(r['H2'], n)
    ef, best, asym = herm_part_analysis(T)
    sc = abs(T).max()
    print(n, td, 'dim', r['dim_cell'], 'deg', r['deg'], 'gap %.4f' % r['gap'], 'spread', r['code_spread'], 'orth', r['code_orth'], 'imp', r['code_impurity'],
          'first', r['first_order'], 'leak', r['V_leaves_exc'], '|T| %.4e' % sc, 'Rs', ef, 'crawford', best, 'asym %.3f' % asym,
          'mu/T', abs(muA+muB).max()/sc, 'Re mu/T', abs((muA+muB).real).max()/sc, 'NG/T', res/abs(h[:n,n:]).max(), 'Ncomm', ncomm, 'e0', e0.real, '%.1fs' % (time.time()-t0))
