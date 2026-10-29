"""(4) Frozen pin selection (PROOF Section 6), referee checks.
 (a) Exact arithmetic for r_pin = 2 cos(pi/16) - B*^2/(K - 2B*) with B* = 2^16, K = 2^2560,
     and an independent norm estimate supporting ||B_cell|| <= 2^16 at the schedule (16).
 (b) Legal pinned cell words L_a0 ... L_a(p-2) L1 m are exact eigenvectors of the pinless
     material Hamiltonian at zero port. Inside a ring with a good neighbour cell, at
     finite K and nonzero port, no state of the invariant closure (illegal states
     included) changes a pinned data letter; only the pinned marker can change, through
     the intercell port with the neighbour's R0.
 (c) The quadratic-form bound (22) checked by minimising over (x, y)."""
import sys, os, json, itertools
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13

out = {}
# (a)
Bs = 2 ** 16
K = 2 ** 2560
eps = Fraction(Bs ** 2, K - 2 * Bs)
# rational lower bound on cos(pi/16) from the alternating Taylor series, with pi < 22/7 upper bound
x = Fraction(22, 7 * 16)   # > pi/16, so cos(x) < cos(pi/16)
cos_lb = 1 - x ** 2 / 2 + x ** 4 / 24 - x ** 6 / 720
out['r_pin'] = dict(two_cos_lower_bound=float(2 * cos_lb), log2_eps=float(np.log2(float(Bs ** 2)) - 2560),
                    r_pin_gt_1=bool(2 * cos_lb - eps > 1), r_pin_lower=float(2 * cos_lb - eps))
# norm estimate of a length-15 cell plus its left-marker bond, penalty removed, schedule-like values
p = 14
pm = my13.params(p, D=1000., Delta=4 - 2 * np.cos(np.pi / 16), g=2. ** -20, td=1e-6, K=0.,
                 musel=2 - 0.11, field=0.)
b = my13.b11(pm, 0.)
nb = float(np.linalg.norm(b, 2))
out['B_cell_estimate'] = dict(bond_norm=nb, bonds_in_cell=p + 1, crude_bound=(p + 1) * nb, below_2_16=(p + 1) * nb < Bs)
print(out, flush=True)

# (b)
res = []
for p in range(2, 7):
    pm = my13.params(p, K=50., musel=.7, field=.03, g=.2, td=.1)
    for lam in (0., .2):
        h, ls = my13.bond(pm, lam, J=1., t=1., mu=0.)
        nev = 0; maxdim = 0; changed = 0
        for bits in itertools.product((0, 1), repeat=p - 1):
            for m in my13.MARKERS:
                pin = tuple(f'L{a}' for a in bits) + ('L1', m)
                # frozen neighbour: all-L terminal good cell ending L0 E (no R0 next to the marker)
                ring = pin + tuple('L0' for _ in range(p)) + ('E',)
                st, H, pen, mb = my13.ring_closure(ring, h, ls)
                maxdim = max(maxdim, len(st)); nev += 1
                # mobile neighbour: vacuum head at the left of the neighbour (R0 next to the pinned marker)
                ring2 = pin + ('R0',) * p + ('E',)
                st2, H2, _, _ = my13.ring_closure(ring2, h, ls, cap=200000)
                changed += sum(w[:p] != pin[:p] for w in st2)
        res.append(dict(p=p, lam=lam, pinned_words=nev, max_closure_dim_neighbour_terminal_L0E=maxdim,
                        states_where_pinned_data_letters_changed_mobile_neighbour=changed))
        print(res[-1], flush=True)
out['pinned_words_frozen'] = res

# (c) min over x^2 + y^2 = 1 of mu x^2 - 2B x y + (K - B) y^2 - mu  versus -B^2/(K - 2B)
chk = []
for mu, B, KK in [(1.889, 2. ** 16, 2. ** 40), (1.889, 10., 100.), (0.5, 3., 9.)]:
    M = np.array([[mu, -B], [-B, KK - B]])
    emin = np.linalg.eigvalsh(M)[0]
    chk.append(dict(mu=mu, B=B, K=KK, exact_min_minus_mu=float(emin - mu), bound=float(-B * B / (KK - 2 * B)),
                    holds=bool(emin - mu >= -B * B / (KK - 2 * B) - 1e-12)))
out['quadratic_form'] = chk
print(chk)
json.dump(out, open(os.path.join(HERE, 'check_frozen13.json'), 'w'), indent=1)
