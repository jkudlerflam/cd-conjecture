"""(8) Floating-point recomputation of the midpoint selectors (C84 eq. (18), b = 13;
C85 eq. (19), b = 15): mu_sel, the density margin min_(l != b)(e_l - e_b) over
l = 1..10^5, e_b, the vacuum gap, and the F3 condition gamma > 5 sqrt(2) lambda."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
a = lambda l: 2 * (1 - np.cos(np.pi / (l + 1)))
dl = lambda l: (l + 1) * a(l) - l * a(l + 1)
out = {}
for tag, b, lam_pow, K_pow in [('C84', 13, -480, 2048), ('C85', 15, -608, 2560)]:
    mu = 2 - (dl(b) + dl(b - 1)) / 2
    ls = np.arange(1, 100001)
    e = (mu - 2 + a(ls)) / ls
    eb = e[b - 1]
    marg = np.min(np.delete(e, b - 1) - eb)
    p = b - 1
    gap = 2 * (np.cos(np.pi / (p + 2)) - np.cos(2 * np.pi / (p + 2)))
    out[tag] = dict(b=b, mu_sel=float(mu), density_margin=float(marg), e_b=float(eb), vacuum_gap=float(gap),
                    argmin_competitor=int(ls[np.argsort(e)[1]]),
                    log2_gamma_minus_log2_5sqrt2_lambda=float(np.log2(marg) - (np.log2(5 * np.sqrt(2)) + lam_pow)))
    print(tag, out[tag])
json.dump(out, open(os.path.join(HERE, 'check_selector_numbers.json'), 'w'), indent=1)
