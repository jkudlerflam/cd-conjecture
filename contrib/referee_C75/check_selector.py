"""(3) Zero-port selector at finite K: full single-cell Hilbert spaces
(every word, illegal ones and any number of excited tokens), and marker-free
rings. Referee rebuild only."""
import itertools, json
import numpy as np
from scipy.sparse.linalg import eigsh
import my15
from my15 import LETTERS

DATA = [l for l in LETTERS if l[0] != 'M']
out = {'cells': [], 'rings': [], 'selector_float': {}}

def lowest(H):
    N = H.shape[0]
    if N <= 3000:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    return float(eigsh(H, k=1, which='SA', tol=1e-12, ncv=40, return_eigenvectors=False)[0])

for Delta, tag in [(4 - 2 * np.cos(np.pi / 17), 'C75_Delta'), (1.5, 'Delta_below_2_hypothesis_violated')]:
    for k in range(0, 5):
        if tag != 'C75_Delta' and k < 2: continue
        words = [w + (m,) for w in itertools.product(DATA, repeat=k) for m in my15.MARKERS]
        for K in [6., 10., 100.]:
            p = my15.default_params(max(k, 1), D=10., Delta=Delta, g=0.01, td=0.3, K=K,
                                    phiA=0.37, psiA=1.1, psiB=0.4, musel=0., field=0.)
            bonds = [(-1, 0)] + [(i, i + 1) for i in range(k)]
            st, H, _ = my15.closure(words, p, bonds, project=False, left='M0')
            assert len(st) == len(words)
            E = lowest(H)
            Ev = -2 * np.cos(np.pi / (k + 2))
            r = dict(case=tag, k=k, K=K, dim=len(st), E=E, upper=Ev, lower=Ev - max(k - 1, 0) / (K - 4),
                     within_bounds=bool(Ev - max(k - 1, 0) / (K - 4) - 1e-10 <= E <= Ev + 1e-10))
            out['cells'].append(r); print(r, flush=True)

for N in [2, 3, 4]:
    words = list(itertools.product(DATA, repeat=N))
    for K in [6., 10., 100.]:
        p = my15.default_params(2, D=10., Delta=4 - 2 * np.cos(np.pi / 17), g=0.01, td=0.3, K=K)
        bonds = [(i, (i + 1) % N) for i in range(N)]
        st, H, _ = my15.closure(words, p, bonds, project=False)
        E = lowest(H)
        r = dict(ring=N, K=K, dim=len(st), E=E, claimed_lower=-N / (K - 4), ok=bool(E >= -N / (K - 4) - 1e-10))
        out['rings'].append(r); print(r, flush=True)

# inherited midpoint selector, float cross-check of the rational inequalities
a = lambda b: 2 * (1 - np.cos(np.pi / (b + 1)))
dl = lambda b: (b + 1) * a(b) - b * a(b + 1)
mu = 2 - (dl(16) + dl(15)) / 2
e = {b: (mu - 2 * np.cos(np.pi / (b + 1))) / b for b in range(1, 2001)}
bmin = min(e, key=e.get)
gaps = sorted(e[b] - e[16] for b in e if b != 16)
out['selector_float'] = dict(mu_sel=mu, argmin_b=bmin, e16=e[16], min_density_gap=gaps[0],
                             gap_over_2m16=gaps[0] * 2 ** 16, inv_K_minus_4_below_2m17=bool(2 ** 1152 - 4 > 2 ** 17))
print(out['selector_float'])
json.dump(out, open('check_selector.json', 'w'), indent=1)
