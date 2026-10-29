"""Payload dictionary in the slow-vacancy regime (numerical, small sizes).

Frozen: three C84 cells (p=2, one payload bit each) on a 9-site ring, K=infinity.
Mobile: the same cells plus one vacancy on a 10-site ring, d=12 bond, K=infinity.
With payload << t_eff the vacancy sits in its q=0 band state and each of the
three inter-cell bonds is open with probability 2/3.  Two predictions:

 (a) port-induced shift of the zero-fermion ground state (the payload-blind
     lambda^2 scalar part, a pure contact term):
       [E_m(N_f=0,lam)-E_m(N_f=0,0)] / [E_f(N_f=0,lam)-E_f(N_f=0,0)] -> 2/3
 (b) the three one-fermion levels measured from the zero-fermion ground state
     (the payload-specific part; the scalar part cancels):
       [E_m(N_f=1)_i-E_m(N_f=0)] / [E_f(N_f=1)_i-E_f(N_f=0)] -> 2/3
 both up to O(t^2) clock distortion and O(payload/t_eff).
Diagnostic couplings (g=0.5, lambda=0.03, t_d=0.3): lambda is small enough that the
port scalar energy (about 3e-5) lies below t_eff, and g, t_d are large so that the
payload band (about 5e-12) is resolvable in double precision.
"""
import os, json
import numpy as np
from scipy.sparse.linalg import eigsh
import bond12 as b
from closure import closure

HERE = os.path.dirname(os.path.abspath(__file__))
L0, L1, R0, R1, E, V = ('L', 0), ('L', 1), ('R', 0), ('R', 1), ('e', 0), b.V


def e0(seed, p, lam, t, k=3):
    """Lowest k levels.  Two-stage shift-invert (coarse ground energy, then
    sigma 1e-8 below it) so that a cluster split by ~1e-12 is resolved."""
    H, _, _, st = closure(seed, [b.bond(p, lam=lam, t=t)] * len(seed), legal_only=True)
    if H.shape[0] <= 3000:
        return np.linalg.eigvalsh(H.toarray())[:k], len(st)
    H = H.tocsc()
    sig = -float(np.max(np.asarray(abs(H).sum(axis=1)))) - 0.5
    g0 = eigsh(H, k=1, sigma=sig, which='LM', tol=1e-14)[0][0]
    ev = eigsh(H, k=k, sigma=g0 - 1e-8, which='LM', tol=1e-15, ncv=40)[0]
    return np.sort(ev), len(st)


g, lam, td = 0.5, 0.03, 0.3
cell0 = [R0, R0, E]; cell1 = [R0, R1, E]
seeds = {('f', 0): cell0 * 3, ('f', 1): cell1 + cell0 * 2,
         ('m', 0): cell0 + [V] + cell0 * 2, ('m', 1): cell1 + [V] + cell0 * 2}
out = []
for t, w in [(0.1, 0.), (0.07, 0.), (0.05, 0.), (0.07, -0.02)]:
    p = b.params(2, K=0., D=4., g=g, td=td, w=w)
    E_ = {}; dims = {}
    for key, s in seeds.items():
        tt = 0. if key[0] == 'f' else t
        E_[key + ('lam',)], dims[key] = e0(s, p, lam, tt, k=(3 if key[1] == 1 else 1))
        if key[1] == 0:
            E_[key + ('0',)], _ = e0(s, p, 0., tt, k=1)
    shift_f = E_[('f', 0, 'lam')][0] - E_[('f', 0, '0')][0]
    shift_m = E_[('m', 0, 'lam')][0] - E_[('m', 0, '0')][0]
    band_f = E_[('f', 1, 'lam')] - E_[('f', 0, 'lam')][0]   # the three one-fermion levels
    band_m = E_[('m', 1, 'lam')] - E_[('m', 0, 'lam')][0]
    rec = dict(t=t, w=w, g=g, lam=lam, td=td, t_eff_estimate=0.5 * t ** 3,
               dims={'%s%d' % k: v for k, v in dims.items()},
               scalar_shift_frozen=shift_f, scalar_shift_mobile=shift_m, ratio_scalar=shift_m / shift_f,
               one_fermion_levels_frozen=[float(x) for x in band_f],
               one_fermion_levels_mobile=[float(x) for x in band_m],
               ratios_one_fermion=[float(x) for x in band_m / band_f], predicted_ratio=2 / 3)
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open(os.path.join(HERE, 'dictionary_checks.json'), 'w'), indent=1)
