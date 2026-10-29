#!/usr/bin/env python3
"""Numerical checks of the structural lemmas in README.md (stacking, twisting, gauge
triviality for reflection-symmetric bonds, blocking).  Writes data/lemmas_check.json."""
import json
import numpy as np
from core import hnorm, spread, Ring, fit_cas, LOG2
from families import heisenberg, swap

rng = np.random.default_rng(7)
out = {}

sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.diag([.5, -.5])
SS = sum(np.kron(s, s) for s in (sx, sy, sz)).real


def haar_su2():
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    q = q * (np.diag(r) / abs(np.diag(r)))
    return q / np.sqrt(np.linalg.det(q))


# Lemma 1: stacking.  M_stack >= max(M_i); here Heisenberg (x) Heisenberg and Heisenberg (x) random.
def stack(h1, d1, h2, d2):
    h = np.kron(h1, np.eye(d2 * d2)).reshape(d1, d1, d2, d2, d1, d1, d2, d2)
    h = h.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape((d1 * d2) ** 2, (d1 * d2) ** 2)
    g = np.kron(np.eye(d1 * d1), h2).reshape(d1, d1, d2, d2, d1, d1, d2, d2)
    g = g.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape((d1 * d2) ** 2, (d1 * d2) ** 2)
    return h + g

hh = stack(SS, 2, SS, 2)
M_hh = hnorm(hh, 4)
out['stack_heis_heis'] = {'M_stack': M_hh, 'M_sum': 1.0, 'A': np.pi, 'ratio': np.pi / (np.pi * M_hh * 2)}
stk = []
for _ in range(5):
    r = rng.normal(size=(4, 4)); r = r + r.T
    h2 = r / spread(r)
    M1, M2 = hnorm(SS, 2), hnorm(h2, 2)
    Ms = hnorm(stack(SS, 2, h2, 2), 4)
    stk.append({'M1': M1, 'M2': M2, 'M_stack': Ms, 'ge_max': bool(Ms >= max(M1, M2) - 1e-7), 'le_sum': bool(Ms <= M1 + M2 + 1e-7)})
out['stack_random'] = stk

# Lemma 3: twisted Heisenberg h_w = (1 (x) w) h (1 (x) w^dag), w in SU(2): M_w = 1/2 for all w.
tw = []
for _ in range(8):
    w = haar_su2()
    U = np.kron(np.eye(2), w)
    hw = U @ SS @ U.conj().T
    tw.append({'M_w': hnorm(hw, 2), 'spread_half': spread(hw) / 2})
# also with w of finite order so that rings of commensurate L carry no boundary twist: w = exp(i pi sz) (order 2 up to sign)
w = np.diag([1j, -1j])
U = np.kron(np.eye(2), w); hw = (U @ SS @ U.conj().T)
E = [Ring(hw, 2, L).lowest(k=1)[0] for L in (8, 10, 12, 14)]
E0 = [Ring(SS, 2, L).lowest(k=1)[0] for L in (8, 10, 12, 14)]
out['twist'] = {'random_w': tw, 'w_order2_M': hnorm(hw, 2), 'E_twisted_minus_E_heis': [float(a - b) for a, b in zip(E, E0)],
                'note': 'w^L = +-1 for even L, so the ring energies coincide with the Heisenberg ring (up to the sign twist for L = 2 mod 4)'}

# Lemma 2: reflection-symmetric bonds with nondegenerate extremes have M = half spread.
refl = []
P2 = swap(2).real
for _ in range(6):
    r = rng.normal(size=(4, 4)); r = r + r.T; r = r + P2 @ r @ P2
    e = np.linalg.eigvalsh(r)
    refl.append({'M': hnorm(r, 2), 'half_spread': spread(r) / 2, 'nondegenerate': bool(e[1] - e[0] > 1e-6 and e[-1] - e[-2] > 1e-6)})
# and a non-reflection-symmetric control, where the gauge does help
ctrl = []
for _ in range(6):
    r = rng.normal(size=(4, 4)); r = r + r.T
    ctrl.append({'M': hnorm(r, 2), 'half_spread': spread(r) / 2})
out['reflection_symmetric'] = refl
out['no_reflection_control'] = ctrl

# Lemma 4: blocking b = 2 of the Heisenberg chain (d = 4): A' = 2 A, M' from the blocked bond.
hb = stack(SS, 2, np.zeros((4, 4)), 2)  # placeholder shape; build properly below
# blocked bond on sites (1,2 | 3,4): h_{23} + (h_{12} + h_{34})/2
def two_site_on(h, pos, n=4, d=2):
    full = np.eye(d ** n).reshape([d] * (2 * n))
    h4 = h.reshape(d, d, d, d)
    # apply h on axes pos within the n-site space
    op = np.zeros((d ** n, d ** n))
    for a in range(d ** n):
        v = np.zeros(d ** n); v[a] = 1
        psi = v.reshape([d] * n)
        t = np.tensordot(h4, psi, axes=([2, 3], list(pos)))
        t = np.moveaxis(t, [0, 1], list(pos))
        op[:, a] = t.reshape(-1)
    return op
hb = two_site_on(SS, (1, 2)) + 0.5 * (two_site_on(SS, (0, 1)) + two_site_on(SS, (2, 3)))
Mb = hnorm(hb, 4)
Eb = [Ring(hb, 4, L).lowest(k=1)[0] for L in (4, 5, 6, 7)]
Ab = fit_cas([4, 5, 6, 7], Eb)[0]
out['blocking_b2'] = {'M_blocked': Mb, 'A_blocked_fit': Ab, 'A_expected': np.pi, 'ratio': Ab / (np.pi * Mb * 2), 'ratio_expected_le': 0.5 / Mb}

json.dump(out, open('data/lemmas_check.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
