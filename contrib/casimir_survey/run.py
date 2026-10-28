#!/usr/bin/env python3
"""One chain, one ring length: low spectrum by (charge sector, momentum), E_0, k_0, the
half-ring entropy of the ground state, and ||h||. Appends a JSON line to --out.

  python3 run.py --model SU3 --L 12 --out data/calib.jsonl
  python3 run.py --model blbq --par 0.9 --L 12 --out data/blbq.jsonl

Models with a U(1)^r charge list the sectors to search (targets); 'half' momenta means
k = 0..L/2 (reflection- or time-reversal-symmetric bonds), 'all' means k = 0..L-1.
"""
import argparse, json, time
import numpy as np
import core, chains as ch


def colors(d, n=None):
    """U(1)^(d-1) color charges: state a < d-1 carries unit charge a."""
    q = np.zeros((d, d - 1), dtype=int)
    for a in range(d - 1):
        q[a, a] = 1
    return q


def sz_charges(S):
    return np.round(2 * np.arange(S, -S - 1, -1)).astype(int)[:, None]


def equal_colors(d, L):
    return [[L // d] * (d - 1)] if L % d == 0 else []


def build(model, par):
    """Return dict(h, d, charges, targets(L), momenta)."""
    p = par
    if model in ('XX', 'Heisenberg'):
        h = ch.models()[model][0]
        return dict(h=h, d=2, charges=sz_charges(0.5), targets=lambda L: [[L % 2]], momenta='half')
    if model == 'Ising':
        return dict(h=ch.models()['Ising'][0], d=2, charges=None, targets=lambda L: [None], momenta='half')
    if model in ('SU3', 'SU4', 'SU5'):
        d = int(model[2])
        return dict(h=ch.swap(d), d=d, charges=colors(d), targets=lambda L: equal_colors(d, L), momenta='half')
    if model == 'TB':                     # Takhtajan-Babujian spin S = p[0]
        S = p[0]
        return dict(h=ch.takhtajan_babujian(S), d=int(round(2 * S + 1)), charges=sz_charges(S),
                    targets=lambda L: [[int(round(2 * S * L)) % 2]], momenta='half')
    if model == 'blbq':                   # spin-1 cos(t) S.S + sin(t) (S.S)^2
        return dict(h=ch.blbq(p[0]), d=3, charges=sz_charges(1.0), targets=lambda L: [[0]], momenta='half')
    if model == 'sp4':                    # Sp(4) on C^4: cos(t) P5 + sin(t) P10
        import models as cm
        m = cm.sp4()
        h = np.cos(p[0]) * m['hs'][0] + np.sin(p[0]) * m['hs'][1]
        return dict(h=h, d=4, charges=m['charges'], targets=lambda L: [[L % 2, 0]], momenta='half')
    if model == 'spin32':                 # spin 3/2: sum_J a_J P_J, a = p (4 numbers)
        Ps = ch.spin_projectors(1.5)
        h = sum(c * P for c, P in zip(p, Ps))
        return dict(h=h, d=4, charges=sz_charges(1.5), targets=lambda L: [[L % 2]], momenta='half')
    if model == 'orbital':                # spin-orbital: p = (c_AB, c_A, c_B)
        import models as cm
        m = cm.spin_orbital()
        h = sum(c * x for c, x in zip(p, m['hs']))
        return dict(h=h.astype(complex), d=4, charges=m['charges'], targets=lambda L: [[L % 2, L % 2]], momenta='half')
    if model == 'son':                    # SO(N) vector: p = (N, theta), cos P + sin K
        import models as cm
        N = int(p[0]); m = cm.son_vector(N); r = N // 2
        h = np.cos(p[1]) * m['hs'][0] + np.sin(p[1]) * m['hs'][1]
        tg = (lambda L: [[0] * r]) if N % 2 else (lambda L: [[L % 2] + [0] * (r - 1)])
        return dict(h=h, d=N, charges=m['charges'], targets=tg, momenta='half')
    if model == 'symrep':                 # integrable SU(N) chain in Sym^k, p = (N, k)
        N, k = int(p[0]), int(p[1])
        h, d, q = ch.sym_rep_chain(N, k)
        return dict(h=h, d=d, charges=q, targets=lambda L: [[k * L // N] * (N - 1)] if (k * L) % N == 0 else [],
                    momenta='half')
    if model == 'graded':                 # s * P_g + sum_a mu_a n_a, p = (s, parities..., mus...)
        d = (len(p) - 1) // 2
        par = [int(x) for x in p[1:1 + d]]; mu = p[1 + d:]
        h = p[0] * ch.graded_swap(par) + sum(mu[a] * ch.onsite(ch.number_op(d, a), d) for a in range(d))
        return dict(h=h, d=d, charges=colors(d), targets=None, momenta='half')
    raise ValueError(model)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--par', type=float, nargs='*', default=[])
    ap.add_argument('--L', type=int, nargs='+', required=True)
    ap.add_argument('--nlev', type=int, default=2)
    ap.add_argument('--ks', default=None, help="'near' = only k in {0,1,L/2-1,L/2,L/2+1,L-1}")
    ap.add_argument('--targets', default=None, help='JSON list of charge targets, overrides the default')
    ap.add_argument('--noent', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    m = build(a.model, a.par)
    hn = core.hnorm(m['h'], m['d'])
    for L in a.L:
        t0 = time.time()
        tg = json.loads(a.targets) if a.targets else m['targets'](L)
        if a.ks == 'near':
            ks = sorted({0, 1, L - 1, L // 2 - 1, L // 2, (L // 2 + 1) % L})
        elif m['momenta'] == 'half':
            ks = list(range(L // 2 + 1))
        else:
            ks = None
        r = core.ring(m['h'], m['d'], L, charges=m['charges'], targets=tg, ks=ks, nlev=a.nlev, entropy=not a.noent)
        rec = dict(model=a.model, par=a.par, L=L, hnorm=hn, E0=r['E0'], k0=r['k0'],
                   target=(None if r['target'] is None else list(np.atleast_1d(r['target']).tolist())),
                   gap1=core.gap_at(r, 1), S_half=(r['S'][0] if 'S' in r else None),
                   spec=[[list(t) if t is not None else None, k, e] for (t, k), e in r['spec'].items()],
                   time=time.time() - t0)
        with open(a.out, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        print(a.model, a.par, L, 'E0=%.10f k0=%d gap1*L/2pi=%.5f S=%s t=%.1fs' % (
            rec['E0'], rec['k0'], rec['gap1'] * L / (2 * np.pi), rec['S_half'], rec['time']), flush=True)


if __name__ == '__main__':
    main()
