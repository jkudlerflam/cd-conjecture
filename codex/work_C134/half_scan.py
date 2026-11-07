"""C134 half-filled checks. Large mode is for SLURM, never run on import.

J=1 and M=1/2. For odd L, the doubled ring has n=L, not 2*floor(L/2).
Floating-point diagnostics are not certified energy enclosures.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from bethe_scan import solve, functions
from check_c129 import energy as ed_energy

DELTAS = [-.999, -.99, -.9, 0., .5, .9, .99, 1.]


def constants(delta):
    if not -1 < delta <= 0:
        return None
    g = np.arccos(delta)
    s = np.sin(g)
    q = 2*g/np.pi-1
    B = 1/(1-q)
    V1 = 2*abs(1/np.tan(g))/np.pi
    V2 = q*(8/s**4+4/s**2)
    D = 1+B*V1
    D1 = B*V2
    b, a = np.pi/g, np.pi*s
    Klarge = a*(1.5*b+.5*D1/D+32*D**2/b+8*D*B)
    K = max((4*g*D+1)**2, Klarge)
    return dict(gamma=g, q=q, B=B, D=D, D1=D1, K=K, C=2.5*K,
                Q_bound=5*K, cutoff_length=4*g*D)


def compact_constant(delta_gap):
    alpha = np.arccos(1-delta_gap)
    s = np.sqrt(delta_gap*(2-delta_gap))
    B = np.pi/(2*alpha)
    D = 1+1/(alpha*s)
    D1 = B*(8/s**4+4/s**2)
    K = (4*np.pi*D+1)**2+np.pi*(3+.5*D1+32*D**2+8*D*B)
    return dict(delta_gap=delta_gap, K=K, C=2.5*K, Q_bound=5*K)


def sech(x):
    u = np.exp(-np.abs(x))
    return 2*u/(1+u*u)


def bulk(delta):
    if delta == 1:
        return .25-np.log(2)
    g = np.arccos(delta)
    def integrand(x):
        t = np.tanh(x)
        bare = -np.sin(g)**2*(1-t*t)/(2*(t*t+np.sin(g/2)**2*(1-t*t)))
        return bare*sech(np.pi*x/g)/(2*g)
    return delta/4+2*quad(integrand, 0, 40, epsabs=3e-14, epsrel=3e-14)[0]


def root_checks(L, delta, roots):
    c = constants(delta)
    if c is None:
        return {}
    g, n = c['gamma'], len(roots)
    p, theta, dp, dt, bare, _ = functions(delta)
    xs = np.unique(np.r_[np.linspace(-max(8., abs(roots[0])+2),
                                    max(8., abs(roots[-1])+2), 1001), roots])
    cum_rho = .25+np.arctan(np.sinh(np.pi*xs/g))/(2*np.pi)
    before = np.searchsorted(roots, xs, side='left')/L
    after = np.searchsorted(roots, xs, side='right')/L
    max_G = max(np.max(np.abs(before-cum_rho)), np.max(np.abs(after-cum_rho)))
    sigma = (dp(xs)-dt(xs[:, None]-roots[None, :]).sum(axis=1)/L)/(2*np.pi)
    rho = sech(np.pi*xs/g)/(2*g)
    result = dict(cdf_scaled=float(L*max_G), cdf_bound=c['B'],
                  density_scaled=float(L*np.max(abs(sigma-rho))), density_bound=c['D'],
                  minimum_sigma=float(np.min(sigma)))
    assert L*max_G <= c['B']+1e-9
    assert L*np.max(abs(sigma-rho)) <= c['D']+1e-9
    if L > c['cutoff_length']:
        A = g/np.pi*np.arccosh(L/(4*g*c['D']))
        grid = np.linspace(-A, A, 4001)
        sig = (dp(grid)-dt(grid[:, None]-roots[None, :]).sum(axis=1)/L)/(2*np.pi)
        rh = sech(np.pi*grid/g)/(2*g)
        b, a = np.pi/g, np.pi*np.sin(g)
        fprime = a*b*np.tanh(b*grid)*rh/sig
        variation = abs(fprime[0])+np.sum(abs(np.diff(fprime)))+abs(fprime[-1])
        varbound = 12*a*b+4*a*c['D1']/c['D']
        result.update(cutoff=A, sampled_variation=float(variation), variation_bound=varbound)
        assert variation <= varbound+1e-8
    return result


def run(output, lengths, deltas, check_ed=False, check_roots=False):
    cache = {}
    failures = 0
    rows = []
    def get(L, delta):
        key = (L, delta)
        if key not in cache:
            cache[key] = solve(L, L//2, delta, return_roots=True)
        return cache[key]
    with Path(output).open('w') as handle:
        for delta in deltas:
            e0 = bulk(delta)
            for L in lengths:
                try:
                    E, diag, roots = get(L, delta)
                    E2, diag2, _ = get(2*L, delta)
                    Q = float(2*L*(E2-2*E))
                    row = dict(status='ok', L=L, n=L//2, n2=L, delta=delta,
                               E=float(E), E2=float(E2), Q=Q, bulk=e0,
                               lower_casimir=float(L*(L*e0-E)), diag=diag, diag2=diag2,
                               correction_indicator=2*L*(diag2['newton_energy_correction']+
                                                         2*diag['newton_energy_correction']))
                    c = constants(delta)
                    if c is not None:
                        row['bound'] = c
                        assert Q <= c['Q_bound']+1e-8
                        assert row['lower_casimir'] <= c['K']+1e-8
                        if L % 2 == 0:
                            assert abs(row['lower_casimir']) <= c['K']+1e-8
                    if check_ed and L <= 12:
                        row['ed_error'] = abs(float(E)-ed_energy(L, L//2, delta))
                        row['spinflip_error'] = abs(ed_energy(L, L//2, delta)-
                                                   ed_energy(L, L-L//2, delta))
                        assert row['ed_error'] < 2e-9
                        assert row['spinflip_error'] < 2e-9
                    if check_roots:
                        row['root_checks'] = root_checks(L, delta, roots)
                    if delta == 0:
                        # Analytic cancellation-free reference for both parities.
                        expected = (2*L*np.tan(np.pi/(4*L))/np.cos(np.pi/(2*L))
                                    if L % 2 == 0 else 0.)
                        row['free_Q_error'] = abs(Q-expected)
                        assert abs(Q-expected) < 2e-7
                except Exception as exc:
                    failures += 1
                    row = dict(status='failed', L=L, delta=delta, error=str(exc))
                rows.append(row)
                handle.write(json.dumps(row)+'\n')
                handle.flush()
            print(f'Delta={delta} complete', flush=True)
    ok = [r for r in rows if r['status'] == 'ok']
    summary = dict(cases=len(rows), failures=failures, max_L=max(lengths),
                   max_Bethe_length=2*max(lengths),
                   max_ed_error=max([r.get('ed_error', 0) for r in ok], default=0),
                   max_spinflip_error=max([r.get('spinflip_error', 0) for r in ok], default=0),
                   max_correction_indicator=max([r['correction_indicator'] for r in ok], default=0),
                   max_Q=max([r['Q'] for r in ok], default=None),
                   per_delta=[dict(delta=d, max_Q=max([r['Q'] for r in ok if r['delta']==d], default=None),
                                   last=[dict(L=r['L'], Q=r['Q']) for r in ok if r['delta']==d][-2:])
                              for d in deltas])
    Path(output+'.summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['small', 'moderate', 'large'], required=True)
    parser.add_argument('--index', type=int, choices=range(len(DELTAS)))
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.mode == 'small':
        run(args.output, list(range(2, 13)), DELTAS, True, True)
    elif args.mode == 'moderate':
        run(args.output, [16, 17, 32, 33, 64, 65, 128, 129], DELTAS, False, True)
    else:
        if args.index is None:
            parser.error('Large mode needs --index for a SLURM array job.')
        lengths = sorted([v for L in [256, 512, 1024, 2048, 4096] for v in [L, L+1]])
        run(args.output, lengths, [DELTAS[args.index]])
