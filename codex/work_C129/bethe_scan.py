"""Numerical logarithmic Bethe equations at the analytically minimizing flux.

J=1, S=sigma/2. The physical minimizing flux is L*pi modulo 2*pi.
The centered quantum numbers encode this flux for either parity of L.
No large scan runs on import. Use --small locally, the SLURM file for large work.
Failures are output explicitly and are never replaced by interpolated energies.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from check_c129 import energy as ed_energy, bond_M

DELTAS = [-.999, -.99, -.9, -.5, 0., .5, .9, 1., 1.5, 3.]


def functions(delta):
    if -1 < delta < 1:
        z = np.arccos(delta)
        def phase(x, angle):
            return 2*np.arctan(np.tanh(x)/np.tan(angle/2))
        def deriv(x, angle):
            t, u = np.tan(angle/2), np.tanh(x)
            return 2*t*(1-u*u)/(t*t+u*u)
        p = lambda x: phase(x, z)
        t = lambda x: phase(x, 2*z)
        dp = lambda x: deriv(x, z)
        dt = lambda x: deriv(x, 2*z)
        bare = lambda x: -np.sin(z)**2/(2*(np.sinh(x)**2+np.sin(z/2)**2))
        seed = lambda I, L, n: z/np.pi*np.arcsinh(np.tan(2*np.pi*I/L))
    elif delta == 1:
        p = lambda x: 2*np.arctan(2*x)
        t = lambda x: 2*np.arctan(x)
        dp = lambda x: 4/(1+4*x*x)
        dt = lambda x: 2/(1+x*x)
        bare = lambda x: -.5/(x*x+.25)
        seed = lambda I, L, n: np.arcsinh(np.tan(2*np.pi*I/L))/np.pi
    elif delta > 1:
        z = np.arccosh(delta)
        def phase(x, scale):
            periods = np.floor(x/np.pi+.5)
            y = x-periods*np.pi
            return 2*np.arctan(np.tan(y)/np.tanh(scale/2))+2*np.pi*periods
        p = lambda x: phase(x, z)
        t = lambda x: phase(x, 2*z)
        dp = lambda x: np.sinh(z)/(np.sinh(z/2)**2+np.sin(x)**2)
        dt = lambda x: np.sinh(2*z)/(np.sinh(z)**2+np.sin(x)**2)
        bare = lambda x: -np.sinh(z)**2/(2*(np.sinh(z/2)**2+np.sin(x)**2))
        seed = lambda I, L, n: np.pi*I/(L-n)
    else:
        raise ValueError('Real-root solver is restricted to Delta > -1')
    return p, t, dp, dt, bare, seed


def solve(L, n, delta, tol=4e-14, iterative_threshold=256):
    if n == 0:
        return L*delta/4, dict(residual=0., steps=0, newton_energy_correction=0.)
    I = np.arange(1, n+1, dtype=float)-(n+1)/2
    p, theta, dp, dt, bare, seed = functions(delta)
    x = seed(I, L, n)
    def residual(y):
        return p(y)-theta(y[:, None]-y[None, :]).sum(axis=1)/L-2*np.pi*I/L
    def correction(y, rhs):
        W = dt(y[:, None]-y[None, :])/L
        diagonal = dp(y)-W.sum(axis=1)
        if n <= iterative_threshold:
            return np.linalg.solve(np.diag(diagonal)+W, rhs)
        op = LinearOperator((n, n), matvec=lambda v: diagonal*v+W@v)
        jd = diagonal+np.diag(W)
        pre = LinearOperator((n, n), matvec=lambda v: v/jd)
        result, info = gmres(op, rhs, M=pre, rtol=1e-11, atol=1e-16,
                            restart=80, maxiter=100)
        if info:
            raise RuntimeError(f'GMRES did not converge: {info}')
        return result
    for step in range(100):
        r = residual(x)
        norm = float(np.max(np.abs(r)))
        if norm < tol:
            break
        dx = correction(x, -r)
        for halves in range(35):
            trial = x+dx*2.**(-halves)
            trial = (trial-trial[::-1])/2
            if np.max(np.abs(residual(trial))) < norm:
                x = trial
                break
        else:
            raise RuntimeError(f'Newton line search failed, residual={norm}')
    else:
        raise RuntimeError(f'Newton iteration limit, residual={norm}')
    if n > 1 and np.min(np.diff(x)) <= 0:
        raise RuntimeError('Roots are not strictly ordered')
    # Sum in long double to reduce subtraction error in the doubling defect.
    E = np.longdouble(L)*np.longdouble(delta)/4 + np.sum(bare(x).astype(np.longdouble))
    dx = correction(x, -residual(x))
    Ecorr = float(abs(np.sum((bare(x+dx)-bare(x)).astype(np.longdouble))))
    return E, dict(residual=norm, steps=step, newton_energy_correction=Ecorr,
                   largest_root=float(np.max(np.abs(x))))


def write_row(handle, row):
    handle.write(json.dumps(row)+'\n')
    handle.flush()


def small(output):
    worst, error, cases, failures = (-np.inf, None), 0., 0, 0
    with Path(output).open('w') as f:
        for L in range(2, 11):
            for n in range(1, L//2+1):
                for delta in DELTAS:
                    try:
                        E, diag = solve(L, n, delta)
                        reference = ed_energy(L, n, delta)
                        diff = abs(float(E)-reference)
                        error = max(error, diff)
                        if diff > 2e-9:
                            raise RuntimeError(f'ED mismatch {diff}')
                        E2, diag2 = solve(2*L, 2*n, delta)
                        Q = float(L*(E2-2*E)/bond_M(1, delta))
                        cases += 1
                        if Q > worst[0]:
                            worst = (Q, [L, n, delta])
                        write_row(f, dict(L=L, n=n, delta=delta, E=float(E), E2=float(E2),
                                         Q=Q, ed_error=diff, diag=diag, diag2=diag2, status='ok'))
                    except Exception as exc:
                        failures += 1
                        write_row(f, dict(L=L, n=n, delta=delta, status='failed', error=str(exc)))
            print(f'Bethe vs ED L={L} complete', flush=True)
    summary = dict(successful_cases=cases, failures=failures, maximum_ed_error=error,
                   largest_Q=worst[0], at=worst[1], ED_maxL=10, Bethe_maxL=20)
    print(json.dumps(summary, indent=2))
    Path(output+'.summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if failures:
        raise SystemExit(1)


def large(output, index):
    delta0 = DELTAS[index] if index < len(DELTAS) else None
    with Path(output).open('w') as f:
        for L in [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
            ns = sorted(set([1, 2, 3, 4, L//2, L//2-1, int(np.sqrt(L))]
                            + [max(1, int(round(L*r))) for r in [.01, .02, .05, .1, .2, .3, .4, .49]]))
            ds = [delta0] if delta0 is not None else [-1+c/L**2 for c in [0.1, 1, 10, 100]]
            for delta in ds:
                for n in ns:
                    if not 1 <= n <= L//2:
                        continue
                    try:
                        E, diag = solve(L, n, delta)
                        E2, diag2 = solve(2*L, 2*n, delta)
                        M = bond_M(1, delta)
                        Q = float(L*(E2-2*E)/M)
                        write_row(f, dict(status='ok', L=L, n=n, delta=delta, Q=Q,
                                         E=float(E), E2=float(E2), diag=diag, diag2=diag2,
                                         correction_indicator=L*(diag2['newton_energy_correction']+
                                                                 2*diag['newton_energy_correction'])/M))
                    except Exception as exc:
                        write_row(f, dict(status='failed', L=L, n=n, delta=delta, error=str(exc)))
            print(f'L={L} complete', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--small', action='store_true')
    parser.add_argument('--index', type=int, choices=range(11))
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.small:
        small(args.output)
    elif args.index is not None:
        large(args.output, args.index)
    else:
        parser.error('Choose --small or --index. Large scans belong on the cluster.')
