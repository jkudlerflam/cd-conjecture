"""C112 checks. Large runs require --cluster and a SLURM allocation.

Only NumPy and SciPy are required. No symmetry restriction is made for
generic bonds. XXZ uses its balanced sector and labels that restriction.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, eigvalsh
from scipy.optimize import minimize
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import logsumexp


def gradient(w):
    return float(np.sum((np.roll(w, -1)-w)**2/w))


def profile(length, shape, eps, mode=2):
    x = np.arange(length, dtype=float)
    if shape == 'cosine':
        f = np.cos(2*np.pi*mode*(x+.5)/length)
    else:
        sigma = length/10
        f = sum(np.exp(-.5*((x+.5-length/2+j*length)/sigma)**2)
                for j in [-2, -1, 0, 1, 2])
    return 1+eps*f


def cft_response(w, A):
    """Spectral continuum integral for a sampled smooth periodic profile."""
    length = len(w)
    freq = 2*np.pi*np.fft.fftfreq(length)
    wp = np.fft.ifft(1j*freq*np.fft.fft(w)).real
    G = float(np.sum(wp**2/w))
    ell = float(np.sum(1/w))
    return A*G/(24*np.pi)+np.pi*A/6*(1/ell-w.sum()/length**2)


def spin_bond(delta):
    x = np.array([[0, 1], [1, 0]])/2
    y = np.array([[0, -1j], [1j, 0]])/2
    z = np.diag([1., -1.])/2
    return (np.kron(x, x)+np.kron(y, y)+delta*np.kron(z, z)).real


def basis(d, length, balanced=False):
    powers = d**np.arange(length, dtype=np.int64)
    if balanced:
        assert length % d == 0
        def assignments(remaining, color, code):
            if color == d-1:
                yield code+sum(color*powers[j] for j in remaining)
                return
            for sites in itertools.combinations(remaining, length//d):
                yield from assignments(tuple(j for j in remaining if j not in sites),
                                       color+1, code+sum(color*powers[j] for j in sites))
        codes = np.array(sorted(assignments(tuple(range(length)), 0, 0)))
    else:
        codes = np.arange(d**length, dtype=np.int64)
    digits = (codes[:, None]//powers) % d
    return codes, digits, powers


def bond_operators(h, length, balanced=False):
    d = int(round(np.sqrt(h.shape[0])))
    codes, digits, powers = basis(d, length, balanced)
    n = len(codes)
    ops = []
    for j in range(length):
        k = (j+1) % length
        inputs = d*digits[:, j]+digits[:, k]
        rr, cc, vv = [], [], []
        for output in range(d*d):
            values = h[output, inputs]
            ix = np.flatnonzero(np.abs(values) > 1e-14)
            targets = codes[ix]+(output//d-digits[ix, j])*powers[j]+(output % d-digits[ix, k])*powers[k]
            rows = np.searchsorted(codes, targets)
            if len(rows) and (rows.max() >= n or np.any(codes[rows] != targets)):
                raise ValueError('Bond does not preserve the selected sector.')
            rr.extend(rows)
            cc.extend(ix)
            vv.extend(values[ix])
        ops.append(coo_matrix((vv, (rr, cc)), shape=(n, n)).tocsr())
    return ops


def ground(ops, w):
    H = sum(a*b for a, b in zip(w, ops))
    n = H.shape[0]
    if n <= 64:
        ev, u = eigh(H.toarray(), subset_by_index=[0, 1])
    else:
        v0 = np.cos(np.arange(n)*np.sqrt(2))+np.sin(np.arange(n)*np.sqrt(3))
        ev, u = eigsh(H, k=2, which='SA', tol=2e-12, v0=v0, maxiter=20000)
        ix = np.argsort(ev)
        ev, u = ev[ix], u[:, ix]
    residual = np.linalg.norm(H@u[:, 0]-ev[0]*u[:, 0])
    return float(ev[0]), float(ev[1]-ev[0]), float(residual)


def hermitian_basis(d):
    out = []
    for j in range(d):
        for k in range(j+1, d):
            a = np.zeros((d, d), dtype=complex)
            a[j, k] = a[k, j] = 1/np.sqrt(2)
            out.append(a)
            a = np.zeros((d, d), dtype=complex)
            a[j, k], a[k, j] = 1j/np.sqrt(2), -1j/np.sqrt(2)
            out.append(a)
    for j in range(1, d):
        a = np.diag([1.]*j+[-float(j)]+[0.]*(d-j-1))/np.sqrt(j*(j+1))
        out.append(a)
    return np.array(out)


def gauge_optimize(h):
    """Return a feasible primal representative and a dual lower bound for M.

    The dual is explicitly projected to equal partial traces and normalized
    in trace norm. The interval is subject to floating-point roundoff.
    """
    d = int(round(np.sqrt(len(h))))
    eye = np.eye(d)
    gens = np.array([np.kron(o, eye)-np.kron(eye, o) for o in hermitian_basis(d)])
    p = np.zeros(len(gens))
    for beta in [8., 32., 128., 512., 2048.]:
        def obj(x):
            mat = h+np.einsum('a,aij->ij', x, gens)
            e, u = eigh(mat)
            pos = np.exp(beta*e-logsumexp(beta*e))
            neg = np.exp(-beta*e-logsumexp(-beta*e))
            dual = (u*((pos-neg)/2))@u.conj().T
            val = (logsumexp(beta*e)+logsumexp(-beta*e))/(2*beta)
            jac = np.einsum('ij,aji->a', dual, gens).real
            return val, jac
        fit = minimize(obj, p, jac=True, method='BFGS', options={'gtol': 2e-9, 'maxiter': 250})
        p = fit.x
    opt = h+np.einsum('a,aij->ij', p, gens)
    e, u = eigh(opt)
    pos = np.exp(beta*e-logsumexp(beta*e))
    neg = np.exp(-beta*e-logsumexp(-beta*e))
    dual = (u*((pos-neg)/2))@u.conj().T
    tensor = dual.reshape(d, d, d, d)
    r = np.trace(tensor, axis1=1, axis2=3)-np.trace(tensor, axis1=0, axis2=2)
    dual -= (np.kron(r, eye)-np.kron(eye, r))/(2*d)
    dual /= np.sum(np.abs(eigvalsh(dual)))
    lower = float(np.trace(h@dual).real)
    upper = float((e[-1]-e[0])/2)
    assert lower <= upper+1e-10
    return opt, lower, upper


def random_checks(cluster=False):
    rng = np.random.default_rng(1122026)
    rows = []
    for d in [2, 3, 4]:
        lengths = {2: [12, 16], 3: [7, 9], 4: [6, 7]}[d] if cluster else {2: [8], 3: [5], 4: [4]}[d]
        for draw in range(24 if cluster else 3):
            z = rng.normal(size=(d*d, d*d))+1j*rng.normal(size=(d*d, d*d))
            h = (z+z.conj().T)/2
            h -= np.trace(h)*np.eye(d*d)/(d*d)
            h /= (eigvalsh(h)[-1]-eigvalsh(h)[0])/2
            opt, ml, mu = gauge_optimize(h)
            for representative, bond in [('raw', h), ('optimized', opt)]:
                for length in lengths:
                    ops = bond_operators(bond, length)
                    e0, gap, res0 = ground(ops, np.ones(length))
                    for shape in ['cosine', 'gaussian']:
                        w = profile(length, shape, .15, 1)
                        ew, gapw, resw = ground(ops, w)
                        deficit = e0*w.mean()-ew
                        G = gradient(w)
                        rows.append(dict(d=d, draw=draw, L=length, dimension=d**length,
                                         representative=representative, shape=shape, epsilon=.15,
                                         mode=1, M_lower=ml, M_upper=mu,
                                         uniform_energy=e0, weighted_energy=ew,
                                         deficit=deficit, G=G, ratio_lower=deficit/(mu*G),
                                         ratio_upper=deficit/(ml*G), gap=gap, weighted_gap=gapw,
                                         residual=max(res0, resw),
                                         bond_real=bond.real.tolist(), bond_imag=bond.imag.tolist()))
                        assert deficit >= -2e-9
    return rows


def xxz_checks(cluster=False):
    rows = []
    for delta in [-1.01, -1., -.9999, -.99, -.9, 0., .5, 1., 1.5]:
        A = np.pi/2 if delta == 1 else (np.pi*np.sqrt(1-delta**2)/(2*np.arccos(delta)) if -1 < delta < 1 else 0.)
        M = (eigvalsh(spin_bond(delta))[-1]-eigvalsh(spin_bond(delta))[0])/2
        for length in ([16, 20, 24] if cluster else [8, 10, 12]):
            ops = bond_operators(spin_bond(delta), length, balanced=True)
            e0, gap, r0 = ground(ops, np.ones(length))
            if delta <= -1:
                e0 = min(e0, delta*length/4)
            for shape, mode in [('cosine', 1), ('cosine', 2), ('gaussian', 2)]:
                for eps in [.02, .08]:
                    w = profile(length, shape, eps, mode)
                    ew, gapw, rw = ground(ops, w)
                    if delta <= -1:
                        ew = min(ew, delta*w.sum()/4)
                    deficit = e0*w.mean()-ew
                    prediction = cft_response(w, A)
                    rows.append(dict(delta=delta, L=length, dimension=ops[0].shape[0],
                                     sector='balanced, plus polarized energies for Delta <= -1',
                                     shape=shape, mode=mode, epsilon=eps, A=A, M=M,
                                     deficit=deficit, G=gradient(w), ratio=deficit/(M*gradient(w)),
                                     cft_prediction=prediction,
                                     cft_quotient=deficit/prediction if abs(prediction)>1e-13 else None,
                                     uniform_gap_in_sector=gap, residual=max(r0, rw)))
                    assert deficit >= -2e-9
    return rows


def xx_checks(cluster=False):
    rows = []
    for length in ([256, 512, 1024] if cluster else [32, 64, 128]):
        # For L divisible by four the spin-ring half-filled vacuum has AP fermions.
        for shape, mode in [('cosine', 1), ('cosine', 2), ('cosine', 4), ('gaussian', 2)]:
            for eps in [.02, .08]:
                w = profile(length, shape, eps, mode)
                mat = np.diag(w[:-1]/2, 1)+np.diag(w[:-1]/2, -1)
                mat[0, -1] = mat[-1, 0] = -w[-1]/2
                e = eigvalsh(mat)
                ew = e[:length//2].sum()
                e0 = -1/np.sin(np.pi/length)
                deficit = float(e0*w.mean()-ew)
                pred = cft_response(w, 1.)
                rows.append(dict(L=length, shape=shape, mode=mode, epsilon=eps,
                                 deficit=deficit, G=gradient(w), cft_prediction=pred,
                                 cft_quotient=deficit/pred if abs(pred)>1e-13 else None))
    return rows


def exact_counterexample():
    rows = []
    length = 4096
    for d in [2, 3, 4]:
        for delta in [1., .1, .01, .001]:
            eps = delta*.01
            for shape in ['cosine', 'gaussian']:
                w = profile(length, shape, eps, 2)
                bar = (w+np.roll(w, 1))/2
                dw = w-np.roll(w, 1)
                # Rationalized expression avoids subtracting nearly equal roots.
                deficit = float(np.sum(dw**2/(np.hypot(delta*bar, dw)+delta*bar)))
                rows.append(dict(d=d, L=length, delta=delta, t=1., epsilon=eps,
                                 shape=shape, M=delta, A=0., gap=2*delta,
                                 deficit=deficit, G=gradient(w), ratio=deficit/(delta*gradient(w)),
                                 limiting_ratio=1/(2*delta**2)))
    # Independent many-body verification of the analytic spectrum.
    checks = []
    for d, length in [(2, 4), (3, 3), (4, 3)]:
        delta, t = .2, .7
        a = np.diag([-delta]+[delta]*(d-1))
        x = np.zeros((d, d)); x[0, 1] = x[1, 0] = 1
        h = (np.kron(a, np.eye(d))+np.kron(np.eye(d), a))/2+t*(np.kron(x, np.eye(d))-np.kron(np.eye(d), x))
        w = profile(length, 'cosine', .1, 1)
        en, _, residual = ground(bond_operators(h, length), w)
        bar = (w+np.roll(w, 1))/2
        exact = -np.hypot(delta*bar, t*(w-np.roll(w, 1))).sum()
        checks.append(dict(d=d, L=length, energy_error=float(abs(en-exact)), residual=residual))
        assert abs(en-exact) < 1e-11
    return dict(scalars=rows, independent_ed=checks)


def additional_checks(cluster=False):
    rows = []
    for d in [3, 4]:
        h = np.eye(d*d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d*d, d*d)
        for length in ([12, 15] if d == 3 else [8, 12]) if cluster else ([6, 9] if d == 3 else [4, 8]):
            ops = bond_operators(h, length, balanced=True)
            e0, _, r0 = ground(ops, np.ones(length))
            for shape in ['cosine', 'gaussian']:
                w = profile(length, shape, .02, 2)
                G = gradient(w)
                if G < 1e-20:
                    continue
                ew, _, rw = ground(ops, w)
                deficit = e0*w.mean()-ew
                A = 2*np.pi*(d-1)/d
                pred = cft_response(w, A)
                rows.append(dict(model='SU_swap', d=d, L=length, shape=shape, mode=2,
                                 dimension=ops[0].shape[0], sector='balanced colors',
                                 epsilon=.02, A=A, M=1., deficit=deficit, G=G,
                                 ratio=deficit/G, cft_prediction=pred,
                                 cft_quotient=deficit/pred if abs(pred)>1e-13 else None,
                                 residual=max(r0, rw)))
    # Longitudinal divergence at the XXZ ferromagnetic point.
    z = np.diag([.5, -.5])
    for length in ([16, 20, 24] if cluster else [6, 8, 10, 12]):
        h = spin_bond(-1)+np.kron(z, np.eye(2))-np.kron(np.eye(2), z)
        ops = bond_operators(h, length, balanced=True)
        for eps in [.001, .002]:
            w = profile(length, 'cosine', eps, 1)
            ew, _, residual = ground(ops, w)
            deficit = -length*w.mean()/4-ew
            limit = length/(4*(length-1)*np.sin(np.pi/length)**2)
            rows.append(dict(model='XXZ_minus1_divergence', d=2, L=length,
                             shape='cosine', mode=1, epsilon=eps, t=1., M=.5,
                             deficit=deficit, G=gradient(w), ratio=deficit/(.5*gradient(w)),
                             second_order_ratio=limit, residual=residual,
                             sector='balanced; second-order minimizer over the ground multiplet'))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cluster', action='store_true')
    ap.add_argument('--suite', choices=['all', 'xx', 'xxz', 'random', 'exact', 'additional'], default='all')
    ap.add_argument('--output', default='small_results.json')
    args = ap.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Heavy calculations require a SLURM allocation.')
    start = time.perf_counter()
    out = dict(status='numerical', cluster=args.cluster, suite=args.suite)
    funcs = dict(xx=xx_checks, xxz=xxz_checks, random=random_checks, additional=additional_checks)
    for name, fn in funcs.items():
        if args.suite in ['all', name]:
            out[name] = fn(args.cluster)
            Path(args.output).write_text(json.dumps(out, indent=2)+'\n')
    if args.suite in ['all', 'exact']:
        out['exact'] = exact_counterexample()
    out['seconds'] = time.perf_counter()-start
    Path(args.output).write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({key: len(value) for key, value in out.items() if isinstance(value, list)}))
    print('Seconds:', out['seconds'])


if __name__ == '__main__':
    main()
