"""Literal vacuum length words and tail-sensitive finite-volume response.

One BLAS thread, no thermodynamic extrapolation, no finite-gap bound claimed
uniform in volume. Writes only the adjacent small JSON report.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

from itertools import product
from pathlib import Path
import importlib.util
import json
import time

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'clock_checks', HERE.parent/'mobile_vacuum_clock_20261117'/'checks.py')
clock = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clock)
mobile = clock.mobile
mp.mp.dps = 70


def cell(k, q=.25):
    beta = q*q
    Q = mp.matrix(k+1)
    for j in range(k):
        Q[j, j] += mp.mpf(q)**2
        Q[j+1, j+1] += 1
        Q[j, j+1] = Q[j+1, j] = -mp.mpf(q)
    Q[k, k] += mp.mpf(q)**2
    es, vs = mp.eigsy(Q)
    v = np.array([float(vs[j, 0]) for j in range(k+1)])
    v *= np.sign(v[0])
    b = [np.zeros((k+1, k+1)) for _ in range(k+1)]
    for j in range(k):
        e = np.zeros(k+1); e[j:j+2] = [q, -1]
        b[j] += np.outer(e, e)
    b[k-1][-1, -1] += beta/2
    b[k][-1, -1] += beta/2
    a = np.array([v@x@v for x in b])
    xi = np.stack([x@v-mean*v for x, mean in zip(b, a)])
    return dict(E=float(es[0]), v=v, b=b, a=a, xi=xi,
                p=float(v[-1]**2), tau=beta**(k+1)/sum(beta**j for j in range(k+1)))


def literal_words(ks):
    words = []
    for heads in product(*(range(k+1) for k in ks)):
        word = []
        for k, head in zip(ks, heads):
            word += ([('L', 0)]*k+[('m',)] if head == k else
                     [('L', 0)]*head+[('v',)]+[('R', 0)]*(k-head-1)+[('e', 0)])
        words.append(tuple(mobile.frozen.INDEX[x] for x in word))
    return words


def literal_audit():
    out = []
    for n in (1, 2):
        p = mobile.frozen.config(n, field=0., lam=0.)
        mu = p['mu']
        for ks in ([n]*(n+2), [n+1]*(n+1)):
            words = literal_words(ks)
            actual = mobile.material_sources(words, mobile.frozen.literal_bond(p)[0])
            expected = clock.sources(ks, .25, mu)
            error = max(float(sp.linalg.norm(a-b)) for a, b in zip(actual, expected))
            assert error < 1e-12
            out.append(dict(n=n, word=ks, occupied_length=sum(k+1 for k in ks),
                            material_dimension=len(words), source_error=error))
        en, ep = cell(n)['E'], cell(n+1)['E']
        diff = (n+1)*(mu+ep)-(n+2)*(mu+en)
        assert abs(diff-(n+1)*ep) < 1e-14
    return out


def response(N, V):
    gaps, H, contacts = clock.bath(N, V)
    energies, vectors = la.eigh(H.toarray().real)
    omega = vectors[:, 0]
    alpha = mobile.alpha(N, V)
    X = np.stack([a@omega-alpha*omega for a in contacts], axis=1)
    Z = vectors[:, 1:].T@X
    G = (Z.T/(energies[1:]-energies[0]))@Z
    assert abs(energies[0]-mobile.free_energy(N, V)) < 1e-11
    assert np.max(np.abs(omega@X)) < 1e-12
    assert la.eigvalsh(G)[0] > -1e-12
    assert np.max(np.abs(G-np.roll(np.roll(G, 1, axis=0), 1, axis=1))) < 1e-11
    return G, dict(N=N, V=V, dimension=len(gaps), alpha=alpha,
                   finite_bath_gap=float(energies[1]-energies[0]),
                   response_norm=float(la.norm(G, 2)))


def word_second(k, count, mu, G):
    info = cell(k)
    source = info['a'].copy(); source[-2:] += mu/2
    s = np.tile(source, count)
    R = float(np.sum(G[:k+1, :k+1]*(info['xi']@info['xi'].T)))
    return dict(scalar=-float(s@G@s), clock=-count*R, single_clock_R=R)


def coefficient(n, V):
    N = (n+1)*(n+2)
    G, meta = response(N, V)
    mu = -(n+2)*cell(n)['E']
    a = word_second(n, n+2, mu, G)
    b = word_second(n+1, n+1, mu, G)
    parts = {key: b[key]-a[key] for key in ('scalar', 'clock')}
    second = sum(parts.values())
    margin = (n+1)*cell(n+1)['E']
    first = meta['alpha']*margin
    ratio = b['single_clock_R']/a['single_clock_R']
    assert ratio <= 32/81+1e-10
    assert parts['clock'] > 0
    # This sign is an observation about the actual free bath, not assumed.
    assert parts['scalar'] < 0 and second > 0
    mu_cross = mu+margin
    ac = word_second(n, n+2, mu_cross, G)
    bc = word_second(n+1, n+1, mu_cross, G)
    cross_second = sum(bc[key]-ac[key] for key in ('scalar', 'clock'))
    return dict(n=n, mu=mu, frozen_margin=margin, first_coefficient=first,
                A_second=a, B_second=b, second_difference=parts,
                second_coefficient=second, relative_second_per_J=second/first,
                clock_response_ratio=ratio,
                frozen_mu_cross=mu_cross,
                mobile_mu_cross_linear_shift=cross_second/meta['alpha'], **meta)


def form_factor_checks():
    reports = []
    theta = np.linspace(.02, 2*np.pi-.02, 401)
    for q in (.25, .125):
        beta = q*q
        previous = None
        for k in range(1, 10):
            info = cell(k, q)
            phases = np.exp(1j*np.outer(theta, np.arange(k+1)))
            form = np.sum(abs(phases@info['xi'])**2, axis=1)
            scale = beta**2*info['p']*abs(np.exp(1j*theta)-1)**2
            lower, upper = (7/16)**2*scale, (35/36)**2*scale
            assert np.min(form/lower) > 1-1e-7
            assert np.max(form/upper) < 1+1e-7
            assert 25/32*q**(2*k) <= info['p'] <= q**(2*k)
            ratio = None if previous is None else float(np.max(form/previous))
            if ratio is not None:
                assert ratio <= (512/81)*q*q*(1+1e-7)
            reports.append(dict(k=k, q=q, E=info['E'], endpoint_weight=info['p'],
                                max_successive_form_ratio=ratio,
                                max_upper_fraction=float(np.max(form/upper)),
                                min_lower_fraction=float(np.min(form/lower)),
                                zero_mode_residual=float(la.norm(info['xi'].sum(axis=0)))))
            previous = form
    return reports


def tilted_candidates(N, count, mu, J, upper, s=.6):
    # Positive clocks and C_m<=1: E_word >= E_free+J*mu*count.
    vmax = max(0, int(np.ceil((upper+2*N-J*mu*count)/s)))+1
    return [v for v in range(vmax+1)
            if mobile.free_energy(N, v)+s*v+J*mu*count <= upper+1e-11]


def finite_scan(coefficients):
    records = []
    for n, V, Js in ((1, 6, (1/320, 1/160, 1/80, .05, .2, 1.)),
                     (2, 2, (1/80, .05, .2)), (2, 4, (1/80,))):
        N = (n+1)*(n+2)
        info = next(x for x in coefficients if x['n'] == n and x['V'] == V)
        mu = info['mu']
        words = ([n]*(n+2), [n+1]*(n+1))
        sources = [clock.sources(ks, .25, mu) for ks in words]
        for J in Js:
            energies, residuals, dimensions, candidates = [], [], [], []
            for ks, hs in zip(words, sources):
                H = clock.assemble(hs, V, J)[1]
                energy, psi, residual = clock.minimum(H)
                energies.append(energy); residuals.append(residual); dimensions.append(H.shape[0])
                candidates.append(tilted_candidates(N, len(ks), mu, J, energy+.6*V))
            difference = energies[1]-energies[0]
            first = J*info['first_coefficient']
            second = J*J*info['second_coefficient']
            assert difference > 0
            certified = candidates == [[V], [V]]
            if n == 1 and J <= 1/80:
                assert certified
            records.append(dict(n=n, N=N, V=V, rho=N/(N+V), J_over_t=J,
                                energies=energies, difference=difference,
                                first_order_margin=first, ratio=difference/first,
                                second_order_prediction=first+second,
                                after_second_remainder=difference-first-second,
                                dimensions=dimensions, eigenpair_residual=max(residuals),
                                s=.6, both_tilted_minima_certified=certified,
                                possible_tilted_V=candidates))
    return records


def crossing_checks(coefficients):
    n, V, N = 1, 6, 6
    info = next(x for x in coefficients if x['n'] == n and x['V'] == V)
    result = []
    for J in (1/80, .05):
        def difference(mu):
            energies = [clock.minimum(clock.assemble(clock.sources(ks, .25, mu), V, J)[1])[0]
                        for ks in ([1, 1, 1], [2, 2])]
            return energies[1]-energies[0]
        center = info['frozen_mu_cross']
        root = brentq(difference, center-2e-4, center+2e-4, xtol=2e-10)
        candidates = []
        for ks in ([1, 1, 1], [2, 2]):
            e = clock.minimum(clock.assemble(clock.sources(ks, .25, root), V, J)[1])[0]
            candidates.append(tilted_candidates(N, len(ks), root, J, e+.6*V))
        assert root > center > info['mu']
        result.append(dict(n=n, N=N, V=V, J_over_t=J, frozen_mu=info['mu'],
                           frozen_crossing=center, observed_crossing=root,
                           linear_prediction=center+J*info['mobile_mu_cross_linear_shift'],
                           crossing_residual=abs(difference(root)),
                           possible_tilted_V=candidates))
    return result


def physical_checks():
    p = mobile.frozen.config(1, field=0., lam=0.)
    _, bond = mobile.literal_mobile(p, J=1/80, t=1., mu=0.)
    results = []
    for ks in ([1, 1, 1], [2, 2]):
        fixed = literal_words(ks)
        orbit = sorted({w[j:]+w[:j] for w in fixed for j in range(6)})
        for V in (1, 2):
            words = mobile.physical_basis(orbit, V)
            H = mobile.physical_matrix(words, bond)
            ep, _, rp = clock.minimum(H)
            er, _, rr = clock.minimum(clock.assemble(clock.sources(ks, .25, p['mu']), V, 1/80)[1])
            assert abs(ep-er) < 1e-10
            results.append(dict(word=ks, N=6, V=V, physical_dimension=len(words),
                                physical_minus_rooted=ep-er, residual=max(rp, rr)))
    return results


def scalar_sandwich():
    n, N, V, J, s, q = 1, 6, 6, 1/80, .6, .25
    mu = -(n+2)*cell(n)['E']
    gaps, H, contacts = clock.bath(N, V)
    results = []
    for ks in ([1, 1, 1], [2, 2]):
        def scalar(values):
            source = sp.csr_matrix(H.shape)
            offset = 0
            for k, value in zip(ks, values):
                marker = offset+k
                source += value*(contacts[marker-1]+contacts[marker])/2
                offset += k+1
            return clock.minimum(H+J*source)[0]+s*V
        theta = .5
        tau = np.array([cell(k)['tau'] for k in ks])
        delta = np.array([(1+q*q-2*q*np.cos(np.pi/(k+1)))/257 for k in ks])
        kappa = theta*delta/(theta*delta+q*q)
        f0 = scalar(np.full(len(ks), mu))
        fmiddle = scalar(mu+kappa*tau/(1-theta))
        lower = theta*f0+(1-theta)*fmiddle
        upper = scalar(mu+tau)
        frozen_source = np.concatenate([cell(k)['a']+np.r_[np.zeros(k-1), mu/2, mu/2]
                                        for k in ks])
        frozen_upper = clock.minimum(H+J*sum((a*v for a, v in zip(contacts, frozen_source)),
                                             sp.csr_matrix(H.shape)))[0]+s*V
        actual = clock.minimum(clock.assemble(clock.sources(ks, q, mu), V, J)[1])[0]+s*V
        # All three scalar energies and the full energy use the actual tilted minimizer.
        for value in (f0, fmiddle, upper, frozen_upper, actual):
            assert tilted_candidates(N, len(ks), mu, J, value, s) == [V]
        assert lower <= actual+1e-11 <= upper+2e-11
        assert actual <= frozen_upper+1e-11
        results.append(dict(word=ks, N=N, V=V, J_over_t=J, lower=lower,
                            exact=actual, upper=upper, frozen_upper=frozen_upper,
                            interval_width=upper-lower))
    assert results[1]['lower'] > results[0]['frozen_upper']
    return results


def main():
    start = time.time()
    out = dict(status='finite numerical checks; no uniform mobile length-selection theorem',
               versions=dict(numpy=np.__version__, scipy=scipy.__version__, mpmath=mp.__version__),
               BLAS_threads=1)
    jobs = [('literal_audit', literal_audit), ('form_factor_bounds', form_factor_checks),
            ('second_order', lambda: [coefficient(n, V) for n, V in
                                     ((1, 4), (1, 6), (1, 8), (2, 2), (2, 4), (3, 2))]),
            ('finite_scan', lambda: finite_scan(out['second_order'])),
            ('crossings', lambda: crossing_checks(out['second_order'])),
            ('physical_gluing', physical_checks), ('scalar_sandwich', scalar_sandwich)]
    for name, run in jobs:
        out[name] = run()
        print(name, 'passed', flush=True)
        (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    out['runtime_seconds'] = time.time()-start
    (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print('all checks passed in', round(out['runtime_seconds'], 2), 'seconds')


if __name__ == '__main__':
    main()
