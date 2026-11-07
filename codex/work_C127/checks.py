"""C127 small local checks. Larger transfer matrices require --cluster under SLURM.
All generated files stay next to this script. No external project imports.
"""
from pathlib import Path
import argparse
import json
import os
import time
import numpy as np
from scipy.linalg import expm, eigh, eigvalsh
from scipy.integrate import quad
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent
sx = np.array([[0., 1.], [1., 0.]])/2
sy = np.array([[0., -1j], [1j, 0.]])/2
sz = np.diag([.5, -.5])

def bond(delta, reflected=False):
    signs = [-1, 1, -delta] if reflected else [1, 1, delta]
    return sum(c*np.kron(a, a) for c, a in zip(signs, [sx, sy, sz])).real

def boson_log(a, S):
    # Tail <= S exp[-a(N+1)]/[(1-exp(-a))(1-exp(-a(N+1)))].
    N = max(1, int(np.ceil((40 + max(0., -np.log(a)))/a)))
    vals = a*np.arange(1, N+1)
    ans = -S*np.log(-np.expm1(-vals)).sum()
    tail = S*np.exp(-a*(N+1))/(-np.expm1(-a))/(-np.expm1(-a*(N+1)))
    return float(ans), float(tail)

def boson_sectors(a, m, r=3):
    # r species of each momentum sign. Reference momenta reduced modulo m.
    N = max(1, int(np.ceil((40 + max(0., -np.log(a)))/a)))
    k = np.arange(1, N+1)
    weights = np.zeros(m); weights[0] = 1.
    scale = 0.
    # Positive cyclic convolution avoids Fourier cancellation in tiny sectors.
    for kk in k:
        z = np.exp(-a*kk)
        ell = m//np.gcd(m, int(kk))
        for sign in [-1, 1]:
            for _ in range(r):
                new = sum(z**j*np.roll(weights, sign*int(kk)*j) for j in range(ell))
                new /= -np.expm1(-a*kk*ell)
                norm = new.sum()
                weights = new/norm
                scale += np.log(norm)
    return [float(scale+np.log(c)) if c>0 else None for c in weights]

def xx_pressure(beta):
    if beta < 50:
        pc = 2/np.pi*quad(lambda k: np.logaddexp(beta*np.cos(k)/2, -beta*np.cos(k)/2),
                         0, np.pi/2, epsabs=2e-12)[0]
        gap = 2/np.pi*quad(lambda k: -np.log(np.tanh(beta*np.cos(k)/2)),
                          0, np.pi/2, epsabs=2e-12)[0]
    else:
        pc = beta/np.pi+2/(np.pi*beta)*quad(
            lambda y: np.log1p(np.exp(-y))/np.sqrt(1-(y/beta)**2), 0, 40, epsabs=1e-12)[0]
        gap = 2/(np.pi*beta)*quad(
            lambda y: -np.log(np.tanh(y/2))/np.sqrt(1-(y/beta)**2), 0, 40, epsabs=1e-12)[0]
    return pc, gap

def xx_logz(L, beta):
    # Physical spin trace with both Jordan-Wigner parity sectors.
    sectors = []
    for shift, parity in [(.5, 0), (0., 1)]:
        e = np.cos(2*np.longdouble(np.pi)*(np.arange(L, dtype=np.longdouble)+shift)/L)
        even, odd = np.longdouble(0), np.longdouble(-np.inf)
        for x in e:
            even, odd = np.logaddexp(even, odd-beta*x), np.logaddexp(odd, even-beta*x)
        sectors.append((even, odd)[parity])
    return float(np.logaddexp(*sectors))

def xx_cross(L, beta, N=20000):
    _, delta = xx_pressure(beta)
    ko = np.arcsinh(np.pi*(2*np.arange(N)+1)/beta)
    ke = np.arcsinh(2*np.pi*(np.arange(N)+1)/beta)
    ap = 4*np.log1p(np.exp(-L*ko)).sum()
    am = 4*np.log(-np.expm1(-L*ko)).sum()
    ao = np.logaddexp(ap, am)-np.log(2)
    ae = np.log(2)-L*delta+4*np.log1p(np.exp(-L*ke)).sum()
    ans = float(np.logaddexp(ao, ae))
    # Each sector's omitted log product is bounded using asinh u >= log(2u).
    tail = np.exp(np.log(4)+L*np.log(beta/(4*np.pi))+(1-L)*np.log(N-.5)-np.log(L-1))
    return ans, float(tail)

def rotation_bases(m):
    dim = 2**(2*m)
    paths = np.array(np.unravel_index(np.arange(dim), (2,)*(2*m))).T
    shifted = np.roll(paths, -2, axis=1)
    perm = np.ravel_multi_index(shifted.T, (2,)*(2*m))
    seen = set()
    orbits = []
    for a in range(dim):
        if a in seen: continue
        orb, b = [], a
        while b not in seen:
            orb.append(b); seen.add(b); b = int(perm[b])
        orbits.append(orb)
    bases = []
    for q in range(m):
        cols = []
        for orb in orbits:
            ell = len(orb)
            if (q*ell) % m: continue
            v = np.zeros(dim, complex)
            v[orb] = np.exp(-2j*np.pi*q*np.arange(ell)/m)/np.sqrt(ell)
            cols.append(v)
        bases.append(np.column_stack(cols))
    return paths, perm, bases

def transfer_spectrum(delta, beta, m, cache):
    paths, perm, bases = cache
    dim = len(paths)
    g = expm(-beta*bond(delta, True)/m).reshape(2, 2, 2, 2)
    mats = []
    for parity in [0, 1]:
        K = np.ones((dim, dim))
        for j in range(parity, 2*m, 2):
            jj = (j+1) % (2*m)
            K *= g[paths[:, jj, None], paths[None, :, jj], paths[:, j, None], paths[None, :, j]]
        mats.append(K)
    Ke, Ko = mats
    comm = max(np.linalg.norm(K[np.ix_(perm, perm)]-K)/np.linalg.norm(K) for K in mats)
    blocks = []
    mine = np.inf
    minu = np.inf
    for B in bases:
        e, W = eigh(B.conj().T@Ke@B, check_finite=False)
        mine = min(mine, float(e[0]/max(np.max(abs(e)), 1e-300)))
        X = (W*np.sqrt(np.maximum(e, 0)))@W.conj().T
        U = X@(B.conj().T@Ko@B)@X
        vals = eigvalsh(U, check_finite=False)
        blocks.append(vals)
    top = max(max(v) for v in blocks)
    gaps = []
    for vals in blocks:
        minu = min(minu, float(min(vals)/top))
        vals = vals[vals > 1e-13*top]
        gaps.append(-.5*np.log(vals/top))
    trace4 = float(sum(np.sum(v**2) for v in blocks))
    return gaps, top, dict(commutator=comm, min_Ke_relative=mine,
                          min_U_relative=minu, logZ4=float(np.log(trace4)),
                          retained=sum(map(len, gaps)), omitted=dim-sum(map(len, gaps)))

def physical_trotter(delta, beta, m, L=4):
    h = bond(delta)
    dim = 2**L
    Hs = [np.zeros((dim, dim)), np.zeros((dim, dim))]
    powers = 2**np.arange(L)
    for state in range(dim):
        digits = (state//powers) % 2
        for site in range(L):
            j = (site+1) % L
            a, b = int(digits[site]), int(digits[j])
            for aa in range(2):
                for bb in range(2):
                    target = state+(aa-a)*powers[site]+(bb-b)*powers[j]
                    Hs[site % 2][target, state] += h[2*aa+bb, 2*a+b]
    Q = expm(-beta*Hs[0]/m)@expm(-beta*Hs[1]/m)
    trot = np.log(np.trace(np.linalg.matrix_power(Q, m)))
    exact = logsumexp(-beta*eigvalsh(Hs[0]+Hs[1]))
    return float(trot), float(exact)

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster', action='store_true')
    args = parser.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Large transfer checks require SLURM')
    start = time.time()
    bconst = -np.log1p(-np.exp(-1))
    Kxx = np.log(3)+2*bconst
    out = dict(parameters=dict(S=6, V_over_M=4*np.e, K_candidate=float(np.log(3)),
                               XX_proved_S=2, XX_proved_K=float(Kxx)))
    exact = []
    for L in [2, 4, 6, 16, 64, 256]:
        for ratio in [.01, .1, 1., 4., 32.]:
            beta = L*ratio
            pc, _ = xx_pressure(beta)
            logz = xx_logz(L, beta)-L*pc
            cross, tail = xx_cross(L, beta)
            prod, bt = boson_log(2*np.pi*L/(2*np.e*beta), 2)
            sharp = np.log(3)+beta*(1/np.sin(np.pi/L)-L/np.pi)
            exact.append(dict(L=L, beta=beta, logZ=logz, product_value=cross,
                              product_tail=tail, product_error=logz-cross,
                              sharp_margin=float(sharp-logz), boson_margin=float(Kxx+prod-logz),
                              boson_tail=bt, zeroK_margin=float(prod-logz)))
    out['XX_exact'] = exact
    mats = []
    for beta in [1., 10., 100., 1000.]:
        for k in [1, 2, 10, 100, 1000]:
            w = 2*np.pi*k/beta
            rate = np.arcsinh(np.pi/beta)+np.arcsinh((2*k-1)*np.pi/beta)
            alpha = w/(2*np.e) if w<=2*np.e else np.log(w/2)
            mats.append(dict(beta=beta, k=k, omega=w, energy=float(rate),
                             low_frequency=w<2*np.e, LR_margin=float(rate-alpha),
                             global_linear_margin=float(rate-w/(2*np.e))))
    out['XX_matsubara'] = mats
    rank = []
    for beta in [10., 100., 1000.]:
        k = 1
        ns = np.arange(1, 501)
        nu = (2*ns+1)*np.pi/beta
        other = nu-2*np.pi*k/beta
        e = np.arcsinh(nu)+np.arcsinh(other)
        weights = np.exp(-e)/(beta*np.sqrt((1+nu**2)*(1+other**2)))
        rank.append(dict(beta=beta, k=k, distinct_count=len(np.unique(e)),
                         count_below_one=int(sum(e<=1)), first_8_energies=e[:8].tolist(),
                         first_8_reflected_weights=weights[:8].tolist()))
    out['XX_insertion_support'] = rank
    qtm = []
    mvals = [2, 3, 4, 5, 6] if args.cluster else [2, 3, 4]
    caches = {m: rotation_bases(m) for m in mvals}
    for delta in [0., .5, 1., 2.]:
        M = float(np.ptp(eigvalsh(bond(delta)))/2)
        for beta in [.25, 1., 2., 8.]:
            for m in mvals:
                gaps, top, checks = transfer_spectrum(delta, beta, m, caches[m])
                trot, phy = physical_trotter(delta, beta, m)
                checks.update(contraction_error=checks['logZ4']-trot, trotter_error=trot-phy)
                spec = []
                for q, e in enumerate(gaps):
                    qs = q if q<=m//2 else q-m
                    omega = 2*np.pi*abs(qs)/beta
                    spec.append(dict(q=q, q_signed=qs, omega=omega, energies=e.tolist(),
                                     min_gap=float(min(e)) if len(e) else None,
                                     finite_m_cone_margin=float(min(e)-omega/(4*np.e*M)) if len(e) and q else None))
                traces = []
                for L in [2, 4, 8, 16, 64]:
                    sec = np.array([logsumexp(-L*e) if len(e) else -np.inf for e in gaps])
                    z = float(logsumexp(sec))
                    a = 2*np.pi*L/(4*np.e*M*beta)
                    ref, tail = boson_log(a, 6)
                    bs = boson_sectors(a, m)
                    sector_margins = [float(np.log(3)+b-s) if b is not None and np.isfinite(s) else None
                                      for b, s in zip(bs, sec)]
                    # Discarded U eigenvalues <= 1e-13 top contribute at most dim*1e-13^(L/2).
                    omitted_trace_bound = checks['omitted']*1e-13**(L/2)
                    traces.append(dict(L=L, logZ=z, zeroK_margin=ref-z,
                                       candidate_margin=float(np.log(3)+ref-z), boson_tail=tail,
                                       sector_logZ=sec.tolist(), sector_candidate_margins=sector_margins,
                                       discarded_trace_bound=omitted_trace_bound))
                qtm.append(dict(delta=delta, M=M, beta=beta, m=m, checks=checks,
                                spectrum=spec, traces=traces))
    out['QTM'] = qtm
    analytic = []
    for a in np.geomspace(.001, 100, 101):
        p, tail = boson_log(a, 2)
        t = 2*np.pi/a
        analytic.append(dict(a=float(a), product_upper_margin=float(np.pi**2/(3*a)-p),
                             product_lower_margin=float(p-(bconst*t/np.pi-2*bconst)), tail=tail))
    out['boson_inequalities'] = analytic
    theta = []
    for a in np.geomspace(.001, 1000, 101):
        N = int(np.ceil(np.sqrt(100*a)))+1
        logtheta = float(np.log1p(2*np.exp(-np.arange(1,N+1)**2/(2*a)).sum()))
        theta.append(dict(a=float(a), logtheta=logtheta, margin=float(4*a-logtheta)))
    out['theta_inequality'] = theta
    # Gibbs variational comparison in a noncommuting two-dimensional example.
    H0 = np.diag([0., 1.])
    H = np.array([[.3, .4], [.4, .2]])
    gp = []
    for L in [.1, 1., 10.]:
        A, B = expm(-L*H), expm(-L*H0)
        z, z0 = np.trace(A), np.trace(B)
        diff = np.log(z/z0)
        low = -L*np.trace((B/z0)@(H-H0))
        high = -L*np.trace((A/z)@(H-H0))
        gp.append(dict(L=L, lower_margin=float(diff-low), upper_margin=float(high-diff)))
    out['Gibbs_comparison'] = gp
    potts = []
    for d in [2, 3, 5]:
        for beta in [4., 8., 16., 24.]:
            L = 2*max(1, int(np.exp(beta/2)/2))
            ratio = -np.expm1(-beta)/(1+(d-1)*np.exp(-beta))
            logz = np.log1p((d-1)*np.exp(L*np.log(ratio)))
            prod, _ = boson_log(2*np.pi*L/(2*np.e*beta), 2*(d*d-1))
            potts.append(dict(d=d, beta=beta, L=L, logZ=float(logz),
                              logd=float(np.log(d)), zeroK_margin=float(prod-logz)))
    out['classical_K_necessity'] = potts
    nonlinear = []
    for L in [2, 4, 8, 16]:
        for beta in [.1, 1., 10., 100.]:
            V = 2*np.e
            N = 20000
            w = 2*np.pi*np.arange(1, N+1)/beta
            alpha = np.where(w<=V, w/V, np.log(w/2))
            value = -6*np.log(-np.expm1(-L*alpha)).sum()
            tail = 6*(beta/np.pi)**L*N**(1-L)/(L-1)/(-np.expm1(-L*alpha[-1]))
            upper = 6*beta*V/(2*np.pi*L)*(np.pi**2/6+2/np.expm1(2))
            nonlinear.append(dict(L=L, beta=beta, value=float(value), tail=float(tail),
                                  certified_margin=float(upper-value-tail)))
    out['nonlinear_product'] = nonlinear
    out['runtime_seconds'] = time.time()-start
    dest = ROOT/('checks_cluster.json' if args.cluster else 'checks.json')
    dest.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    summary = dict(runtime=out['runtime_seconds'], qtm_cases=len(qtm),
                   max_contraction_error=max(abs(x['checks']['contraction_error']) for x in qtm),
                   min_qtm_candidate_margin=min(t['candidate_margin'] for x in qtm for t in x['traces']),
                   min_XX_sharp_margin=min(x['sharp_margin'] for x in exact),
                   min_XX_boson_margin=min(x['boson_margin'] for x in exact),
                   max_XX_product_error=max(abs(x['product_error']) for x in exact),
                   rank_examples=rank)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    run()
