"""One explicit member of the 9-letter family: logical Bloch data in arbitrary precision, its 2n Fermi points,
separations, whitened slopes and the error allowance they imply, and the cell-length selector."""
import mpmath as mp, json, sys
from scalar9 import blocks, Hs
from analyze_mp import roots_mp

def selector(b):
    a = lambda x: 2*(1-mp.cos(mp.pi/(x+1)))
    delta = lambda x: (x+1)*a(x) - x*a(x+1)
    nu = 2 - (delta(b) + delta(b-1))/2
    e = lambda x: (nu - 2 + a(x))/x
    gap = min(e(x) - e(b) for x in range(1, 40*b) if x != b)
    return nu, e(b), gap, 3*0 + (2*mp.cos(mp.pi/(b+1)) - 2*mp.cos(2*mp.pi/(b+1)))

def member(n, a, D, eta, fgam, dps):
    mp.mp.dps = dps
    a = mp.mpf(a); D = mp.mpf(D); eta = mp.mpf(eta); fgam = mp.mpf(fgam)
    d0, d1 = mp.mpf(1), a + D
    for k in range(2, n+2): d0, d1 = d1, (a+D)*d1 - d0
    beta = 1/d1; td = eta*beta
    T, muL, muR = blocks(n, a=a, D=D, td=td)
    token_gap = min(mp.eigh(Hs(n, a, D, ph, td, ps))[0][0] for ph, ps in [(mp.pi/(2*(n+1)), mp.pi), (0, mp.pi/2)])
    assert token_gap > 0, 'token sector is not gapped above the cell vacuum'
    mu = muL + muR
    th = mp.arg(sum(T[i,i] for i in range(n)))
    Th = T*mp.expj(-th); S = (Th + Th.H)/2
    es, U = mp.eigh(S)
    assert es[0] > 0, 'Hermitian part of the hopping is not positive'
    f = fgam*es[0]
    M = mu.copy()
    for i in range(n): M[i,i] += f
    zs = roots_mp(M, T)
    off = max(abs(abs(z)-1) for z in zs)
    ks = sorted(mp.arg(z) for z in zs)
    Wm = U*mp.diag([1/mp.sqrt(e) for e in es])*U.H
    slopes, resid = [], []
    for k in ks:
        z = mp.expj(k)
        h = Wm*(M + T*z + T.H/z)*Wm; hp = Wm*(1j*T*z - 1j*T.H/z)*Wm
        ev, V = mp.eigh((h + h.H)/2)
        i0 = min(range(n), key=lambda i: abs(ev[i])); u = V[:, i0]
        slopes.append((u.H*hp*u)[0].real); resid.append(abs(ev[i0]))
    twopi = 2*mp.pi; m = 2*n
    d2 = min((ks[(i+1)%m]-ks[i]) % twopi for i in range(m))
    sums = sorted((ks[i]+ks[j]) % twopi for i in range(m) for j in range(i, m))
    d4 = min((sums[(i+1)%len(sums)]-sums[i]) % twopi for i in range(len(sums)))
    minslope = min(abs(s) for s in slopes)
    nu, eb, gap, clockgap = selector(n+1)
    scT = max(abs(T[i,j]) for i in range(n) for j in range(n))
    return dict(n=n, cell_length=n+1, a=str(a), D=str(D), eta=str(eta), fgam=str(fgam), dps=dps,
        Delta=mp.nstr(a - 2*mp.cos(mp.pi/(n+2)), 20), phiA_over_pi=mp.nstr(mp.mpf(1)/(2*(n+1)), 20), psiA_over_pi='1', psiB_over_pi='0.5',
        token_gap_above_vacuum=mp.nstr(token_gap, 12), alpha4=mp.nstr((mp.mpf(2)/(n+2))**2*mp.sin(mp.pi/(n+2))**4, 15),
        log10_gamma_min=mp.nstr(mp.log10(es[0]), 12), log10_field=mp.nstr(mp.log10(f), 12) if f > 0 else None, log10_beta=mp.nstr(mp.log10(beta), 12), log10_td=mp.nstr(mp.log10(td), 12),
        phase_T_over_pi=mp.nstr(th/mp.pi, 15), log10_maxT=mp.nstr(mp.log10(scT), 12),
        hermitian_hopping_eigs_over_maxT=[mp.nstr(e/scT, 8) for e in es],
        whitened_antihermitian_norm=mp.nstr(mp.norm(Wm*((Th - Th.H)/2)*Wm, 2), 8),
        whitened_mu_norm=mp.nstr(mp.norm(Wm*mu*Wm, 2), 8),
        fermi_points=2*n if off < mp.mpf(10)**(-dps//3) else 'NOT ALL ON CIRCLE', max_dist_from_unit_circle=mp.nstr(off, 5),
        roots=[mp.nstr(k, 25) for k in ks], right_movers=sum(1 for s in slopes if s > 0), left_movers=sum(1 for s in slopes if s < 0),
        root_separation=mp.nstr(d2, 12), pair_sum_separation=mp.nstr(d4, 12), unordered_pairs=len(sums),
        min_whitened_slope=mp.nstr(minslope, 8), max_whitened_root_residual=mp.nstr(max(resid), 5),
        whitened_error_allowance=mp.nstr(d4*minslope/4, 8),
        selector=dict(nu=mp.nstr(nu, 20), e_b=mp.nstr(eb, 12), density_gap=mp.nstr(gap, 10), clock_gap=mp.nstr(clockgap, 10)))
if __name__ == '__main__':
    n, a, D, eta, fgam, dps = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])
    r = member(n, a, D, eta, fgam, dps)
    json.dump(r, open(f'member_n{n}_eta{eta}_f{fgam}_dps{dps}.json', 'w'), indent=1)
    print(json.dumps({k:v for k,v in r.items() if k not in ('roots','hermitian_hopping_eigs_over_maxT')}, indent=1))
