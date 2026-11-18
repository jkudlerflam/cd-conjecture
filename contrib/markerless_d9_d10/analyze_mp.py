"""High-precision band analysis of the 9-letter logical chain."""
import mpmath as mp, itertools, json, sys, time
from scalar9 import blocks

def roots_mp(mu, T):
    """2n roots z of det(T z^2 + mu z + T^dag) = 0"""
    n = T.rows
    Ti = T**-1
    C = mp.zeros(2*n, 2*n)
    for i in range(n): C[i, n+i] = 1
    A = -(Ti*T.H); B = -(Ti*mu)
    for i in range(n):
        for j in range(n): C[n+i, j] = A[i,j]; C[n+i, n+j] = B[i,j]
    ev = mp.eig(C, left=False, right=False)
    return list(ev)

def analyse(n, a, D, eta, fgam=0, phiA=None, psiA=None, psiB=None, verbose=True):
    """eta = t_d/beta with beta the express-rail amplitude; fgam = uniform field in units of the smallest
    eigenvalue of the Hermitian part of the hopping."""
    a = mp.mpf(a); D = mp.mpf(D)
    # beta: amplitude from X_n to X_0 through the gapped rail ~ 1/det of the express block
    d0, d1 = mp.mpf(1), a + D
    for k in range(2, n+2): d0, d1 = d1, (a+D)*d1 - d0
    beta = 1/d1
    td = mp.mpf(eta)*beta
    T, muL, muR = blocks(n, a=a, D=D, td=td, phiA=phiA, psiA=psiA, psiB=psiB)
    mu = muL + muR
    herm_mu = max(abs(mu[i,j]-mp.conj(mu[j,i])) for i in range(n) for j in range(n))
    th = mp.arg(sum(T[i,i] for i in range(n)))
    Th = T*mp.expj(-th)
    S = (Th + Th.H)/2; Aa = (Th - Th.H)/(2j)
    es, _ = mp.eigh(S)
    scT = max(abs(T[i,j]) for i in range(n) for j in range(n))
    out = dict(n=n, a=float(a), D=float(D), eta=float(eta), log10beta=float(mp.log10(beta)),
               log10T=float(mp.log10(scT)), eigS_over_T=[float(e/scT) for e in es],
               antiherm=float(max(abs(Aa[i,j]) for i in range(n) for j in range(n))/scT),
               mu_over_T=float(max(abs(mu[i,j]) for i in range(n) for j in range(n))/scT), herm_mu=float(herm_mu/scT))
    out['definite'] = bool(es[0] > 0)
    f = mp.mpf(fgam)*abs(es[0])
    M = mu.copy()
    for i in range(n): M[i,i] += f
    zs = roots_mp(M, T)
    off = max(abs(abs(z)-1) for z in zs)
    ks = sorted(float(mp.arg(z)) for z in zs if abs(abs(z)-1) < mp.mpf(10)**(-mp.mp.dps//3))
    out['on_circle'] = len(ks); out['max_off_circle'] = float(off)
    if len(ks) == 2*n:
        twopi = 2*float(mp.pi); m = 2*n
        d2 = min((ks[(i+1)%m]-ks[i])%twopi for i in range(m))
        sums = sorted((ks[i]+ks[j])%twopi for i in range(m) for j in range(i, m))
        d4 = min((sums[(i+1)%len(sums)]-sums[i])%twopi for i in range(len(sums)))
        out.update(d2=d2, d4=d4, ks=ks, npairs=len(sums))
    return out

if __name__ == '__main__':
    mp.mp.dps = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    n, a, D = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    etas = [mp.mpf(x) for x in sys.argv[4].split(',')]
    for eta in etas:
        for fg in (0, 0.3):
            t0 = time.time(); r = analyse(n, a, D, eta, fg)
            ks = r.pop('ks', None)
            print(json.dumps({k:(v if not isinstance(v,list) else ['%.2e'%x for x in v]) for k,v in r.items()}), 'fgam=%s'%fg, '%.1fs'%(time.time()-t0), flush=True)
