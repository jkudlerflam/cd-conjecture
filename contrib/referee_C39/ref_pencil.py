"""Referee high-precision check of the C39 limiting reference, n = 1..8.

Independent of certify.py, pencil.py and the Lyapunov route: the Gram matrix
Gamma(w) of PROOFS eq. (7) is evaluated from the closed eigenmode formula

  Gamma_ij = sum alpha^(i)_{m1 m1'} alpha^(j)_{m2 m2'} J(l_m1, l_m1', l_m2, l_m2'),
  J = (l1+l1'+l2+l2') / [(l1+l2)(l1'+l2')(l1'+l2)(l1+l2')],

with l_m = a - 2w cos(m pi/(n+1)) the eigenvalues of A(w)=aI-w adj(P_n).
All arithmetic is decimal at a user-set precision. Roots of det h(k) are
bracketed by exact sign changes of a degree-2n polynomial (2n brackets for
2n roots, so each bracket holds exactly one root) and then bisected.
"""
import sys, json, time
from decimal import Decimal as Dm, getcontext
from fractions import Fraction as Fr
from math import factorial
from itertools import combinations_with_replacement

def setprec(p):
    getcontext().prec = p

def D(x):
    if isinstance(x, Fr):
        return Dm(x.numerator) / Dm(x.denominator)
    return Dm(x)

_pi_cache = {}
def pi():
    p = getcontext().prec
    if p not in _pi_cache:
        _pi_cache[p] = 16 * atan_series(Dm(1) / 5) - 4 * atan_series(Dm(1) / 239)
    return _pi_cache[p]

def atan_series(x):
    # valid |x| < 1, used for small |x|
    eps = Dm(10) ** (-getcontext().prec - 5)
    s = Dm(0); term = x; k = 0; x2 = x * x
    while True:
        add = term / (2 * k + 1)
        s += add if k % 2 == 0 else -add
        if abs(add) < eps: break
        term *= x2; k += 1
    return s

def atan(x):
    if abs(x) <= Dm('0.4'):
        return atan_series(x)
    sgn = 1 if x > 0 else -1
    y = (x - sgn) / (1 + sgn * x)
    return sgn * pi() / 4 + atan_series(y)

def sincos(x):
    eps = Dm(10) ** (-getcontext().prec - 5)
    # reduce
    twopi = 2 * pi()
    x = x - twopi * int(x / twopi)
    s = Dm(0); c = Dm(0); term = Dm(1); k = 0
    while True:
        if k % 4 == 0: c += term
        elif k % 4 == 1: s += term
        elif k % 4 == 2: c -= term
        else: s -= term
        k += 1; term = term * x / k
        if abs(term) < eps and k > 4: break
    return s, c

def tan(x):
    s, c = sincos(x); return s / c

def det(M):
    M = [row[:] for row in M]; n = len(M); d = Dm(1)
    for j in range(n):
        p = max(range(j, n), key=lambda i: abs(M[i][j]))
        if M[p][j] == 0: return Dm(0)
        if p != j: M[p], M[j] = M[j], M[p]; d = -d
        piv = M[j][j]; d *= piv
        for i in range(j + 1, n):
            f = M[i][j] / piv
            if f:
                Mi = M[i]; Mj = M[j]
                for k in range(j + 1, n): Mi[k] -= f * Mj[k]
    return d

def gram(n, w, a=4):
    """Gamma(w) by the eigenmode formula, Decimal."""
    w = D(w); a = Dm(a); P = pi()
    lam = []; phi = []
    for m in range(1, n + 1):
        s, c = sincos(m * P / (n + 1)); lam.append(a - 2 * w * c)
        row = []
        for j in range(n):
            sj, _ = sincos(m * (j + 1) * P / (n + 1)); row.append(sj * (Dm(2) / (n + 1)).sqrt())
        phi.append(row)
    # alpha^(j)_{m m'} = phi_m(0) phi_m(j) phi_m'(j) phi_m'(0)
    al = [[[phi[m][0] * phi[m][j] * phi[mp][j] * phi[mp][0] for mp in range(n)] for m in range(n)] for j in range(n)]
    J = {}
    for m1 in range(n):
        for m1p in range(n):
            for m2 in range(n):
                for m2p in range(n):
                    l1, l1p, l2, l2p = lam[m1], lam[m1p], lam[m2], lam[m2p]
                    J[m1, m1p, m2, m2p] = (l1 + l1p + l2 + l2p) / ((l1 + l2) * (l1p + l2p) * (l1p + l2) * (l1 + l2p))
    # contract
    G = [[Dm(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            s = Dm(0)
            for (m1, m1p, m2, m2p), val in J.items():
                s += al[i][m1][m1p] * al[j][m2][m2p] * val
            G[i][j] = G[j][i] = s
    return G

def C0(n, a=4):
    return [[Fr(factorial(2*i+2*j+2), factorial(2*i+1)*factorial(2*j+1)*(2*a)**(2*i+2*j+3)) for j in range(n)] for i in range(n)]

def cholesky(A):
    n = len(A); L = [[Dm(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = s.sqrt() if i == j else s / L[j][j]
    return L

def tri_inv_lower(L):
    n = len(L); X = [[Dm(0)] * n for _ in range(n)]
    for c in range(n):
        for i in range(n):
            s = (Dm(1) if i == c else Dm(0)) - sum(L[i][k] * X[k][c] for k in range(i))
            X[i][c] = s / L[i][i]
    return X

def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def transpose(A):
    return [list(r) for r in zip(*A)]

def jacobi(M, tol_exp=None):
    n = len(M); A = [r[:] for r in M]; V = [[Dm(int(i == j)) for j in range(n)] for i in range(n)]
    tol = Dm(10) ** (-(getcontext().prec - 20))
    for sweep in range(200):
        off = sum((A[i][j] ** 2 for i in range(n) for j in range(n) if i != j), Dm(0))
        if off.sqrt() < tol * Dm(10) ** (-0): break
        for p in range(n):
            for q in range(p + 1, n):
                if A[p][q] == 0: continue
                # skip negligible relative to both diagonals (graded matrices)
                if abs(A[p][q]) < tol * (abs(A[p][p] * A[q][q])).sqrt(): A[p][q] = A[q][p] = Dm(0); continue
                th = (A[q][q] - A[p][p]) / (2 * A[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + (th * th + 1).sqrt())
                c = 1 / (t * t + 1).sqrt(); s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq
    return [A[i][i] for i in range(n)], V

def reference(n, w, v, u, eta, a=4):
    """Return A=C(w), B=D_v C(wv) D_v, E=diag(u^4j), r_j, x_j (x_jj=1), theta_j, s_j."""
    w, v, u, eta = D(w), D(v), D(u), D(eta)
    Gw = gram(n, w, a); Gwv = gram(n, w * v, a)
    A = [[Gw[i][j] / w ** (2 * i + 2 * j) for j in range(n)] for i in range(n)]
    B = [[Gwv[i][j] / w ** (2 * i + 2 * j) for j in range(n)] for i in range(n)]
    E = [u ** (4 * j) for j in range(n)]
    L = cholesky(A); Li = tri_inv_lower(L)
    M = matmul(matmul(Li, B), transpose(Li))
    ev, V = jacobi(M)
    order = sorted(range(n), key=lambda i: -ev[i])
    r = [ev[i] for i in order]
    LiT = transpose(Li)
    X = []
    for idx, i in enumerate(order):
        y = [V[k][i] for k in range(n)]
        x = [sum(LiT[p][k] * y[k] for k in range(n)) for p in range(n)]
        x = [z / x[idx] for z in x]
        X.append(x)
    theta = [-atan(eta * rj) for rj in r]
    s = []
    for j in range(n):
        x = X[j]
        num = sum(E[i] * x[i] ** 2 for i in range(n))
        den = 2 * (1 + (eta * r[j]) ** 2).sqrt() * sum(x[i] * A[i][k] * x[k] for i in range(n) for k in range(n))
        s.append(num / den)
    return dict(A=A, B=B, E=E, r=r, X=X, theta=theta, s=s)

def Pdet(ref, eta, mu, x):
    A, B, E = ref['A'], ref['B'], ref['E']; n = len(A); eta = D(eta); mu = D(mu)
    M = [[2 * A[i][j] * (x * x - 1) + 4 * eta * B[i][j] * x + (mu * E[i] * (1 + x * x) if i == j else 0) for j in range(n)] for i in range(n)]
    return det(M)

def roots(ref, eta, mu, tol_digits):
    """2n roots k of det h(k). Cluster + near pi/2, cluster - near -pi/2."""
    n = len(ref['A']); mu = D(mu); th = ref['theta']; s = ref['s']
    out = []
    for sgn in (1, -1):
        ys = [(th[j] - sgn * s[j] * mu) / 2 for j in range(n)]  # k = sgn*pi/2 + 2y
        xmap = (lambda t: (1 + t) / (1 - t)) if sgn == 1 else (lambda t: (t - 1) / (1 + t))
        order = sorted(range(n), key=lambda j: ys[j])
        srt = [ys[j] for j in order]
        bounds = [srt[0] - Dm('0.05')] + [(srt[i] + srt[i + 1]) / 2 for i in range(n - 1)] + [srt[-1] + Dm('0.05')]
        f = lambda y: Pdet(ref, eta, mu, xmap(tan(y)))
        for i, j in enumerate(order):
            lo, hi = bounds[i], bounds[i + 1]
            flo, fhi = f(lo), f(hi)
            if flo == 0 or fhi == 0 or (flo > 0) == (fhi > 0):
                raise RuntimeError(f'no sign change n={n} j={j} sgn={sgn}')
            tol = Dm(10) ** (-tol_digits)
            # Illinois (bracketing regula falsi), bracket always kept
            side = 0; it = 0
            while hi - lo > tol and it < 20000:
                it += 1
                if it % 7 == 0:
                    c = (lo + hi) / 2
                else:
                    c = (lo * fhi - hi * flo) / (fhi - flo)
                    if not (lo < c < hi): c = (lo + hi) / 2
                fc = f(c)
                if fc == 0: lo = hi = c; break
                if (fc > 0) == (fhi > 0):
                    hi, fhi = c, fc
                    if side == 1: flo /= 2
                    side = 1
                else:
                    lo, flo = c, fc
                    if side == -1: fhi /= 2
                    side = -1
            y = (lo + hi) / 2
            out.append(((j, sgn), sgn * pi() / 2 + 2 * y, (hi - lo)))
    return out

def margins(rts):
    P = pi(); twopi = 2 * P
    def md(z):
        return z - twopi * (z / twopi).to_integral_value(rounding='ROUND_FLOOR')
    def circ(z):
        z = md(z)
        return min(z, twopi - z)
    ks = [k for _, k, _ in rts]; m = len(ks)
    d2 = min(circ(ks[i] - ks[j]) for i in range(m) for j in range(i))
    sums = [ks[i] + ks[j] for i, j in combinations_with_replacement(range(m), 2)]
        # distinct unordered pairs: min circular distance among all sums
    best = None
    zs = [md(z) for z in sums]
    zs_sorted = sorted(zs)
    for i in range(len(zs_sorted)):
        dlt = md(zs_sorted[(i + 1) % len(zs_sorted)] - zs_sorted[i])
        if best is None or dlt < best: best = dlt
    return d2, best

def check_lemma(n, ref, v, eta):
    """r_j / (kappa_j v^{4j}) and eigenvector tails vs (14)."""
    A = ref['A']; C0f = C0(n); v = D(v)
    detA = det(A)
    def dsub(M, lo, hi):
        if hi <= lo: return Dm(1)
        return det([[M[i][j] for j in range(lo, hi)] for i in range(lo, hi)])
    C0d = [[D(z) for z in row] for row in C0f]
    c = [dsub(C0d, 0, k) * dsub(A, k, n) / detA for k in range(n + 1)]
    kappa = [c[j + 1] / c[j] for j in range(n)]
    ratio = [ref['r'][j] / (kappa[j] * v ** (4 * j)) for j in range(n)]
    # limiting eigenvector: tail = A[j:,j:]^{-1} e0 normalized
    tail_err = []
    for j in range(n):
        T = [[A[p][q] for q in range(j, n)] for p in range(j, n)]
        m = n - j
        # solve T y = e0 by Cramer-free elimination
        Mx = [T[p][:] + [Dm(int(p == 0))] for p in range(m)]
        for cidx in range(m):
            pv = max(range(cidx, m), key=lambda p: abs(Mx[p][cidx])); Mx[cidx], Mx[pv] = Mx[pv], Mx[cidx]
            for p in range(m):
                if p != cidx:
                    f = Mx[p][cidx] / Mx[cidx][cidx]
                    Mx[p] = [Mx[p][t] - f * Mx[cidx][t] for t in range(m + 1)]
        y = [Mx[p][m] / Mx[p][p] for p in range(m)]
        y = [z / y[0] for z in y]
        x0 = [Dm(0)] * j + y
        tail_err.append(max(abs(ref['X'][j][p] - x0[p]) for p in range(n)))
    return kappa, ratio, tail_err

def fl(z):
    return float(z) if abs(z) > Dm('1e-300') else ('%.6e' % z)

def sci(z):
    return format(Dm(z), '.6e')

def run_point(n, w, v, u, eta, mu, prec, tol_digits):
    setprec(prec); t0 = time.time()
    ref = reference(n, w, v, u, eta)
    kappa, ratio, tail = check_lemma(n, ref, v, eta)
    rts = roots(ref, eta, mu, tol_digits)
    d2, d4 = margins(rts)
    th = ref['theta']; s = ref['s']
    bound17 = min(abs(th[-1]) / 4, D(mu) * s[-1] / 4)
    res = dict(n=n, w=str(w), v=str(v), u=str(u), eta=str(eta), mu=str(mu),
               theta=[sci(z) for z in th], slopes=[sci(z) for z in s],
               center_ratio=[sci(th[j + 1] / th[j]) for j in range(n - 1)],
               slope_ratio=[sci(s[j + 1] / s[j]) for j in range(n - 1)],
               r_over_kappa_v4j=[sci(z) for z in ratio], tail_err=[sci(z) for z in tail],
               delta2=sci(d2), delta4=sci(d4), bound17=sci(bound17), delta4_ge_bound17=bool(d4 >= bound17),
               lacunary_ok=bool(all(abs(th[j + 1]) <= abs(th[j]) / 32 for j in range(n - 1)) and all(s[j + 1] < s[j] / 32 for j in range(n - 1)) and abs(th[0]) < Dm('0.25')),
               seconds=round(time.time() - t0, 1))
    return res

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'fixed':
        out = []
        for n in [int(z) for z in sys.argv[2].split(',')]:
            prec = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
            mu = Fr(1, 10 ** (8 * n + 10))
            res = run_point(n, Fr(1, 4), Fr(1, 1024), Fr(1, 64), Fr(1, 100), mu, prec, tol_digits=min(prec - 60, 50 + 16 * n))
            print(json.dumps(res), flush=True); out.append(res)
        json.dump(out, open(f'fixed_{sys.argv[2].replace(",", "_")}.json', 'w'), indent=1)


def limiting_slopes(n, A, kappa0, u, eta):
    """s_j^0 from the limiting vectors (14); r_0(0)=kappa_0, r_j(0)=0 for j>=1."""
    u = D(u); eta = D(eta); out = []
    for j in range(n):
        m = n - j
        Mx = [[A[p][q] for q in range(j, n)] + [Dm(int(p == j))] for p in range(j, n)]
        for c in range(m):
            pv = max(range(c, m), key=lambda p: abs(Mx[p][c])); Mx[c], Mx[pv] = Mx[pv], Mx[c]
            for p in range(m):
                if p != c:
                    f = Mx[p][c] / Mx[c][c]; Mx[p] = [Mx[p][t] - f * Mx[c][t] for t in range(m + 1)]
        y = [Mx[p][m] / Mx[p][p] for p in range(m)]; y = [z / y[0] for z in y]
        x = [Dm(0)] * j + y
        num = sum(u ** (4 * i) * x[i] ** 2 for i in range(n))
        rj = kappa0 if j == 0 else Dm(0)
        den = 2 * (1 + (eta * rj) ** 2).sqrt() * sum(x[i] * A[i][k] * x[k] for i in range(n) for k in range(n))
        out.append(num / den)
    return out


def search(n, prec, mu_step=4, mu_max_exp=2000):
    setprec(prec); t0 = time.time(); w = Fr(1, 4); log = {}
    Gw = gram(n, D(w)); A = [[Gw[i][j] / D(w) ** (2 * i + 2 * j) for j in range(n)] for i in range(n)]
    C0d = [[D(z) for z in row] for row in C0(n)]
    detA = det(A)
    sub = lambda M, lo, hi: det([[M[i][j] for j in range(lo, hi)] for i in range(lo, hi)]) if hi > lo else Dm(1)
    c = [sub(C0d, 0, k) * sub(A, k, n) / detA for k in range(n + 1)]
    kappa = [c[j + 1] / c[j] for j in range(n)]
    m = 0
    while not (Fr(1, 2 ** m) < 1 and D(Fr(1, 2 ** m)) * kappa[0] < Dm(1) / 8): m += 1
    eta = Fr(1, 2 ** m); log['eta'] = str(eta); log['kappa'] = [sci(z) for z in kappa]
    m = 1
    while True:
        s0 = limiting_slopes(n, A, kappa[0], Fr(1, 2 ** m), eta)
        if all(s0[j + 1] < s0[j] / 64 for j in range(n - 1)): break
        m += 1
    u = Fr(1, 2 ** m); log['u'] = str(u)
    m = 1
    while True:
        v = Fr(1, 2 ** m)
        ref = reference(n, w, v, u, eta); th = ref['theta']; s = ref['s']
        ok = abs(th[0]) < Dm('0.25') and all(abs(th[j + 1]) <= abs(th[j]) / 32 for j in range(n - 1)) \
            and all(s[j + 1] < s[j] / 32 for j in range(n - 1)) and len(set(ref['r'])) == n
        if ok: break
        m += 1
    log['v'] = str(v); log['theta'] = [sci(z) for z in th]; log['slopes'] = [sci(z) for z in s]
    log['center_ratio'] = [sci(th[j + 1] / th[j]) for j in range(n - 1)]
    log['slope_ratio'] = [sci(s[j + 1] / s[j]) for j in range(n - 1)]
    m = mu_step; tried = []
    while m <= mu_max_exp:
        mu = Fr(1, 2 ** m)
        est = -int(min(abs(th[-1]), D(mu) * s[-1]).log10()) + 40
        try:
            rts = roots(ref, eta, mu, est)
            d2, d4 = margins(rts)
            passed = d4 > Dm(10) ** (-(est - 10))
        except RuntimeError as e:
            passed = False; d2 = d4 = None
        b17 = min(abs(th[-1]) / 4, D(mu) * s[-1] / 4)
        tried.append(dict(mu_exp2=m, passed=bool(passed), delta4=None if d4 is None else sci(d4), bound17=sci(b17),
                          ge17=None if d4 is None else bool(d4 >= b17)))
        if passed and d4 >= b17: break
        m += mu_step
    log['mu_trials'] = tried; log['mu'] = f'2^-{m}'
    log['delta2'] = sci(d2); log['delta4'] = sci(d4); log['seconds'] = round(time.time() - t0, 1)
    log['n'] = n
    return log


if __name__ == '__main__' and sys.argv[1] == 'search':
    n = int(sys.argv[2]); prec = int(sys.argv[3])
    res = search(n, prec)
    print(json.dumps(res), flush=True)
    json.dump(res, open(f'search_n{n}.json', 'w'), indent=1)
