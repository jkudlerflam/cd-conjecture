"""Bands, Fermi points, centres, velocities and resonance margins of
h(k) = mu + eps + T e^{ik} + T^dag e^{-ik}, from high-precision T and mu.

Whitening: with R = e^{-i Phi} T, R_s = Herm(R) = L L^dag > 0 and W = L^-1,
    W h W^dag = 2 cos q I - 2 sin q A + B,   q = k + Phi,
    A = W R_a W^dag,  B = W (mu + eps) W^dag,  R_a = (R - R^dag)/(2i).
Congruence keeps the zeros and the inertia of h, and the whitened problem is
well conditioned, so roots are found in double precision and then refined by
Newton's method on det in high precision.  All margins are computed from the
refined roots in Decimal arithmetic.
"""
import itertools, math
from decimal import Decimal, localcontext, getcontext, ROUND_FLOOR
import numpy as np
from scipy.linalg import eig as geig
from hpmat import (CM, ZERO, ONE, dzeros, dpi, dcos, dsin, cholesky_hp,
                   lower_inverse_hp, lu_det_logderiv, todec, log10abs)


def herm_parts(R):
    RH = R.H()
    Rs = CM((R.r + RH.r) / 2, (R.i + RH.i) / 2)
    X = R - RH                       # anti-Hermitian; R_a = X/(2i) = -i X / 2
    Ra = CM(X.i / 2, -X.r / 2)
    return Rs, Ra


def extreme_eigs(Rs, W, iters=60):
    """lambda_min(R_s) by power iteration on R_s^-1 = W^dag W, lambda_max by power
    iteration on R_s, both in Decimal.  Graded spectra are fine here."""
    n = Rs.shape[0]
    WH = W.H()
    x = CM(np.array([Decimal(1 + k % 3) for k in range(n)], dtype=object), dzeros(n))

    def mv(Mx, v):
        return CM(np.dot(Mx.r, v.r) - np.dot(Mx.i, v.i), np.dot(Mx.r, v.i) + np.dot(Mx.i, v.r))

    def nrm(v):
        return (sum((a * a for a in v.r), ZERO) + sum((a * a for a in v.i), ZERO)).sqrt()
    res = []
    for op in ('min', 'max'):
        v = x.copy()
        lam = None
        for _ in range(iters):
            y = mv(WH, mv(W, v)) if op == 'min' else mv(Rs, v)
            s = nrm(y)
            v = CM(y.r / s, y.i / s)
        y = mv(Rs, v)
        rq = sum((a * b for a, b in zip(v.r, y.r)), ZERO) + sum((a * b for a, b in zip(v.i, y.i)), ZERO)
        res.append(rq)
    return res[0], res[1]


def whiten(T, mu, eps, phi='pi'):
    n = T.shape[0]
    if phi == 'pi':
        R = -T
        Phi = math.pi
    else:
        tr = complex(sum(T.tocomplex()[k, k] for k in range(n)))
        Phi = math.atan2(tr.imag, tr.real)
        c, s = Decimal(math.cos(Phi)), Decimal(math.sin(Phi))
        R = T.cscale(c, -s)
    Rs, Ra = herm_parts(R)
    L, piv = cholesky_hp(Rs)
    W = lower_inverse_hp(L)
    WH = W.H()
    A = W @ Ra @ WH
    Mt = W @ mu @ WH
    F = W @ WH
    B = Mt + F.rmul(eps)
    lmin, lmax = extreme_eigs(Rs, W)
    return dict(Phi=Phi, R=R, Rs=Rs, Ra=Ra, W=W, A=A, B=B, Mt=Mt, F=F, piv=piv,
                lmin=lmin, lmax=lmax)


def float_roots(Af, Bf):
    n = Af.shape[0]
    I = np.eye(n)
    P2, P1, P0 = I + 1j * Af, Bf, I - 1j * Af
    Z = np.zeros((n, n))
    Amat = np.block([[Z, I], [-P0, -P1]])
    Bmat = np.block([[I, Z], [Z, P2]])
    z = geig(Amat, Bmat, right=False)
    return z


def htilde(A, B, q, deriv=False):
    """whitened symbol and its q-derivative, Decimal (q Decimal)"""
    c, s = dcos(q), dsin(q)
    n = A.shape[0]
    Hr = B.r - 2 * s * A.r
    Hi = B.i - 2 * s * A.i
    for k in range(n):
        Hr[k, k] += 2 * c
    if not deriv:
        return Hr, Hi
    Dr = -2 * c * A.r
    Di = -2 * c * A.i
    for k in range(n):
        Dr[k, k] += -2 * s
    return Hr, Hi, Dr, Di


def refine_root(A, B, q0, digits, maxit=12):
    q = Decimal(q0)
    tol = Decimal(10) ** (-digits)
    for it in range(maxit):
        Hr, Hi, Dr, Di = htilde(A, B, q, deriv=True)
        (dr, di), (tr, ti) = lu_det_logderiv(Hr, Hi, Dr, Di)
        step = tr / (tr * tr + ti * ti)          # Re(1/tr), tr should be real
        # (the returned |ti/tr| is taken at the final, over-converged step, where
        #  h is singular to working precision, so it is not a meaningful diagnostic)
        q = q - step
        if abs(step) < tol:
            break
    return q, it + 1, abs(ti / tr) if tr != 0 else None


def dmod(x, m):
    """floor modulo for Decimals (Decimal % keeps the sign of x)"""
    return x - m * (x / m).to_integral_value(rounding=ROUND_FLOOR)


def circ(x, twopi):
    """distance of x to 2 pi Z"""
    y = dmod(x, twopi)
    return min(y, twopi - y)


def analyse(T, mu, eps, phi='pi', rdigits=60, verbose=False):
    n = T.shape[0]
    wh = whiten(T, mu, eps, phi)
    A, B = wh['A'], wh['B']
    Af, Bf = A.tocomplex(), B.tocomplex()
    Af = (Af + Af.conj().T) / 2
    Bf = (Bf + Bf.conj().T) / 2
    z = float_roots(Af, Bf)
    unit_dev = float(np.max(np.abs(np.abs(z) - 1)))
    q0 = np.sort(np.mod(np.angle(z), 2 * np.pi))
    twopi = 2 * dpi()
    # high-precision refinement on rounded whitened matrices
    with localcontext() as ctx:
        ctx.prec = rdigits + 25
        Ar = CM(np.array([[+x for x in r] for r in A.r], dtype=object), np.array([[+x for x in r] for r in A.i], dtype=object))
        Br = CM(np.array([[+x for x in r] for r in B.r], dtype=object), np.array([[+x for x in r] for r in B.i], dtype=object))
        roots, its, imres = [], [], []
        for qq in q0:
            r, it, im = refine_root(Ar, Br, qq, rdigits)
            roots.append(dmod(r, twopi))
            its.append(it)
            imres.append(float(im) if im is not None else None)
    order = sorted(range(len(roots)), key=lambda k: roots[k])
    roots = [roots[k] for k in order]
    m2 = len(roots)
    newton_shift = max(abs(float(roots[k]) - q0[k]) for k in range(m2))
    # inertia on each arc, float (whitened matrices are O(1))
    nus, gaps = [], []
    for k in range(m2):
        a, b = float(roots[k]), float(roots[(k + 1) % m2])
        if b <= a:
            b += 2 * np.pi
        qm = (a + b) / 2
        h = 2 * np.cos(qm) * np.eye(n) - 2 * np.sin(qm) * Af + Bf
        ev = np.linalg.eigvalsh(h)
        nus.append(int(np.sum(ev < 0)))
        gaps.append(float(np.min(np.abs(ev))))
    # nus[k] is the inertia on the arc after roots[k]
    steps = [(nus[k] - nus[k - 1]) for k in range(m2)]   # step at roots[k]
    ok_steps = all(abs(s) == 1 for s in steps)
    k0 = nus.index(0)                 # arc after roots[k0] has nu = 0
    seq = [(k0 + 1 + j) % m2 for j in range(m2)]       # roots in circular order after the empty arc
    up = seq[:n]
    dn = seq[n:]
    unimodal = all(steps[k] == 1 for k in up) and all(steps[k] == -1 for k in dn)
    # unwrap positions along the circle starting at the empty arc
    base = roots[seq[0]]
    pos = {}
    for k in seq:
        x = roots[k]
        while x < base:
            x += twopi
        pos[k] = x
    bands = []
    for mb in range(1, n + 1):
        kl, kr = up[mb - 1], dn[n - mb]
        pl, pr = pos[kl], pos[kr]
        bands.append(dict(kL=kl, kR=kr, pL=pl, pR=pr, centre=(pl + pr) / 2, kF=(pr - pl) / 2))
    # margins (all in Decimal)
    P = roots
    d2, w2 = None, None
    for i in range(m2):
        dd = circ(P[(i + 1) % m2] - P[i], twopi)
        if d2 is None or dd < d2:
            d2, w2 = dd, (i, (i + 1) % m2)
    pairs = [(i, j) for i in range(m2) for j in range(i, m2)]
    sums = sorted((dmod(P[i] + P[j], twopi), (i, j)) for i, j in pairs)
    d4, w4 = None, None
    gl = []
    for k in range(len(sums)):
        a, b = sums[k], sums[(k + 1) % len(sums)]
        dd = circ(b[0] - a[0], twopi)
        gl.append((dd, a[1], b[1]))
        if d4 is None or dd < d4:
            d4, w4 = dd, (a[1], b[1])
    small_gaps = sorted(gl, key=lambda x: x[0])[:3 * n + 5]
    th = [b['centre'] for b in bands]
    kf = [b['kF'] for b in bands]
    dPT, wPT = None, None
    dF, dU2 = None, None
    for a, b in itertools.combinations(range(n), 2):
        dd = circ(2 * th[a] - 2 * th[b], twopi)
        if dPT is None or dd < dPT:
            dPT, wPT = dd, (a + 1, b + 1)
        e1 = circ(2 * kf[a] - 2 * kf[b], twopi)
        e2 = circ(2 * kf[a] + 2 * kf[b], twopi)
        dF = e1 if dF is None else min(dF, e1)
        dU2 = e2 if dU2 is None else min(dU2, e2)
    dU1 = min(circ(4 * x, twopi) for x in kf)
    # velocities: whitened null vector, original units via W W^dag
    Fd = wh['F']
    vel, velw = [], []
    for k in range(m2):
        qf = float(P[k])
        h = 2 * np.cos(qf) * np.eye(n) - 2 * np.sin(qf) * Af + Bf
        hp = -2 * np.sin(qf) * np.eye(n) - 2 * np.cos(qf) * Af
        ev, U = np.linalg.eigh(h)
        u = U[:, np.argmin(np.abs(ev))]
        vw = float(np.real(np.vdot(u, hp @ u)))
        ud = CM.fromcomplex(u[:, None])
        den = (ud.H() @ Fd @ ud).r[0, 0]
        velw.append(vw)
        vel.append(Decimal(abs(vw)) / den)
    Phi = Decimal(wh['Phi'])
    kroots = [dmod(r - Phi, twopi) for r in P]
    kroots = [x - twopi if x > twopi / 2 else x for x in kroots]
    centres_k = [dmod(b['centre'] - Phi + twopi / 2, twopi) - twopi / 2 for b in bands]
    out = dict(n=n, Phi=wh['Phi'], unit_circle_dev=unit_dev, newton_shift=newton_shift,
               newton_iters=max(its), newton_imag_ratio=max(x for x in imres if x is not None),
               inertia_steps_ok=ok_steps, unimodal=unimodal, nus=nus,
               min_abs_eig_on_arcs=min(gaps),
               chol_min_pivot=min(wh['piv']), lmin_Rs=wh['lmin'], lmax_Rs=wh['lmax'],
               normA=float(np.linalg.norm(Af, 2)), normB=float(np.linalg.norm(Bf, 2)),
               loop_whitened=float(np.linalg.norm(wh['Mt'].tocomplex(), 2)),
               field_whitened=float(np.linalg.norm(wh['F'].tocomplex(), 2)) * float(eps),
               roots_q=P, roots_k=kroots, centres_k=centres_k, kF=kf,
               bands=bands, delta2=d2, w2=w2, delta4=d4, w4=w4, small_gaps=small_gaps, dPT=dPT, wPT=wPT,
               dU1=dU1, dU2=dU2, dF=dF, vel=vel, vel_whitened=velw,
               vmin=min(vel), vmax=max(vel), vw_min=min(abs(x) for x in velw))
    return out, wh
