"""High-precision complex matrices on Python's decimal module (stdlib only).

A complex matrix is a CM holding two numpy object arrays of Decimal (real and
imaginary parts).  Products use numpy's object dot, so each entry of a product
is an exact sum of exactly rounded Decimal products at the context precision.
"""
from decimal import Decimal, getcontext, localcontext
import numpy as np

ZERO = Decimal(0)
ONE = Decimal(1)


def setprec(p):
    getcontext().prec = p


def dzeros(shape):
    return np.full(shape, ZERO, dtype=object)


def todec(x):
    """float / complex / ndarray -> Decimal object array (exact conversion)."""
    a = np.asarray(x)
    out = np.empty(a.shape, dtype=object)
    for idx, v in np.ndenumerate(a):
        out[idx] = Decimal(float(v))
    return out


def tofloat(a):
    return np.array(a, dtype=float)


class CM:
    __slots__ = ('r', 'i')

    def __init__(self, r, i=None):
        self.r = r
        self.i = dzeros(r.shape) if i is None else i

    @staticmethod
    def zeros(shape):
        return CM(dzeros(shape), dzeros(shape))

    @staticmethod
    def eye(m):
        r = dzeros((m, m))
        for k in range(m):
            r[k, k] = ONE
        return CM(r)

    @staticmethod
    def fromcomplex(z):
        z = np.asarray(z, dtype=complex)
        return CM(todec(z.real), todec(z.imag))

    def tocomplex(self):
        return tofloat(self.r) + 1j * tofloat(self.i)

    @property
    def shape(self):
        return self.r.shape

    def copy(self):
        return CM(self.r.copy(), self.i.copy())

    def __add__(self, o):
        return CM(self.r + o.r, self.i + o.i)

    def __sub__(self, o):
        return CM(self.r - o.r, self.i - o.i)

    def __neg__(self):
        return CM(-self.r, -self.i)

    def __matmul__(self, o):
        # Gauss three-product form
        k1 = np.dot(self.r + self.i, o.r)
        k2 = np.dot(self.r, o.i - o.r)
        k3 = np.dot(self.i, o.r + o.i)
        return CM(k1 - k3, k1 + k2)

    def H(self):
        return CM(self.r.T.copy(), -self.i.T)

    def T(self):
        return CM(self.r.T.copy(), self.i.T.copy())

    def conj(self):
        return CM(self.r.copy(), -self.i)

    def emul(self, o):
        """elementwise product with another CM"""
        return CM(self.r * o.r - self.i * o.i, self.r * o.i + self.i * o.r)

    def rmul(self, a):
        """elementwise / scalar product with a real Decimal (array)"""
        return CM(self.r * a, self.i * a)

    def cscale(self, cr, ci):
        """multiply by the complex scalar cr + i ci"""
        return CM(self.r * cr - self.i * ci, self.r * ci + self.i * cr)

    def __getitem__(self, k):
        return CM(self.r[k], self.i[k])

    def __setitem__(self, k, v):
        self.r[k] = v.r
        self.i[k] = v.i

    def absmax(self):
        a = np.abs(self.tocomplex())
        return float(a.max()) if a.size else 0.0


# ------------------------------------------------------------- scalars
def dpi():
    with localcontext() as ctx:
        ctx.prec += 5
        three = Decimal(3)
        lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            n, na = n + na, na + 8
            d, da = d + da, da + 32
            t = (t * n) / d
            s += t
    return +s


def dcos(x):
    with localcontext() as ctx:
        ctx.prec += 5
        x = Decimal(x)
        twopi = 2 * dpi()
        x = x - twopi * (x / twopi).to_integral_value()
        i, lasts, s, fact, num, sign = 0, 0, Decimal(1), 1, Decimal(1), 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i - 1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
    return +s


def dsin(x):
    with localcontext() as ctx:
        ctx.prec += 5
        x = Decimal(x)
        twopi = 2 * dpi()
        x = x - twopi * (x / twopi).to_integral_value()
        i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i - 1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
    return +s


def dexpi(phi):
    """(cos phi, sin phi) as Decimals, exact for multiples of pi/2 given as
    strings 'pi/2', 'pi', '0', '-pi/2'."""
    exact = {'0': (ONE, ZERO), 'pi/2': (ZERO, ONE), 'pi': (-ONE, ZERO), '-pi/2': (ZERO, -ONE)}
    if isinstance(phi, str) and phi in exact:
        return exact[phi]
    return dcos(phi), dsin(phi)


# ------------------------------------------------------------- eigen
def eigh_hp(A, maxit=20, verbose=False):
    """Hermitian eigendecomposition at the current Decimal precision.
    Float eigh start, then Ogita-Aishima Newton refinement (quadratic
    convergence for simple eigenvalues).  Returns (lam: object array, V: CM)
    with columns of V the eigenvectors, and the final residual norms."""
    m = A.shape[0]
    w, Vf = np.linalg.eigh(A.tocomplex())
    V = CM.fromcomplex(Vf)
    I = CM.eye(m)
    prec = getcontext().prec
    target = Decimal(10) ** (-(prec - 8))
    for it in range(maxit):
        VH = V.H()
        R = I - VH @ V
        S = VH @ (A @ V)
        lam = np.array([S.r[k, k] / (ONE - R.r[k, k]) for k in range(m)], dtype=object)
        den = lam[None, :] - lam[:, None]
        for k in range(m):
            den[k, k] = ONE
        num = CM(S.r + R.r * lam[None, :], S.i + R.i * lam[None, :])
        E = CM(num.r / den, num.i / den)
        for k in range(m):
            E.r[k, k] = R.r[k, k] / 2
            E.i[k, k] = ZERO
        V = V + V @ E
        # size of the correction in log10 (float of Decimal handles tiny values)
        emax = max(max(abs(x) for x in E.r.ravel()), max(abs(x) for x in E.i.ravel()))
        if verbose:
            print('  eig it', it, float(emax))
        if emax < target:
            break
    VH = V.H()
    S = VH @ (A @ V)
    lam = np.array([S.r[k, k] for k in range(m)], dtype=object)
    orth = max(abs(x) for x in (VH @ V - I).r.ravel())
    off = S.copy()
    for k in range(m):
        off.r[k, k] = ZERO
    offmax = max(max(abs(x) for x in off.r.ravel()), max(abs(x) for x in off.i.ravel()))
    return lam, V, dict(iterations=it + 1, orth=orth, offdiag=offmax)


# ------------------------------------------------------------- dense algebra
def cholesky_hp(A):
    """A = L L^H for Hermitian positive definite A (CM).  Returns L and the
    real pivots (squared diagonal of L).  Raises if a pivot is not positive."""
    n = A.shape[0]
    Lr, Li = dzeros((n, n)), dzeros((n, n))
    piv = []
    for j in range(n):
        s = A.r[j, j] - sum((Lr[j, k] ** 2 + Li[j, k] ** 2 for k in range(j)), ZERO)
        piv.append(s)
        if s <= 0:
            raise ValueError('not positive definite at pivot %d' % j)
        d = s.sqrt()
        Lr[j, j] = d
        for i in range(j + 1, n):
            # A_ij - sum_k L_ik conj(L_jk)
            sr = A.r[i, j] - sum((Lr[i, k] * Lr[j, k] + Li[i, k] * Li[j, k] for k in range(j)), ZERO)
            si = A.i[i, j] - sum((Li[i, k] * Lr[j, k] - Lr[i, k] * Li[j, k] for k in range(j)), ZERO)
            Lr[i, j] = sr / d
            Li[i, j] = si / d
    return CM(Lr, Li), piv


def lower_inverse_hp(L):
    """inverse of a lower-triangular CM with real positive diagonal"""
    n = L.shape[0]
    W = CM.zeros((n, n))
    for j in range(n):
        W.r[j, j] = ONE / L.r[j, j]
        for i in range(j + 1, n):
            sr, si = ZERO, ZERO
            for k in range(j, i):
                ar, ai = L.r[i, k], L.i[i, k]
                br, bi = W.r[k, j], W.i[k, j]
                sr += ar * br - ai * bi
                si += ar * bi + ai * br
            W.r[i, j] = -sr / L.r[i, i]
            W.i[i, j] = -si / L.r[i, i]
    return W


def lu_det_logderiv(Hr, Hi, Dr, Di):
    """For square complex matrices H, Dh given as Decimal object arrays, return
    (det H, tr(H^-1 Dh)) by Gaussian elimination with partial pivoting."""
    n = Hr.shape[0]
    # augmented complex matrix [H | Dh]
    Ar = np.concatenate([Hr, Dr], axis=1).copy()
    Ai = np.concatenate([Hi, Di], axis=1).copy()
    detr, deti = ONE, ZERO
    for k in range(n):
        mags = [Ar[i, k] * Ar[i, k] + Ai[i, k] * Ai[i, k] for i in range(k, n)]
        p = k + int(np.argmax(np.array([float(x) if x != 0 else 0.0 for x in mags])))
        if p != k:
            Ar[[k, p]] = Ar[[p, k]]
            Ai[[k, p]] = Ai[[p, k]]
            detr, deti = -detr, -deti
        pr, pi_ = Ar[k, k], Ai[k, k]
        detr, deti = detr * pr - deti * pi_, detr * pi_ + deti * pr
        nrm = pr * pr + pi_ * pi_
        ir, ii = pr / nrm, -pi_ / nrm
        for i in range(k + 1, n):
            fr = Ar[i, k] * ir - Ai[i, k] * ii
            fi = Ar[i, k] * ii + Ai[i, k] * ir
            if fr == 0 and fi == 0:
                continue
            rowr, rowi = Ar[k, k:], Ai[k, k:]
            Ar[i, k:] = Ar[i, k:] - (rowr * fr - rowi * fi)
            Ai[i, k:] = Ai[i, k:] - (rowr * fi + rowi * fr)
    # back substitution for X = H^-1 Dh, only the trace is needed
    m = Dr.shape[1]
    Xr, Xi = dzeros((n, m)), dzeros((n, m))
    for i in range(n - 1, -1, -1):
        sr = Ar[i, n:].copy()
        si = Ai[i, n:].copy()
        for k in range(i + 1, n):
            sr = sr - (Ar[i, k] * Xr[k] - Ai[i, k] * Xi[k])
            si = si - (Ar[i, k] * Xi[k] + Ai[i, k] * Xr[k])
        pr, pi_ = Ar[i, i], Ai[i, i]
        nrm = pr * pr + pi_ * pi_
        Xr[i] = (sr * pr + si * pi_) / nrm
        Xi[i] = (si * pr - sr * pi_) / nrm
    trr = sum((Xr[i, i] for i in range(n)), ZERO)
    tri = sum((Xi[i, i] for i in range(n)), ZERO)
    return (detr, deti), (trr, tri)


def log10abs(x):
    x = abs(x)
    if x == 0:
        return float('-inf')
    return float(x.log10())
