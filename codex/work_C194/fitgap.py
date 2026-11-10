"""C194 registered extrapolation of the pair-breaking gap (RESULT.md, operational definition 2).
Delta(L) = D + a/L + b/L^2 by least squares.  err = sqrt(standard error^2 + spread^2), where spread is the largest
|D - D_sub| over the two sub-fits (drop smallest L, keep b; all L, b = 0).  With three sizes the full fit is exact,
so its standard error is replaced by the spread alone."""
import numpy as np

def _lsq(Ls, d, quad):
    x = 1 / np.asarray(Ls, float); y = np.asarray(d, float)
    A = np.stack([np.ones_like(x), x] + ([x * x] if quad else []), 1)
    c, res, rk, sv = np.linalg.lstsq(A, y, rcond=None)
    dof = len(y) - A.shape[1]
    if dof > 0:
        r = y - A @ c; s2 = r @ r / dof; cov = s2 * np.linalg.inv(A.T @ A); se = np.sqrt(cov[0, 0])
    else:
        se = 0.0
    return c, se

def fit_gap(Ls, d):
    Ls = list(Ls); d = list(d)
    o = np.argsort(Ls); Ls = [Ls[i] for i in o]; d = [d[i] for i in o]
    c, se = _lsq(Ls, d, len(Ls) >= 3)
    subs = []
    if len(Ls) >= 4:
        subs.append(_lsq(Ls[1:], d[1:], True)[0][0])
    if len(Ls) >= 3:
        subs.append(_lsq(Ls, d, False)[0][0])
    spread = max([abs(c[0] - s) for s in subs], default=0.0)
    return dict(D=float(c[0]), a=float(c[1]), b=float(c[2]) if len(c) > 2 else 0.0, se=float(se), spread=float(spread),
                err=float(np.hypot(se, spread)) if (se or spread) else float('nan'))
