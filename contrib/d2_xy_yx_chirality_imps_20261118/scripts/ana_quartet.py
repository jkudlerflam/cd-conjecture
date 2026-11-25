"""Momentum-zero exponents of the quartet, pair and density correlators from meas_quartet.py output (2026-11-24, Claude Opus 5.5).

The k = 0 part is taken by a boxcar average over W = 34 sites (the wedge's pitch, so the incommensurate harmonics at
multiples of 2 pi/34 average out), then averaged over the reference sites.  For the pair the ims frame shifts momentum
by pi, so ims k = 0 is the original-frame pi pair of REPORT 9.216.  For the quartet both frames agree.

Fits on r in [r1, r2]:  (a) one power law |G| = B r^-p;  (b) A r^-2 + B r^-p (a current-like r^-2 piece plus a scalar);
local slopes -dln|G|/dln r between r and 1.5 r.
usage: python ana_quartet.py measq/T1_L384.json ...
"""
import json, sys
import numpy as np
from scipy.optimize import curve_fit

W = 34


KAPPA = None   # incommensurate wavevector, fitted per file from the pair correlator


def harmonics(kappa):
    qs = []
    for m in (0, 1):
        for n in range(-2, 3):
            q = (m * np.pi + n * kappa) % (2 * np.pi)
            if all(abs((q - p + np.pi) % (2 * np.pi) - np.pi) > 1e-6 for p in qs):
                qs.append(q)
    return np.array(qs)


def k0(r, v, q0=0.0, Wd=48, step=4):
    """local least-squares demodulation: in a Hann window of width Wd around each centre, fit v(r) = sum_q c_q e^{iqr}
    over the harmonic set {m pi + n kappa}, and return the coefficient at q0 (default the k = 0 part)."""
    qs = harmonics(KAPPA)
    k = int(np.argmin(np.abs((qs - q0 + np.pi) % (2 * np.pi) - np.pi)))
    rc, c = [], []
    for a in range(0, len(r) - Wd, step):
        rr, vv = r[a:a + Wd].astype(float), v[a:a + Wd]
        w = np.sqrt(np.hanning(Wd + 2)[1:-1])
        X = np.exp(1j * np.outer(rr - rr.mean(), qs)) * w[:, None]
        coef = np.linalg.lstsq(X, vv * w, rcond=None)[0]
        rc.append(rr.mean()); c.append(coef[k] * np.exp(1j * qs[k] * rr.mean()))
    return np.array(rc), np.array(c)


def fit_kappa(blk, i0):
    b = blk[str(i0)]; v = np.array(b["md"]); v = v[:, 0] + 1j * v[:, 1]
    r = np.arange(len(v)) + b["j0"] - i0; s = r >= 16
    q = np.linspace(0.15, 0.27, 2401) * np.pi
    F = np.abs(np.exp(-1j * np.outer(q, r[s])) @ (v[s] * np.hanning(s.sum())))
    return q[np.argmax(F)]


def curve(blk, refs, lab, q0=0.0):
    """ref-averaged k = q0 curve for correlator label lab, on a common r grid."""
    out = []
    for i0 in refs:
        if str(i0) not in blk:
            continue
        b = blk[str(i0)]
        v = np.array(b[lab]); v = v[:, 0] + 1j * v[:, 1]
        r = np.arange(len(v)) + b["j0"] - i0
        out.append(k0(r, v, q0))
    n = min(len(c) for _, c in out)
    return out[0][0][:n], np.mean([c[:n] for _, c in out], axis=0)


def quad(blk, refs):
    """k = 0 correlation matrix of the Hermitian quadratures O1 = (Q + Q^dag)/2, O2 = (Q - Q^dag)/2i; returns r and
    its two eigenvalues (sorted by magnitude at large r).  A U(1)-charged vertex fills both quadratures equally; a
    neutral admixture with one fixed phase is rank one, so the smaller eigenvalue isolates the vertex."""
    c = {lab: curve(blk, refs, lab) for lab in ("md", "dm", "mm", "dd")}
    r = c["md"][0]; md, dm, mm, dd = (c[k][1] for k in ("md", "dm", "mm", "dd"))
    C11 = ((mm + md + dm + dd) / 4).real
    C22 = (-(mm - md - dm + dd) / 4).real
    C12 = ((mm - md + dm - dd) / 4j).real
    C21 = ((mm + md - dm - dd) / 4j).real
    off = (C12 + C21) / 2
    tr, det = C11 + C22, C11 * C22 - off ** 2
    disc = np.sqrt(np.maximum(tr ** 2 / 4 - det, 0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    big = np.where(np.abs(l1) >= np.abs(l2), l1, l2); small = np.where(np.abs(l1) >= np.abs(l2), l2, l1)
    return r, big, small, np.abs(C12 - C21) / 2


def fits(r, g, L):
    res = {}
    for r1, r2 in ((24, L / 8), (32, L / 6), (40, L / 4), (64, L / 4)):
        s = (r >= r1) & (r <= r2)
        if s.sum() < 8:
            continue
        x, y = np.log(r[s]), np.log(np.abs(g[s]))
        p1 = -np.polyfit(x, y, 1)[0]
        try:
            f = lambda rr, A, B, p: A * rr ** -2.0 + B * rr ** -p
            po, pc = curve_fit(f, r[s], g[s], p0=(g[s][-1] * r[s][-1] ** 2, 0.0, 3.7), maxfev=20000)
            p2, e2, A2, B2 = po[2], np.sqrt(abs(pc[2, 2])), po[0], po[1]
        except Exception:
            p2 = e2 = A2 = B2 = np.nan
        res[(r1, int(r2))] = (p1, p2, e2, A2, B2)
    return res


def local(r, g):
    out = []
    for rc in (24, 32, 48, 64, 96, 128):
        a, b = np.argmin(abs(r - rc)), np.argmin(abs(r - 1.5 * rc))
        if b > a and r[b] > 1.4 * rc:
            out.append((rc, -np.log(abs(g[b] / g[a])) / np.log(r[b] / r[a])))
    return out


for fn in sys.argv[1:]:
    m = json.load(open(fn)); L = m["L"]; refs = m["refs"]
    KAPPA = fit_kappa(m["P"]["2"], refs[0])
    print("== %s  L=%d chi=%d refs=%s  kappa/pi=%.4f" % (fn, L, m["chi"], refs, KAPPA / np.pi))
    # density
    zz = []
    for i0 in refs:
        v = np.array(m["ZZ"][str(i0)]); r = np.arange(1, len(v) + 1); zz.append(k0(r, v))
    n = min(len(c) for _, c in zz); rz = zz[0][0][:n]; gz = np.mean([c[:n] for _, c in zz], axis=0)
    print("  ZZ k0: G(32,64,96)=%s  local %s" % (" ".join("%.2e" % gz[np.argmin(abs(rz - x))] for x in (32, 64, 96)),
                                                " ".join("%d:%.2f" % t for t in local(rz, gz))))
    for kind in ("P", "Q"):
        for key, blk in m[kind].items():
            for lab in ("md", "mm", "mz"):
                r, g = curve(blk, refs, lab)
                # use the component with the larger weight (real or imaginary)
                gr = g.real if np.sum(g.real ** 2) >= np.sum(g.imag ** 2) else g.imag
                ft = fits(r, gr, L)
                sgn = "+" if np.median(gr[(r > 32) & (r < L / 4)]) > 0 else "-"
                print("  %s %-4s %-2s sign %s |G|(32,64,96)=%s  local %s" % (
                    kind, key, lab, sgn, " ".join("%.2e" % abs(gr[np.argmin(abs(r - x))]) for x in (32, 64, 96)),
                    " ".join("%d:%.2f" % t for t in local(r, gr))))
                for w, (p1, p2, e2, A2, B2) in ft.items():
                    print("        window %-9s one-power p=%.3f   A/r^2+B/r^p: p=%.3f+-%.3f A=%.2e B=%.2e" % (w, p1, p2, e2, A2, B2))
            r, big, small, asym = quad(blk, refs)
            for name, g in (("eig+", big), ("eig-", small)):
                sel = (r > 32) & (r < L / 4)
                print("  %s %-4s %-4s sign %s |G|(32,64,96)=%s  local %s   (antisym part/|eig-| at r~64: %.2f)" % (
                    kind, key, name, "+" if np.median(g[sel]) > 0 else "-",
                    " ".join("%.2e" % abs(g[np.argmin(abs(r - x))]) for x in (32, 64, 96)),
                    " ".join("%d:%.2f" % t for t in local(r, g)), asym[np.argmin(abs(r - 64))] / abs(small[np.argmin(abs(r - 64))])))
                for w, (p1, p2, e2, A2, B2) in fits(r, g, L).items():
                    print("        window %-9s one-power p=%.3f   A/r^2+B/r^p: p=%.3f+-%.3f A=%.2e B=%.2e" % (w, p1, p2, e2, A2, B2))
