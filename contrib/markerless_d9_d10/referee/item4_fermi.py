"""Item 4: Fermi points of h(k) = mu^L + mu^R + f + T e^{ik} + T^dag e^{-ik} by a method
different from the author's companion-matrix eigenvalue solve.

  1. d(k) = det h(k) is real (h(k) Hermitian).  Count its sign changes on a uniform k grid,
     then bisect each bracket.  z^n d is a polynomial of degree 2n in z = e^{ik}, so 2n sign
     changes prove exactly 2n simple zeros on the unit circle (given correct signs).
  2. Chirality from the inertia of h(k): N_-(k) = number of negative eigenvalues, evaluated
     at the midpoints between consecutive zeros by an LDL^dag count (not by eigenvectors).
     A zero where N_- drops by one is a right-mover (eigenvalue goes - -> + with k).
  3. Root separation and pair-sum separation recomputed from the bisected zeros.
  4. Independent second count: winding number of z -> det(T z^2 + M z + T^dag) on circles
     |z| = 1 -+ delta (argument principle) -- how many of the 2n zeros lie strictly inside.

usage: item4_fermi.py blocks.json f_gamma [grid]
   f_gamma = number            : field = f_gamma * (smallest eigenvalue of the Hermitian part of e^{-i theta} T)
   f_gamma = beta2:2^-83       : field = beta^2 * 2^-83  (the certificate's field in units lambda^2 alpha^4 g^2)
"""
import sys, os, json, time
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))


def load(fn):
    d = json.load(open(fn))
    mp.mp.dps = d["dps"]
    def mat(L):
        n = len(L)
        M = mp.zeros(n, n)
        for i in range(n):
            for j in range(n):
                M[i, j] = mp.mpc(mp.mpf(L[i][j][0]), mp.mpf(L[i][j][1]))
        return M
    return d, mat(d["T"]), mat(d["muL"]), mat(d["muR"])


def inertia_neg(Hm):
    """number of negative eigenvalues of a Hermitian mp matrix via Gaussian elimination with
    symmetric diagonal pivoting (Sylvester's law of inertia)."""
    n = Hm.rows
    A = Hm.copy()
    neg = 0
    idx = list(range(n))
    for s in range(n):
        # pivot: largest remaining diagonal entry in magnitude
        p = max(range(s, n), key=lambda i: abs(A[i, i].real))
        if p != s:
            for c in range(n):
                A[s, c], A[p, c] = A[p, c], A[s, c]
            for r in range(n):
                A[r, s], A[r, p] = A[r, p], A[r, s]
        piv = A[s, s].real
        if piv == 0:
            raise ZeroDivisionError("zero pivot")
        if piv < 0:
            neg += 1
        for i in range(s + 1, n):
            fac = A[i, s] / piv
            for c in range(s + 1, n):
                A[i, c] -= fac * A[s, c]
            A[i, s] = 0
    return neg


def main(fn, fgam, grid):
    d, T, muL, muR = load(fn)
    n = d["n"]
    mu = muL + muR
    herm = max(abs(mu[i, j] - mp.conj(mu[j, i])) for i in range(n) for j in range(n))
    scT = max(abs(T[i, j]) for i in range(n) for j in range(n))
    th = mp.arg(sum(T[i, i] for i in range(n)))
    Th = T * mp.expj(-th)
    S = (Th + Th.H) / 2
    es = mp.eigh(S, eigvals_only=True)
    gam_min = min(es)
    if str(fgam).startswith("beta2:2^"):
        f = mp.mpf(d["beta"]) ** 2 * mp.mpf(2) ** int(str(fgam)[8:])
    else:
        f = mp.mpf(fgam) * gam_min
    M = (mu + mu.H) / 2
    for i in range(n):
        M[i, i] += f
    TH = T.H
    def hk(k):
        z = mp.expj(k)
        return M + T * z + TH / z
    def dk(k):
        return mp.det(hk(k)).real
    t0 = time.time()
    twopi = 2 * mp.pi
    ks = [-mp.pi + twopi * i / grid for i in range(grid + 1)]
    vals = [dk(k) for k in ks]
    imag_max = max(abs(mp.det(hk(k)).imag) for k in ks[::max(1, grid // 50)])
    brackets = [(ks[i], ks[i + 1]) for i in range(grid) if (vals[i] > 0) != (vals[i + 1] > 0)]
    roots = []
    tol = mp.mpf(10) ** (-(mp.mp.dps - 30))
    for (lo, hi) in brackets:
        flo = dk(lo)
        while hi - lo > tol:
            mid = (lo + hi) / 2
            fm = dk(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        roots.append((lo + hi) / 2)
    m = len(roots)
    # inertia at arc midpoints
    mids = [((roots[i] + roots[(i + 1) % m] + (twopi if i == m - 1 else 0)) / 2) for i in range(m)]
    Nneg = [inertia_neg(hk(k)) for k in mids]
    # root i sits between arc i-1 and arc i
    right = sum(1 for i in range(m) if Nneg[i] - Nneg[i - 1] == -1)
    left = sum(1 for i in range(m) if Nneg[i] - Nneg[i - 1] == +1)
    d2 = min((roots[(i + 1) % m] - roots[i]) % twopi for i in range(m))
    sums = sorted((roots[i] + roots[j]) % twopi for i in range(m) for j in range(i, m))
    d4 = min((sums[(i + 1) % len(sums)] - sums[i]) % twopi for i in range(len(sums)))
    # argument principle on |z| = r for the polynomial p(z) = det(T z^2 + M z + T^dag)
    def winding(r, npts):
        tot = mp.mpf(0)
        prev = None
        for i in range(npts + 1):
            z = r * mp.expj(twopi * i / npts)
            v = mp.det(T * z * z + M * z + TH)
            a = mp.arg(v)
            if prev is not None:
                da = a - prev
                while da > mp.pi: da -= twopi
                while da < -mp.pi: da += twopi
                tot += da
            prev = a
        return tot / twopi
    out = dict(file=os.path.basename(fn), n=n, dps=mp.mp.dps, fgam=str(fgam), grid=grid,
               mu_antihermitian_over_T=mp.nstr(herm / scT, 5), det_imag_over_scale=mp.nstr(imag_max / max(abs(v) for v in vals), 5),
               field_over_gamma_min=mp.nstr(f / gam_min, 12),
               phase_T_over_pi=mp.nstr(th / mp.pi, 15), log10_gamma_min=mp.nstr(mp.log10(gam_min), 12),
               hermitian_part_positive=bool(gam_min > 0),
               sign_changes=m, expected=2 * n, right_movers=right, left_movers=left, Nneg_on_arcs=Nneg,
               roots=[mp.nstr(r, 25) for r in roots], root_separation=mp.nstr(d2, 12), pair_sum_separation=mp.nstr(d4, 12),
               pairs=len(sums), seconds_grid=round(time.time() - t0, 1))
    t0 = time.time()
    try:
        rin = mp.mpf("0.999"); rout = mp.mpf("1.001")
        npts = max(4 * grid, 4000)
        out["winding_inside_r0.999"] = mp.nstr(winding(rin, npts), 8)
        out["winding_inside_r1.001"] = mp.nstr(winding(rout, npts), 8)
    except Exception as e:
        out["winding_error"] = repr(e)
    out["seconds_winding"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    fn, fgam = sys.argv[1], sys.argv[2]
    grid = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    out = main(fn, fgam, grid)
    print(json.dumps(out, indent=1))
    tag = os.path.basename(fn).replace("blocks_", "fermi_").replace(".json", f"_f{str(fgam).replace(':', '_').replace('^', '')}.json")
    json.dump(out, open(os.path.join(HERE, tag), "w"), indent=1)
