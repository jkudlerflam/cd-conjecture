"""C171: independent operator counts for the T / K_I candidates of C159 Proposition C.

Counts the relevant (dimension < 2) operators that are invariant under one-site
translation T1 (unitary) and the antiunitary A (Theta in T), for products of
compact bosons and discrete CFT sectors, as a function of the Luttinger parameters.

Conventions as in C159: V_{p,n} = exp(i(p theta + 2 n phi)), dimension
p^2/(4K) + n^2 K, conformal spin p*n.
T1:  theta -> s theta + b, phi -> s phi + a      (s = +1 shift, s = -1 charge conjugation)
A:   theta -> t_th theta + d, phi -> t_ph phi + c (antiunitary; exactly one of t = -1)
The count of invariant Hermitian real directions equals the complex dimension of the
joint +1 eigenspace of T1 and L = dagger o A (both complex-linear), with
  T1 V_{p,n} = e^{i(p b + 2 n a)} V_{s p, s n},
  L  V_{p,n} = e^{i(p d + 2 n c)} V_{t_th p, t_ph n}.
Currents d_x phi, d_x theta are included (dimension 1, spin +-1) so that couplings like
sigma d_x phi are seen. Spin-0 relevant invariants are tunings (C154 4.1). Relevant
invariants of nonzero spin are listed separately (C154/C159 count sigma J as a tuning,
but uniform currents alone as density/pitch shifts).
"""
import itertools, cmath, math
import numpy as np

ALPHA = 0.7 * math.sqrt(2.0)          # generic incommensurate offset
OMEGA = cmath.exp(2j * math.pi / 3)


def boson(K, s, a, b, tth, tph, c, d, name, pmax=4, nmax=3):
    ops = [("1",)]
    for p in range(-pmax, pmax + 1):
        for n in range(-nmax, nmax + 1):
            if (p, n) != (0, 0):
                ops.append(("V", p, n))
    ops += [("dphi",), ("dth",)]

    def dim(o):
        if o[0] == "1":
            return 0.0
        if o[0] == "V":
            return o[1] ** 2 / (4 * K) + o[2] ** 2 * K
        return 1.0

    def spin(o):
        if o[0] == "V":
            return o[1] * o[2]
        return 1 if o[0] in ("dphi", "dth") else 0

    def T1(o):
        if o[0] == "1":
            return 1.0, o
        if o[0] == "V":
            p, n = o[1], o[2]
            return cmath.exp(1j * (p * b + 2 * n * a)), ("V", s * p, s * n)
        return float(s), o

    def L(o):
        if o[0] == "1":
            return 1.0, o
        if o[0] == "V":
            p, n = o[1], o[2]
            return cmath.exp(1j * (p * d + 2 * n * c)), ("V", tth * p, tph * n)
        return float(tph if o[0] == "dphi" else tth), o

    def label(o):
        if o[0] == "V":
            return f"V{name}({o[1]},{o[2]})"
        return {"1": "", "dphi": f"dphi{name}", "dth": f"dth{name}"}[o[0]]

    return dict(ops=ops, dim=dim, spin=spin, T1=T1, L=L, label=label)


def ising(tau, eta):
    """Ising: sigma with T1 eigenvalue tau (momentum 0 or pi), L eigenvalue eta."""
    ops = [("1",), ("s",), ("e",)]
    dims = {"1": 0.0, "s": 1 / 8, "e": 1.0}
    return dict(ops=ops, dim=lambda o: dims[o[0]], spin=lambda o: 0,
                T1=lambda o: ((tau if o[0] == "s" else 1.0), o),
                L=lambda o: ((eta if o[0] == "s" else 1.0), o),
                label=lambda o: {"1": "", "s": "sigma", "e": "eps"}[o[0]])


def potts():
    """Three-state Potts with T1 acting as the Z3 generator (period-3 order)."""
    ch = {"1": 0, "s": 1, "sb": -1, "Z": 1, "Zb": -1, "e": 0, "X": 0}
    dims = {"1": 0, "s": 2 / 15, "sb": 2 / 15, "Z": 4 / 3, "Zb": 4 / 3, "e": 4 / 5, "X": 14 / 5}
    ops = [(k,) for k in ch]
    return dict(ops=ops, dim=lambda o: dims[o[0]], spin=lambda o: 0,
                T1=lambda o: (OMEGA ** ch[o[0]], o),
                L=lambda o: (1.0, o),   # A: sigma -> sigma^dagger, so dagger o A fixes sigma
                label=lambda o: {"1": "", "s": "sigma", "sb": "sigma*", "Z": "Z", "Zb": "Z*",
                                 "e": "eps", "X": "X"}[o[0]])


def count(sectors, cut=2.0 - 1e-9, keep_spin=True):
    prods = list(itertools.product(*[s["ops"] for s in sectors]))
    prods = [o for o in prods if not all(x[0] == "1" for x in o)]
    dim = lambda o: sum(s["dim"](x) for s, x in zip(sectors, o))
    rel = [o for o in prods if dim(o) < cut]
    idx = {o: i for i, o in enumerate(rel)}

    def mat(key):
        M = np.zeros((len(rel), len(rel)), complex)
        for o in rel:
            ph, img = 1.0, []
            for s, x in zip(sectors, o):
                f, y = s[key](x)
                ph *= f
                img.append(y)
            img = tuple(img)
            if img not in idx:      # image has the same dimension, so it is in rel
                raise RuntimeError("not closed")
            M[idx[img], idx[o]] += ph
        return M

    T1, L = mat("T1"), mat("L")
    comm = np.abs(T1 @ L - L @ T1).max()
    A = np.vstack([T1 - np.eye(len(rel)), L - np.eye(len(rel))])
    u, sv, vh = np.linalg.svd(A)
    null = vh[np.sum(sv > 1e-8):].conj().T
    # read off which operators the invariant directions touch
    out = []
    used = np.abs(null).sum(axis=1) > 1e-8
    for i, o in enumerate(rel):
        if used[i]:
            sp = sum(s["spin"](x) for s, x in zip(sectors, o))
            lab = "*".join(l for l in (s["label"](x) for s, x in zip(sectors, o)) if l)
            out.append((round(dim(o), 3), sp, lab))
    return null.shape[1], comm, sorted(set(out))


def tunings(sectors):
    """Number of invariant relevant directions, split by spin-0 versus spin != 0."""
    n, comm, ops = count(sectors)
    s0 = [o for o in ops if o[1] == 0]
    s1 = [o for o in ops if o[1] != 0 and not o[2].startswith(("dphi", "dth"))]
    return n, comm, s0, s1


# Host: the anomalous half-floating boson of T (C159 conventions):
# T1: phi + pi/2, theta + pi + alpha; Theta: phi -> -phi + pi/2, theta -> theta + pi.
def host(K, name="1"):
    return boson(K, +1, math.pi / 2, math.pi + ALPHA, +1, -1, math.pi / 2, math.pi, name)


def spin0_invariant_count(sectors):
    """Real dimension of the invariant space of relevant spin-0 operators (the tunings)."""
    prods = [o for o in itertools.product(*[s["ops"] for s in sectors])
             if not all(x[0] == "1" for x in o)]
    dim = lambda o: sum(s["dim"](x) for s, x in zip(sectors, o))
    spin = lambda o: sum(s["spin"](x) for s, x in zip(sectors, o))
    cur = lambda o: any(x[0] in ("dphi", "dth") for x in o)
    rel = [o for o in prods if dim(o) < 2 - 1e-9 and spin(o) == 0 and not cur(o)]
    idx = {o: i for i, o in enumerate(rel)}
    def mat(key):
        M = np.zeros((len(rel), len(rel)), complex)
        for o in rel:
            ph, img = 1.0, []
            for s, x in zip(sectors, o):
                f, y = s[key](x); ph *= f; img.append(y)
            M[idx[tuple(img)], idx[o]] += ph
        return M
    if not rel:
        return 0, []
    T1, L = mat("T1"), mat("L")
    A = np.vstack([T1 - np.eye(len(rel)), L - np.eye(len(rel))])
    sv = np.linalg.svd(A, compute_uv=False)
    n = len(rel) - np.sum(sv > 1e-8)
    _, _, vh = np.linalg.svd(A)
    null = vh[len(rel) - n:].conj().T if n else np.zeros((len(rel), 0))
    used = np.abs(null).sum(axis=1) > 1e-8 if n else np.zeros(len(rel), bool)
    labs = sorted({(round(dim(o), 3), "*".join(l for l in (s["label"](x) for s, x in zip(sectors, o)) if l))
                   for i, o in enumerate(rel) if used[i]})
    return int(n), labs


def main():
    Ks = [0.6, 0.9, 1.2, 1.6, 1.85, 1.9, 2.1, 3.0]
    print("== S1 host alone (T): expect only cos 4phi for K < 1/2, nothing for K > 1/2")
    for K in [0.4] + Ks:
        print(f"  K={K:4}: k={spin0_invariant_count([host(K)])}")

    print("\n== S2 Ising (sigma at pi) x host: C159 row, expect k=2 for K<15/8, k=1 above")
    for eta in (+1, -1):
        for K in Ks:
            k, labs = spin0_invariant_count([ising(-1, eta), host(K)])
            print(f"  eta={eta:+d} K={K:4}: k={k} {labs}")

    print("\n== S2' Ising (uniform, Theta-odd) x host: C159 row, expect eps + sigma*dphi (spin 1)")
    for K in (0.9, 1.6):
        n, comm, s0, s1 = tunings([ising(+1, -1), host(K)])
        print(f"  K={K}: spin-0 {s0}; spin!=0 relevant {s1}; [T1,L]={comm:.1e}")

    print("\n== S3 Potts (T1 = Z3 generator) x host: NOT in C159's table")
    for K in Ks:
        k, labs = spin0_invariant_count([potts(), host(K)])
        print(f"  K={K:4}: k={k} {labs}")

    print("\n== S4 host x boson 2 with T1 acting by charge conjugation (NOT in C159's table)")
    grid = [0, math.pi / 2, math.pi, 3 * math.pi / 2]
    for K2 in (0.1, 0.6, 1.0, 1.9, 2.5):
        for K1 in (1.2, 1.6, 1.85, 1.9):
            best = None
            for a2, b2, c2, d2 in itertools.product(grid, repeat=4):
                for tth, tph in ((-1, +1), (+1, -1)):
                    # A^2 = 1 on local operators
                    if tth == -1 and abs(cmath.exp(4j * c2) - 1) > 1e-9: continue
                    if tph == -1 and abs(cmath.exp(2j * d2) - 1) > 1e-9: continue
                    b2s = boson(K2, -1, a2, b2, tth, tph, c2, d2, "2")
                    secs = [host(K1), b2s]
                    _, comm, _ = count([b2s])
                    if comm > 1e-9: continue
                    k, labs = spin0_invariant_count(secs)
                    if best is None or k < best[0]:
                        best = (k, labs, (a2, b2, c2, d2, tth, tph))
            print(f"  K1={K1:4} K2={K2:4}: min k={best[0]}  {best[1]}")

    print("\n== S5 repulsive host x half-floating boson 2 whose theta is Theta-reversed (theta2 at pi)")
    for K1, K2 in ((0.9, 0.2), (0.9, 0.3), (0.7, 0.15), (0.7, 0.2)):
        b2s = boson(K2, +1, 0.37, math.pi, -1, +1, 0.0, 0.0, "2")   # density free (a2 generic)
        k, labs = spin0_invariant_count([host(K1), b2s])
        print(f"  K1={K1} K2={K2}: x1+x2={K1 + 1 / (4 * K2):.2f}  k={k} {labs}")

    print("\n== S6 Lemma L: composite boson 2 with theta2 shift n(pi+alpha) + delta (both half-floating)")
    for n, delta in ((2, 0.0), (2, math.pi), (3, 0.0)):
        for K1, K2 in ((1.2, 0.8), (1.5, 1.0), (1.5, 1.5), (1.8, 1.5), (2.2, 1.5)):
            b1 = math.pi + ALPHA
            b2s = boson(K2, +1, 0.0, n * b1 + delta, +1, -1, 0.0, (n * math.pi) % (2 * math.pi), "2")
            k, labs = spin0_invariant_count([host(K1), b2s])
            print(f"  n={n} delta={delta:.2f} K1={K1} K2={K2}: k={k} {labs}")


if __name__ == "__main__":
    main()
