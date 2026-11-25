"""Quartet (charge-four) correlators on stored open-chain states, for the k = 0 quartet exponent (2026-11-24, Claude Opus 5.5).

Particles are down spins, b^dag = Sm, as in meas_obc.py.  For a pair separation d and a pair-pair offset s,
    P^(d)_i    = Sm_i Sm_{i+d}                                   (charge two)
    Q^(d,s)_i  = Sm_i Sm_{i+d} Sm_{i+s} Sm_{i+s+d}               (charge four, two pairs)
Only parity commutes with H, so <Q> != 0 and both the normal <Q_i Q_j^dag> and the anomalous <Q_i Q_j> are nonzero.
For each reference site i0 and every j with the right operator entirely to the right of the left one we store the
connected correlators
    QQd = <Q_i0 Q_j^dag> - <Q_i0><Q_j^dag>,  QdQ = <Q_i0^dag Q_j> - ...,  QQa = <Q_i0 Q_j> - ...,  QdQd = <Q_i0^dag Q_j^dag> - ...
the same four for P^(d), the cross correlators QZ = <Q_i0 Z_j> - <Q_i0><Z_j> and PZ, and ZZ.
States are in the rotated frame of ims.py (site j rotated by j pi/2 about z).  A product of four Sm picks up the phase
e^{-i pi (4 i + offsets)/2}, which is independent of i, so momentum zero for Q is the same in both frames.  For P the
frame shifts the momentum by pi (meas_obc.py, kpair.py).

usage: python meas_quartet.py file=out_obc/<...>.h5 out=<tag>.json [refs=a,b,..] [ds=2:4,2:5,3:6]   |   test=1 (ED check)
"""
import sys, os, json, time
import numpy as np

DS_DEFAULT = [(2, 4), (2, 5), (3, 6), (2, 1)]


def qterm(op, d, s):
    return [(op, 0), (op, d), (op, s), (op, s + d)]


def pterm(op, d):
    return [(op, 0), (op, d)]


def span(term):
    return max(o for _, o in term)


def measure(psi, L, refs, ds):
    t0 = time.time()
    Zs = np.real(psi.expectation_value("Sigmaz"))
    out = {"L": L, "Z": Zs.tolist(), "refs": refs, "ds": ds, "Q": {}, "P": {}, "ZZ": {}}
    ZZ = np.real(psi.correlation_function("Sigmaz", "Sigmaz", sites1=refs, sites2=list(range(L))))
    for a, i0 in enumerate(refs):
        out["ZZ"][str(i0)] = [float(ZZ[a, j] - Zs[i0] * Zs[j]) for j in range(i0 + 1, L)]

    def block(mk, w):
        """mk(op) -> term; w = span.  Returns expectation lists and connected correlators for every ref."""
        tm, tp = mk("Sm"), mk("Sp")
        ev = {"m": [complex(psi.expectation_value_term([(o, j + k) for o, k in tm])) for j in range(L - w)],
              "p": [complex(psi.expectation_value_term([(o, j + k) for o, k in tp])) for j in range(L - w)]}
        res = {"ev_m": [[z.real, z.imag] for z in ev["m"]]}
        for i0 in refs:
            js = list(range(i0 + w + 1, L - w))
            if not js:
                continue
            r = {}
            for lab, tl, tr, a, b in (("md", tm, tp, "m", "p"), ("dm", tp, tm, "p", "m"),
                                      ("mm", tm, tm, "m", "m"), ("dd", tp, tp, "p", "p")):
                v = psi.term_correlation_function_right(tl, tr, i_L=i0, j_R=js)
                c = [complex(x) - ev[a][i0] * ev[b][j] for x, j in zip(v, js)]
                r[lab] = [[z.real, z.imag] for z in c]
            vz = psi.term_correlation_function_right(tm, [("Sigmaz", 0)], i_L=i0, j_R=js)
            r["mz"] = [[z.real, z.imag] for z in (complex(x) - ev["m"][i0] * Zs[j] for x, j in zip(vz, js))]
            r["j0"] = js[0]
            res[str(i0)] = r
        return res

    for d, s in ds:
        out["Q"]["%d_%d" % (d, s)] = block(lambda op: qterm(op, d, s), s + d)
        print("# Q d=%d s=%d done %.0fs" % (d, s, time.time() - t0), flush=True)
    for d in sorted({d for d, _ in ds}):
        out["P"][str(d)] = block(lambda op: pterm(op, d), d)
        print("# P d=%d done %.0fs" % (d, time.time() - t0), flush=True)
    return out


def test():
    """ED check at L = 14 on the target chain: every stored correlator against dense vectors."""
    import ims, ofd
    from tenpy.algorithms import dmrg
    L = 14
    P = ims.params({"model": "target", "d": "0.1875", "h": "1.783"})
    M = ofd.fmodel(P, L); psi = ofd.init_state(M, P, L, 0, 1)
    dmrg.TwoSiteDMRGEngine(psi, M, {"trunc_params": {"chi_max": 128, "svd_min": 1e-14}, "mixer": True,
                                    "max_sweeps": 40, "min_sweeps": 20, "max_E_err": 1e-14}).run()
    ds = [(2, 4), (2, 1), (1, 3)]
    refs = [0, 1]
    m = measure(psi, L, refs, ds)
    v = psi.get_theta(0, L).to_ndarray().reshape(-1)
    # With conserve='parity' the TeNPy site basis is (down, up): index 0 = down.  (Checked: <Z> from TeNPy and from the
    # dense vector agree only in this basis; correlators with an even number of spin flips cannot tell.)
    sp = np.array([[0, 0], [1, 0.]]); sm = sp.T.copy(); z = np.diag([-1., 1.])
    print("site basis", psi.sites[0].state_labels)
    O = {"Sm": sm, "Sp": sp, "Sigmaz": z}

    def ev(ops):   # ops: list of (name, site)
        w = v.reshape([2] * L)
        for name, k in ops:
            w = np.moveaxis(np.tensordot(O[name], w, axes=([1], [k])), 0, k)
        return np.vdot(v, w.reshape(-1))

    err = 0.0
    for i0 in refs:
        for t, x in enumerate(m["ZZ"][str(i0)]):
            j = i0 + 1 + t
            err = max(err, abs(ev([("Sigmaz", i0), ("Sigmaz", j)]) - ev([("Sigmaz", i0)]) * ev([("Sigmaz", j)]) - x))
    checked = 0; dbg = 0
    print("ZZ err", err)
    for kind, items in (("Q", [((d, s), qterm, s + d) for d, s in ds]), ("P", [((d,), pterm, d) for d in sorted({d for d, _ in ds})])):
        for key, mk, w in items:
            blk = m[kind]["_".join(map(str, key))]
            for i0 in refs:
                if str(i0) not in blk:
                    continue
                r = blk[str(i0)]
                for lab, a, b in (("md", "Sm", "Sp"), ("dm", "Sp", "Sm"), ("mm", "Sm", "Sm"), ("dd", "Sp", "Sp"), ("mz", "Sm", None)):
                    for t, (re, im) in enumerate(r[lab]):
                        j = r["j0"] + t
                        tl = [(o, i0 + k) for o, k in mk(a, *key)]
                        tr = [("Sigmaz", j)] if b is None else [(o, j + k) for o, k in mk(b, *key)]
                        ex = ev(tl + tr) - ev(tl) * ev(tr)
                        e1 = abs(ex - (re + 1j * im))
                        if e1 > 1e-8 and dbg < 12:
                            print("  mismatch", kind, key, i0, lab, j, ex, re + 1j * im, "raw", ev(tl + tr), ev(tl), ev(tr)); dbg += 1
                        err = max(err, e1); checked += 1
    print("TEST meas_quartet L=%d: %d correlators, max |MPS - dense| = %.1e" % (L, checked, err), flush=True)


if __name__ == "__main__":
    import h5py
    from tenpy.tools import hdf5_io
    from tenpy.tools.misc import setup_logging
    kw = dict(a.split("=", 1) for a in sys.argv[1:])
    setup_logging(to_stdout="WARNING")
    if kw.get("test"):
        test(); sys.exit(0)
    fn = kw["file"]
    with h5py.File(fn, "r") as f:
        psi = hdf5_io.load_from_hdf5(f)["psi"]
    L = psi.L
    refs = [int(x) for x in kw["refs"].split(",")] if "refs" in kw else [L // 2 - 17, L // 2 - 8, L // 2]
    ds = [tuple(int(y) for y in x.split(":")) for x in kw["ds"].split(",")] if "ds" in kw else DS_DEFAULT
    m = measure(psi, L, refs, ds)
    m.update(file=fn, chi=int(max(psi.chi)))
    json.dump(m, open(kw["out"], "w"))
    print("# wrote", kw["out"], flush=True)
