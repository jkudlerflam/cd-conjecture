"""Correlators on stored open-chain states, in the format of codex/work_C193/ckpt_meas (C193b) so its frozen
estimators (eta_b, pair exponent, K_c, K_s) run unchanged.  2026-11-23.

Particles are down spins (as in the STATUS candidate, field favouring up): b^dag = Sm, n = (1 - Z)/2.
  G[i, j]   = <Sm_i Sp_j>                 i <= L/2 + 5, all j > i            (one-body, centred pairs used)
  D[d, i0]  = <Sm_i0 Sm_i0+d Sp_j Sp_j+d> - <Sm_i0 Sm_i0+d><Sp_j Sp_j+d>, d = 2, 3, i0 = L/2-3..L/2, j > i0+d
  C[i0]     = <n_i0 n_j> - <n_i0><n_j>,  same i0
States are in the rotated frame of ims.py (site j rotated by j pi/2 about z).  Relative to the STATUS frame (a uniform
rotation) this multiplies Sm_i Sp_j by a phase linear in i - j, i.e. shifts momenta, which the C193 estimators detect
from the data.

usage: python meas_obc.py file=out_obc/<...>.h5 out=meas/<tag>.json      |   python meas_obc.py test=1 (ED check)
"""
import sys, os, json, time
import numpy as np
import h5py
from tenpy.tools import hdf5_io

DS = (2, 3)


def measure(psi, L, imax=None):
    t0 = time.time()
    REFS = list(range(L // 2 - 3, L // 2 + 1))
    imax = min(L // 2 + 5, L - 2) if imax is None else imax
    Zs = np.real(psi.expectation_value("Sigmaz"))
    dens = (1 - Zs) / 2
    G = psi.correlation_function("Sm", "Sp", sites1=list(range(0, imax + 1)), sites2=list(range(L)))
    print("# G done %.0fs" % (time.time() - t0), flush=True)
    ZZ = np.real(psi.correlation_function("Sigmaz", "Sigmaz", sites1=REFS, sites2=list(range(L))))
    C = {i0: [(ZZ[a, j] - Zs[i0] * Zs[j]) / 4 for j in range(i0 + 1, L)] for a, i0 in enumerate(REFS)}
    D, mm, pp = {}, {}, {}
    for d in DS:
        pp[d] = {j: complex(psi.expectation_value_term([("Sp", j), ("Sp", j + d)])) for j in range(L - d)}
        for i0 in REFS:
            mm[(d, i0)] = complex(psi.expectation_value_term([("Sm", i0), ("Sm", i0 + d)]))
            js = list(range(i0 + d + 1, L - d))
            if not js:
                continue
            v = psi.term_correlation_function_right([("Sm", 0), ("Sm", d)], [("Sp", 0), ("Sp", d)], i_L=i0, j_R=js)
            D[(d, i0)] = [complex(x) - mm[(d, i0)] * pp[d][j] for x, j in zip(v, js)]
    print("# pairs done %.0fs" % (time.time() - t0), flush=True)
    return dict(L=L, G=G, C=C, D=D, dens=dens, imax=imax)


def to_json(m, tag, E, extra):
    L = m["L"]; G = m["G"]
    out = {"tag": tag, "cfg": {"P": extra.get("par", 0), "L": L}, "check": {"E_stored": E, "dE": 0.0, "dSbar6": 0.0,
                                                                              "dGb_window": 0.0}, "extra": extra}
    mm = {"L": L, "dens": [float(x) for x in m["dens"]]}
    mm["Gb_rows"] = {str(i): [[float(np.real(G[i, j])) for j in range(i + 1, L)], [float(np.imag(G[i, j])) for j in range(i + 1, L)]]
                     for i in range(m["imax"] + 1)}
    mm["Cnn"] = {str(i0): [float(x) for x in v] for i0, v in m["C"].items()}
    mm["Dpair"] = {"%d_%d" % k: [[v.real, v.imag] for v in vals] for k, vals in m["D"].items()}
    mm["mm"] = {"%d_%d" % k: [0.0, 0.0] for k in m["D"]}      # D is stored already connected
    mm["pp"] = {str(d): {str(j): [0.0, 0.0] for j in range(L)} for d in DS}
    out["meas"] = mm
    return out


def test():
    """ED check at L = 10 on the target chain: G, C and D against dense vectors."""
    import ims, ofd
    from tenpy.algorithms import dmrg
    L = 10
    P = ims.params({"model": "target", "d": "0.1875", "h": "1.783"})
    M = ofd.fmodel(P, L); psi = ofd.init_state(M, P, L, 0, 1)
    dmrg.TwoSiteDMRGEngine(psi, M, {"trunc_params": {"chi_max": 64, "svd_min": 1e-14}, "mixer": True,
                                    "max_sweeps": 40, "min_sweeps": 20, "max_E_err": 1e-14}).run()
    m = measure(psi, L)
    v = psi.get_theta(0, L).to_ndarray().reshape(-1)     # TeNPy basis: index 0 = up
    sp = np.array([[0, 1], [0, 0.]]); sm = sp.T; z = np.diag([1., -1.])
    def op(ops):
        M_ = np.array([[1.]])
        for k in range(L):
            M_ = np.kron(M_, ops.get(k, np.eye(2)))
        return M_
    ev = lambda ops: np.vdot(v, op(ops) @ v)
    err = 0
    for (i, j) in [(0, 9), (2, 7), (4, 5), (1, 8)]:
        err = max(err, abs(ev({i: sm, j: sp}) - m["G"][i, j]))
    for i0 in m["C"]:
        for t, x in enumerate(m["C"][i0][:3]):
            j = i0 + 1 + t
            n = lambda k: (np.eye(2) - z) / 2
            ex = ev({i0: n(i0), j: n(j)}) - ev({i0: n(i0)}) * ev({j: n(j)})
            err = max(err, abs(ex - x))
    for (d, i0), vals in m["D"].items():
        for t, x in enumerate(vals):
            j = i0 + d + 1 + t
            if j + d >= L:
                break
            if j in (i0, i0 + d) or j + d in (i0, i0 + d):
                continue
            ex = ev({i0: sm, i0 + d: sm, j: sp, j + d: sp}) - ev({i0: sm, i0 + d: sm}) * ev({j: sp, j + d: sp})
            err = max(err, abs(ex - x))
    print("TEST meas_obc L=%d max |MPS - dense| over G, C, D = %.1e" % (L, err), flush=True)


if __name__ == "__main__":
    from tenpy.tools.misc import setup_logging
    kw = dict(a.split("=", 1) for a in sys.argv[1:])
    setup_logging(to_stdout="WARNING")
    if kw.get("test"):
        test(); sys.exit(0)
    fn = kw["file"]
    with h5py.File(fn, "r") as f:
        psi = hdf5_io.load_from_hdf5(f)["psi"]
    L = psi.L
    tag = os.path.basename(fn)[:-3]
    npz = fn[:-3] + ".npz"
    E = float(np.load(npz)["E"]) if os.path.exists(npz) else float("nan")
    m = measure(psi, L)
    os.makedirs(os.path.dirname(kw["out"]) or ".", exist_ok=True)
    json.dump(to_json(m, tag, E, {"file": fn, "chi": int(max(psi.chi))}), open(kw["out"], "w"))
    print("# wrote", kw["out"], flush=True)
