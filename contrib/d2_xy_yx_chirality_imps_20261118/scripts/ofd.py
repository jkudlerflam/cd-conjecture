"""Open-chain finite two-site DMRG for the c ~ 2 region (second method).  2026-11-19.

Same Chain model as ims.py (target in the rotated real frame, or the decoupled controls), open boundaries, length L.
Parity is conserved; the run stays in the sector of its initial state (even by default).  Real dtype throughout
(target rot frame and controls are real; the c ~ 2 region is chirality symmetric).

Per chi stage: sweeps until |dE| < tolE and |dS_mid| < tolS between sweeps (TeNPy checks), then saves
  out_obc/<lab>_L<L>_chi<chi>.npz : E, S(l) for every bond, Z(l), trunc error, sweeps, time
and the state <...>.h5 (the previous stage's h5 is deleted to save disk).  Resumable; exits 3 when out of time.

usage: python ofd.py model=target d=0.1875 h=1.783 L=256 chis=64,128,256,384 [par=0|1] [seed=1]
       python ofd.py model=ctrl tA=1 hA=0.3 tB=1 J=0 g=0.3 L=...    (c=2 control, incommensurate)
       python ofd.py test=1 ...   ED check at small L
"""
import sys, os, time, warnings
import numpy as np
import h5py
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
from tenpy.tools import hdf5_io
import ims

warnings.filterwarnings("ignore")


def fmodel(P, L):
    q = dict(P); q.update(L=L, bc_MPS="finite", bc_x="open", eps=0.0, epsB=0.0)
    return ims.Chain(q)


def init_state(M, P, L, par, seed):
    rng = np.random.default_rng(seed)
    if P["kind"] == "target":
        st = list(rng.choice(["up", "down"], size=L, p=[0.6, 0.4]))
        ndown = sum(s == "down" for s in st)
        if ndown % 2 != par:
            st[0] = "down" if st[0] == "up" else "up"
    else:
        a = list(rng.choice(["up", "down"], size=L)); b = list(rng.choice(["up", "down"], size=L))
        if sum(s == "down" for s in a) % 2 != par:
            a[0] = "down" if a[0] == "up" else "up"
        st = ["%s_0 %s_1" % (x, y) for x, y in zip(a, b)]
    return MPS.from_product_state(M.lat.mps_sites(), st, bc="finite", dtype=float)


def hours_left():
    try:
        return (float(os.environ["JOB_END_EPOCH"]) - time.time()) / 3600.0
    except Exception:
        return 1e9


def main(kw):
    P = ims.params(kw)
    L = int(kw.get("L", 64)); par = int(kw.get("par", 0)); seed = int(kw.get("seed", 1))
    chis = [int(x) for x in kw.get("chis", "64,128,256").split(",")]
    tolE = float(kw.get("tolE", 1e-9)); tolS = float(kw.get("tolS", 1e-5))
    M = fmodel(P, L)
    if kw.get("test"):
        from tenpy.algorithms.exact_diag import ExactDiag
        ed = ExactDiag(M); ed.build_full_H_from_mpo(); ed.full_diagonalization()
        Es = np.sort(ed.E)
        psi = init_state(M, P, L, par, seed)
        eng = dmrg.TwoSiteDMRGEngine(psi, M, {"trunc_params": {"chi_max": 64, "svd_min": 1e-14}, "mixer": True,
                                              "max_sweeps": 40, "max_E_err": 1e-13})
        E, _ = eng.run()
        print("TEST L=%d ED E0=%.12f E1=%.12f  DMRG(par=%d)=%.12f" % (L, Es[0], Es[1], par, E), flush=True)
        return
    lab = "%s_L%d_p%d" % (P["lab"], L, par)
    os.makedirs("out_obc", exist_ok=True)
    t0 = time.time()
    psi = None; done = []
    for c in chis:
        if os.path.exists("out_obc/%s_chi%d.npz" % (lab, c)):
            done.append(c)
    todo = [c for c in chis if c not in done]
    if not todo:
        print("# all done", flush=True); open("out_obc/DONE_%s" % lab, "w").close(); return
    for c in [todo[0]] + list(reversed(done)):
        for fn in ("out_obc/%s_chi%d.part.h5" % (lab, c), "out_obc/%s_chi%d.h5" % (lab, c)):
            if psi is None and os.path.exists(fn):
                with h5py.File(fn, "r") as f:
                    psi = hdf5_io.load_from_hdf5(f)["psi"]
                print("# resume from", fn, flush=True)
    first = psi is None
    if first:
        psi = init_state(M, P, L, par, seed)
    prev_h5 = None
    for chi in todo:
        lim = hours_left() - 0.25
        if lim < 0.2:
            sys.exit(3)
        mix = {"amplitude": 1e-4 if first else 1e-6, "decay": 2.0, "disable_after": 12 if first else 4}
        params = {"trunc_params": {"chi_max": chi, "svd_min": 1e-12}, "mixer": True, "mixer_params": mix,
                  "max_E_err": tolE, "max_S_err": tolS, "max_trunc_err": 1.0, "min_sweeps": 14 if first else 4,
                  "max_sweeps": int(kw.get("maxsw", 40)), "max_hours": lim,
                  "lanczos_params": {"N_min": 2, "N_max": 20, "E_tol": 1e-13, "P_tol": 1e-13}}
        first = False
        ts = time.time()
        eng = dmrg.TwoSiteDMRGEngine(psi, M, params)
        E, _ = eng.run()
        st = eng.sweep_stats
        nsw = len(st["sweep"]); dE = abs(st["E"][-1] - st["E"][-2]) if nsw > 1 else np.nan
        S = np.array(psi.entanglement_entropy())
        dS = abs(st["S"][-1] - st["S"][-2]) if nsw > 1 else np.nan
        tr = float(st["max_trunc_err"][-1])
        conv = bool(dE < tolE * max(abs(E), 1) * 10 and dS < tolS * 10)
        if not conv and nsw < params["max_sweeps"] and hours_left() < 0.6:
            fn = "out_obc/%s_chi%d.part.h5" % (lab, chi)
            with h5py.File(fn, "w") as f:
                hdf5_io.save_to_hdf5(f, {"psi": psi})
            print("# out of time at chi=%d after %d sweeps; saved %s" % (chi, nsw, fn), flush=True)
            sys.exit(3)
        Z = np.real(psi.expectation_value("Sigmaz" if P["kind"] == "target" else "Sigmaz0"))
        np.savez("out_obc/%s_chi%d.npz" % (lab, chi), E=E, S=S, Z=Z, trunc=tr, nsw=nsw, dE=dE, dS=dS, conv=conv,
                 chi_real=max(psi.chi), t=time.time() - ts, L=L)
        fn = "out_obc/%s_chi%d.h5" % (lab, chi)
        with h5py.File(fn, "w") as f:
            hdf5_io.save_to_hdf5(f, {"psi": psi})
        for old in (prev_h5, "out_obc/%s_chi%d.part.h5" % (lab, chi)):
            if old and os.path.exists(old):
                os.remove(old)
        prev_h5 = fn
        print("RESULT %s chi=%d conv=%d E=%.12f Smid=%.6f trunc=%.1e sweeps=%d dE=%.1e dS=%.1e chi_real=%d t=%.0fs" %
              (lab, chi, conv, E, S[L // 2 - 1], tr, nsw, dE, dS, max(psi.chi), time.time() - ts), flush=True)
    open("out_obc/DONE_%s" % lab, "w").close()


if __name__ == "__main__":
    from tenpy.tools.misc import setup_logging
    kw = dict(a.split("=", 1) for a in sys.argv[1:])
    setup_logging(to_stdout=kw.get("log", "WARNING"))
    for LL in kw.get("L", "64").split(","):   # L=64,96,128 runs the sizes in sequence
        kw["L"] = LL
        main(kw)
