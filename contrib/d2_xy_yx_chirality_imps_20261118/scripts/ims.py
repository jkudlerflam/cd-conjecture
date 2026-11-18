"""iMPS (iDMRG warm-up + single-site VUMPS) driver for the c=3/2 question.  2026-11-18, v3.

Models (all with Z2 parity prod Z conserved, 2-site unit cell):

 target, frame=rot : H = sum_j [ d (S+S- + h.c.) + (-1)^j (S+S+ + h.c.) + Z Z ] - h sum Z      (real)
                     unitarily equal (site j rotated by j pi/2 about z) to
 target, frame=orig: H = sum_j [ a X_jY_{j+1} + b Y_jX_{j+1} + Z Z ] - h sum Z,  a=(1+d)/2, b=(1-d)/2  (complex)
     + pinning  -eps sum_j (-1)^j kap_j,  kap_j = X_j Y_{j+2} - Y_j X_{j+2}   (frame-dependent sign, see below)
     In the rot frame H is real, so a real MPS is chirality-symmetric (A = complex conjugation there).
     kap_j(orig) = -kap_j(rot).  Everything printed is kap in the frame used; |kstag| is frame independent.

 ctrl : each MPS site carries two spins (A,B) (local dim 4; charge = parity of A only, so B's X is neutral):
        A = XXZ chain  -(tA/2)(XX+YY) + JzA Z Z - hA Z ;  B = TFIM  -J X X - g Z ;  NN in MPS order.
        decoupled  => c = 1 + 1/2 at g = J.  pinning  -eps sum X_B.  (unit cell 2 MPS sites, like the target)

usage: python ims.py model=target h=1.75 d=0.125 chis=32,48,... branch=pin|real [eps0=..] [hs=...] tag=...
"""
import sys, time, warnings
import numpy as np
import h5py
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinHalfSite, Site, GroupedSite
from tenpy.models.lattice import Lattice
from tenpy.linalg.charges import LegCharge
from tenpy.networks.mps import MPS, TransferMatrix
from tenpy.algorithms import dmrg
from tenpy.algorithms.vumps import SingleSiteVUMPSEngine
from tenpy.tools import hdf5_io

warnings.filterwarnings("ignore")


def _isingsite(chinfo):
    """spin-1/2 site whose states all carry charge 0 under the (sublattice-A) parity, so X is neutral."""
    leg = LegCharge.from_qflat(chinfo, [[0], [0]])
    return Site(leg, ["up", "down"], Sigmax=np.array([[0., 1.], [1., 0.]]), Sigmaz=np.diag([1., -1.]),
                Sp=np.array([[0., 1.], [0., 0.]]), Sm=np.array([[0., 0.], [1., 0.]]))


class Chain(CouplingMPOModel):
    def init_sites(self, p):
        sA = SpinHalfSite(conserve="parity")
        if p.get("kind", "target") == "target":
            return sA
        # control: one MPS site = (XXZ spin 'A', Ising spin 'B'); charge = parity of the A spins only.  As in the target,
        # the finite-chi boson sector may form a U(1)-breaking parity cat (odd-sector TM eigenvalue -> 1).
        return GroupedSite([sA, _isingsite(sA.leg.chinfo)], labels=["0", "1"], charges="same")

    def init_terms(self, p):
        kind = p.get("kind", "target")
        L = self.lat.N_sites
        eps = p.get("eps", 0.0)
        if kind == "target":
            d = p.get("d", 0.125); h = p.get("h", 1.75); frame = p.get("frame", "rot")
            for i in range(L):
                self.add_onsite_term(-h, i, "Sigmaz")
                self.add_coupling_term(1.0, i, i + 1, "Sigmaz", "Sigmaz")
                if frame == "rot":
                    s = 1.0 if i % 2 == 0 else -1.0
                    for o1, o2, c in (("Sp", "Sp", s), ("Sm", "Sm", s), ("Sp", "Sm", d), ("Sm", "Sp", d)):
                        self.add_coupling_term(c, i, i + 1, o1, o2)
                else:
                    a, b = (1 + d) / 2, (1 - d) / 2
                    self.add_coupling_term(a, i, i + 1, "Sigmax", "Sigmay")
                    self.add_coupling_term(b, i, i + 1, "Sigmay", "Sigmax")
                epsB = p.get("epsB", 0.0)
                if epsB != 0:  # -epsB (X_i X_{i+2} - Y_i Y_{i+2}) uniform in the rot frame (= staggered in orig): bond order B
                    self.add_coupling_term(-2 * epsB, i, i + 2, "Sp", "Sp"); self.add_coupling_term(-2 * epsB, i, i + 2, "Sm", "Sm")
                if eps != 0:   # -eps (-1)^i (X_i Y_{i+2} - Y_i X_{i+2}),  XY - YX = 2i(S+S- - S-S+)
                    s = -eps * (-1) ** i
                    self.add_coupling_term(2j * s, i, i + 2, "Sp", "Sm")
                    self.add_coupling_term(-2j * s, i, i + 2, "Sm", "Sp")
        else:
            tA = p.get("tA", 1.0); JzA = p.get("JzA", 0.0); hA = p.get("hA", 0.0)
            J = p.get("J", 1.0); g = p.get("g", 1.0)
            for i in range(L):
                self.add_coupling_term(-tA, i, i + 1, "Sp0", "Sm0"); self.add_coupling_term(-tA, i, i + 1, "Sm0", "Sp0")
                if JzA != 0: self.add_coupling_term(JzA, i, i + 1, "Sigmaz0", "Sigmaz0")
                self.add_onsite_term(-hA, i, "Sigmaz0")
                tB = p.get("tB", 0.0)
                if tB != 0:  # c=2 control: second XX chain on the B spins (its own U(1)/parity not conserved)
                    self.add_coupling_term(-tB, i, i + 1, "Sp1", "Sm1"); self.add_coupling_term(-tB, i, i + 1, "Sm1", "Sp1")
                if J != 0:
                    self.add_coupling_term(-J, i, i + 1, "Sigmax1", "Sigmax1")
                self.add_onsite_term(-g, i, "Sigmaz1")
                if eps != 0: self.add_onsite_term(-eps, i, "Sigmax1")


def model(P, eps=0.0, epsB=0.0):
    q = dict(P); q.update(L=2, bc_MPS="infinite", eps=eps, epsB=epsB)
    return Chain(q)


def idmrg(psi, M, chi, sweeps, mixer=True, amp=1e-5):
    par = {"trunc_params": {"chi_max": chi, "svd_min": 1e-14}, "mixer": mixer,
           "mixer_params": {"amplitude": amp, "decay": 1.5, "disable_after": max(sweeps // 2, 4)},
           "max_E_err": 1e-14, "max_S_err": 1e-10, "max_sweeps": sweeps, "min_sweeps": sweeps, "N_sweeps_check": 10,
           "norm_tol": 1e-3, "lanczos_params": {"N_min": 2, "N_max": 30, "E_tol": 1e-14, "P_tol": 1e-14}}
    eng = dmrg.TwoSiteDMRGEngine(psi, M, par)
    E, _ = eng.run()
    return E


LZ = dict(N_min=2, N_max=30, E_tol=1e-15, P_tol=1e-15)


def idmrg_conv(psi, M, P, chi, nmix=20, blk=20, maxblk=80, tolE=1e-11, tolS=5e-6, tolk=5e-6, hours=1e9, minblk=4):
    """Two-site iDMRG at fixed chi: nmix sweeps with mixer (bond expansion), then mixer-free blocks of `blk` sweeps
    until E, S and the order parameter change by less than the tolerances between blocks.
    (TeNPy single-site VUMPS was unstable for the target at chi >= 48 -- diverged from good iDMRG states.)"""
    t0 = time.time()
    E = idmrg(psi, M, chi, nmix) if nmix > 0 else np.nan
    prev = None; nb = 0; conv = False
    lz = {"N_min": 2, "N_max": 20, "E_tol": 1e-13, "P_tol": 1e-13}
    for nb in range(1, maxblk + 1):
        par = {"trunc_params": {"chi_max": chi, "svd_min": 1e-14}, "mixer": False, "max_sweeps": blk,
               "min_sweeps": blk, "N_sweeps_check": blk, "max_E_err": 1e-16, "max_S_err": 1e-16, "norm_tol": 1e9,
               "lanczos_params": lz}
        eng = dmrg.TwoSiteDMRGEngine(psi, M, par)
        E, _ = eng.run()
        S = float(np.mean(psi.entanglement_entropy())); k = float(measure(psi, P)["kstag"])
        if nb >= minblk and prev is not None and abs(E - prev[0]) < tolE and abs(S - prev[1]) < tolS and \
                abs(k - prev[2]) < tolk:
            conv = True
            break
        if prev is not None and (time.time() - t0) / 3600 > hours:
            break
        prev = (E, S, k)
    psi.canonical_form()
    E = float(np.mean(np.real(M.H_MPO.expectation_value(psi))))
    dE = abs(E - prev[0]) if prev else np.nan; dS = abs(S - prev[1]) if prev else np.nan
    return E, nmix + nb * blk, conv, dE, dS


def vumps(psi, M, max_sweeps, tol, max_hours, stol=1e-6):
    opts = {"mixer": False, "max_sweeps": max_sweeps, "min_sweeps": 30, "N_sweeps_check": 10, "max_E_err": 1e-14,
            "max_S_err": stol, "max_split_err": tol, "check_overlap": False, "norm_tol": 1e-14, "max_hours": max_hours,
            "lanczos_params": dict(LZ)}
    eng = SingleSiteVUMPSEngine(psi, M, opts)
    _, psi2 = eng.run()          # norm_tol tiny => returns sweep energy, skips the fragile env re-solve
    psi2.canonical_form()
    E = float(np.mean(np.real(M.H_MPO.expectation_value(psi2))))
    st = eng.sweep_stats
    return E, psi2, eng.sweeps, st["max_split_err"][-1], st["Delta_E"][-1], st["Delta_S"][-1]


def tm_spec(psi, q, n):
    try:
        T = TransferMatrix(psi, psi, charge_sector=([0] if q == 0 else [q]))
        ev, _ = T.eigenvectors(num_ev=n, which="LM")
        ev = np.array(ev)
        return ev[np.argsort(-np.abs(ev))]
    except Exception as e:  # noqa
        print("   tm_spec failed:", repr(e), flush=True)
        return np.zeros(0, complex)


def nnn(psi, i, o1, o2):
    return psi.expectation_value_term([(o1, i), (o2, i + 2)])


def measure(psi, P):
    o = {}
    o["S"] = np.array(psi.entanglement_entropy())
    o["Z"] = np.real(psi.expectation_value("Sigmaz" if P["kind"] == "target" else "Sigmaz0"))
    if P["kind"] == "target":
        kap = np.zeros(2); nx = np.zeros(2); ny = np.zeros(2)
        for i in range(2):
            spm, smp = nnn(psi, i, "Sp", "Sm"), nnn(psi, i, "Sm", "Sp")
            spp, smm = nnn(psi, i, "Sp", "Sp"), nnn(psi, i, "Sm", "Sm")
            kap[i] = np.real(2j * (spm - smp)); nx[i] = np.real(2 * (spp + smm)); ny[i] = np.real(-2j * (spp - smm))
        o["kap"] = kap; o["nx"] = nx; o["ny"] = ny
        o["kstag"] = (kap[0] - kap[1]) / 2; o["kuni"] = (kap[0] + kap[1]) / 2
        o["nxstag"] = (nx[0] - nx[1]) / 2; o["nxuni"] = (nx[0] + nx[1]) / 2
        o["zz"] = np.array([np.real(psi.expectation_value_term([("Sigmaz", i), ("Sigmaz", i + 1)])) for i in range(2)])
        o["pp"] = np.array([np.real(psi.expectation_value_term([("Sp", i), ("Sp", i + 1)]) +
                                    psi.expectation_value_term([("Sm", i), ("Sm", i + 1)])) for i in range(2)])
    else:
        o["mI"] = float(np.real(psi.expectation_value("Sigmax1", [0])[0]))
        o["kstag"] = o["mI"]
    return o


def correlators(psi, P, rmax):
    """connected staggered-chirality (target) or Ising-X (ctrl) correlator from site i0=0,1, and Z, S+S- correlators."""
    r = np.unique(np.concatenate([np.arange(3, 64), np.round(np.geomspace(64, rmax, 60)).astype(int)]))
    r = r[r <= rmax]
    out = {"r": r}
    if P["kind"] == "target":
        terms = [(("Sp", "Sm"), 2j), (("Sm", "Sp"), -2j)]
        for i0 in range(2):
            c = np.zeros(len(r), complex)
            for (a1, b1), c1 in terms:
                for (a2, b2), c2 in terms:
                    c += c1 * c2 * np.array(psi.term_correlation_function_right([(a1, 0), (b1, 2)], [(a2, 0), (b2, 2)],
                                                                                i_L=i0, j_R=list(i0 + r)))
            out["kk%d" % i0] = c
    else:
        out["kk1"] = np.array(psi.correlation_function("Sigmax1", "Sigmax1", sites1=[0], sites2=list(r)))[0]
        out["pm"] = np.array(psi.correlation_function("Sp0", "Sm0", sites1=[0], sites2=list(r)))[0]
        return out
    s2 = list(r)
    out["zz"] = np.array(psi.correlation_function("Sigmaz", "Sigmaz", sites1=[0], sites2=s2))[0]
    out["pm"] = np.array(psi.correlation_function("Sp", "Sm", sites1=[0], sites2=s2))[0]
    out["pp"] = np.array(psi.correlation_function("Sp", "Sp", sites1=[0], sites2=s2))[0]
    return out


def report(psi, E, P, chi, tag, t0, info, docorr, save=True):
    o = measure(psi, P)
    e0 = tm_spec(psi, 0, 16); e1 = tm_spec(psi, 1, 12)
    lam0 = e0[1:] if len(e0) > 1 else np.array([np.nan])
    xi0 = -2 / np.log(np.abs(lam0)); xi1 = -2 / np.log(np.abs(e1)) if len(e1) else np.array([np.nan])
    a0 = np.abs(lam0); a1 = np.abs(e1) if len(e1) else np.array([np.nan])
    good = np.concatenate([a0[a0 < 1 - 1e-7], a1[a1 < 1 - 1e-7]])
    xi = -2 / np.log(np.nanmax(good)) if len(good) else np.nan      # largest finite length (cat partners excluded)
    cat = int(np.sum(a1 > 1 - 1e-7))
    chiact = max(psi.chi)
    extra = ""
    if P["kind"] == "target":
        extra = "kuni=%+.1e nxs=%+.2e nxu=%+.5f dzz=%+.2e dpp=%+.2e" % (o["kuni"], o["nxstag"], o["nxuni"],
                                                                     o["zz"][0] - o["zz"][1], o["pp"][0] + o["pp"][1])
    print("RESULT %s cat=%d chi=%d chiact=%d E=%.13f S=%.8f,%.8f xi=%.3f xi0=%.3f xi1=%.3f kstag=%+.8f Z=%.7f,%.7f %s | %s t=%.0f" %
          (P["lab"], cat, chi, chiact, np.real(E), o["S"][0], o["S"][1], xi, xi0[0], xi1[0], o["kstag"], o["Z"][0], o["Z"][1],
           extra, info, time.time() - t0), flush=True)
    print("   TM0", " ".join("(%.6f,%+.4f)" % (-np.log(abs(x)) / 2, np.angle(x) / np.pi) for x in lam0), flush=True)
    print("   TM1", " ".join("(%.6f,%+.4f)" % (-np.log(abs(x)) / 2, np.angle(x) / np.pi) for x in e1), flush=True)
    dat = dict(E=np.real(E), chi=chi, chiact=chiact, tm0=e0, tm1=e1, spec0=np.array(psi.entanglement_spectrum()[0]),
               spec1=np.array(psi.entanglement_spectrum()[1]), **{k: v for k, v in o.items()},
               **{"P_" + k: v for k, v in P.items() if not isinstance(v, str)})
    fn = "out3/%s_chi%d%s" % (P["lab"], chi, tag)
    if save:
        np.savez(fn + ".npz", **dat)
    if docorr:
        try:
            rm = 6 * xi if np.isfinite(xi) else 200
            c = correlators(psi, P, int(min(max(rm, 200), 4000)))
            dat.update({"c_" + k: v for k, v in c.items()})
            if save:
                np.savez(fn + ".npz", **dat)
        except Exception as e:  # noqa
            print("   correlators failed:", repr(e), flush=True)
    if save:
        with h5py.File(fn + ".h5", "w") as f:
            hdf5_io.save_to_hdf5(f, {"psi": psi, "E": np.real(E)})
    return o


def params(kw, h=None):
    kind = kw.get("model", "target")
    if kind == "target":
        P = dict(kind="target", d=float(kw.get("d", 0.125)), h=float(h if h is not None else kw.get("h", 1.75)),
                 frame=kw.get("frame", "rot"))
        P["lab"] = "tg_d%.5f_h%.5f_%s" % (P["d"], P["h"], P["frame"])
    else:
        P = dict(kind="ctrl", tA=float(kw.get("tA", 1.0)), JzA=float(kw.get("JzA", 0.0)), hA=float(kw.get("hA", 0.0)),
                 J=float(kw.get("J", 1.0)), g=float(h if h is not None else kw.get("g", 1.0)), tB=float(kw.get("tB", 0.0)))
        P["lab"] = "ct_t%.3f_jz%.3f_hA%.3f_J%.3f_g%.5f" % (P["tA"], P["JzA"], P["hA"], P["J"], P["g"])
        if P["tB"] != 0:
            P["lab"] += "_tB%.3f" % P["tB"]
    return P


if __name__ == "__main__":
    from tenpy.tools.misc import setup_logging
    kw = dict(a.split("=", 1) for a in sys.argv[1:])
    chis = [int(x) for x in kw.get("chis", "32,48,64").split(",")]
    branch = kw.get("branch", "pin"); tag = kw.get("tag", "_" + branch)
    eps0 = float(kw.get("eps0", 0.02)); chi0 = int(kw.get("chi0", 16)); eps = float(kw.get("eps", 0.0))
    nd = int(kw.get("nd", 40)); nv = int(kw.get("nv", 3000)); vtol = float(kw.get("vtol", 1e-5))
    maxh = float(kw.get("maxh", 6.0)); docorr = int(kw.get("corr", 1)); init = kw.get("init", "")
    hs = [float(x) for x in kw["hs"].split(",")] if "hs" in kw else None
    seed = float(kw.get("seed", 0.0)); solver = kw.get("solver", "idmrg")
    if "lz" in kw:   # e.g. lz=12,1e-12
        a, b = kw["lz"].split(","); LZ.update(N_max=int(a), E_tol=float(b), P_tol=float(b))
    setup_logging(to_stdout=kw.get("log", "WARNING"))
    t0 = time.time()
    P = params(kw, hs[0] if hs else None)
    # ---- resume: skip (h, chi) stages whose .h5 already exists; restart from the last one
    import os
    todo = [(hh, c) for hh in (hs if hs else [None]) for c in (chis if hs is None else chis[:1])]
    done_files = []
    for hh, c in todo:
        fn = "out3/%s_chi%d%s.h5" % (params(kw, hh)["lab"], c, tag)
        if os.path.exists(fn):
            done_files.append(((hh, c), fn))
        else:
            break
    if done_files:
        init = done_files[-1][1]
        todo = todo[len(done_files):]
        print("# resume: %d stages done, continuing from %s" % (len(done_files), init), flush=True)
    if not todo:
        print("# all stages done", flush=True); open("out3/DONE%s_%s" % (tag, P["lab"]), "w").close(); sys.exit(0)
    import datetime
    def hours_left():
        try:
            return (float(os.environ["JOB_END_EPOCH"]) - time.time()) / 3600.0
        except Exception:
            return 1e9
    part = None
    if todo:
        pf = "out3/%s_chi%d%s.part.h5" % (params(kw, todo[0][0])["lab"], todo[0][1], tag)
        if os.path.exists(pf):
            init = pf; part = pf
            print("# resume from partial", pf, flush=True)
    if init:
        with h5py.File(init, "r") as f:
            psi = hdf5_io.load_from_hdf5(f)["psi"]
        print("# loaded", init, psi.chi, psi.dtype, flush=True)
    else:
        M0 = model(P)
        # target: all up.  ctrl: A sector at half filling (Sz conserved => filling fixed by the initial state), B up
        st = [0, 0] if P["kind"] == "target" else ["up_0 up_1", "down_0 up_1"]
        psi = MPS.from_product_state(M0.lat.mps_sites(), st, bc="infinite",
                                     dtype=(complex if (branch in ("pin", "cplx", "pinB") or P.get("frame") == "orig") else float))
        if branch == "pinB":
            E = idmrg(psi, model(P, 0.0, eps0), chi0, 60, amp=1e-2)
            report(psi, E, P, chi0, tag + "_pinned", t0, "epsB0=%g" % eps0, 0, save=False)
        elif branch == "pin":
            E = idmrg(psi, model(P, eps0), chi0, 60, amp=1e-2)
            report(psi, E, P, chi0, tag + "_pinned", t0, "eps0=%g" % eps0, 0, save=False)
        else:
            E = idmrg(psi, M0, chi0, 60, amp=1e-2)
    for hh, chi in todo:
        P = params(kw, hh)
        M = model(P, eps)
        if True:
            if part is None:
                if seed != 0:   # weak seed of the chirality (hysteresis scans): a few pinned sweeps, then release
                    idmrg(psi, model(P, seed), chi, 10)
                if solver == "vumps":
                    E = idmrg(psi, M, chi, nd)
            part = None
            tl = min(maxh, hours_left() - 0.25)
            if solver == "vumps":
                E, psi, ns, se, dE, dS = vumps(psi, M, nv, vtol, tl)
                ok = (se < vtol and abs(dS) < 1e-6)
            else:
                E, ns, ok, dE, dS = idmrg_conv(psi, M, P, chi, hours=tl, maxblk=int(kw.get("maxblk", 80)),
                                                nmix=int(kw.get("nmix", 20)))
                se = 0.0 if ok else 1.0
            if not ok and hours_left() < 0.4:
                pf = "out3/%s_chi%d%s.part.h5" % (P["lab"], chi, tag)
                with h5py.File(pf, "w") as f:
                    hdf5_io.save_to_hdf5(f, {"psi": psi})
                print("# out of time in stage chi=%d (split=%.1e): saved %s" % (chi, se, pf), flush=True)
                sys.exit(3)
            report(psi, E, P, chi, tag, t0, "%s sweeps=%d conv=%s dE=%.1e dS=%.1e eps=%g dtype=%s" %
                   (solver, ns, ok, dE, dS, eps, psi.dtype.__name__ if hasattr(psi.dtype, "__name__") else psi.dtype), docorr)
    open("out3/DONE%s_%s" % (tag, params(kw, hs[0] if hs else None)["lab"]), "w").close()
    print("# done %.0fs" % (time.time() - t0), flush=True)
