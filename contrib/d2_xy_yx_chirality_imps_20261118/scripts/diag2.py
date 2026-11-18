"""Solver diagnostic at one (h, chi): mode=idmrg (long two-site iDMRG, mixer off after warm-up) or mode=vumps2
(TwoSiteVUMPS).  Prints E, S, kstag, xi every block."""
import sys, time
import numpy as np
import h5py
from tenpy.tools import hdf5_io
from tenpy.algorithms import dmrg
from tenpy.algorithms.vumps import TwoSiteVUMPSEngine
import ims

kw = dict(a.split("=", 1) for a in sys.argv[1:])
P = ims.params(kw); M = ims.model(P)
with h5py.File(kw["init"], "r") as f:
    psi = hdf5_io.load_from_hdf5(f)["psi"]
chi = int(kw["chi"]); mode = kw.get("mode", "idmrg"); nblk = int(kw.get("nblk", 30)); blk = int(kw.get("blk", 20))
t0 = time.time()
lz = {"N_min": 2, "N_max": 20, "E_tol": 1e-13, "P_tol": 1e-13}


def show(tag, E):
    o = ims.measure(psi, P)
    e0 = ims.tm_spec(psi, 0, 4); e1 = ims.tm_spec(psi, 1, 3)
    xi0 = -2 / np.log(abs(e0[1])); xi1 = -2 / np.log(abs(e1[0]))
    print("%s E=%.12f S=%.5f,%.5f k=%+.5f B=%+.5f Z=%.5f xi0=%.1f xi1=%.1f norm=%.1e t=%.0f" %
          (tag, E, o["S"][0], o["S"][1], o["kstag"], o["nxuni"], o["Z"].mean(), xi0, xi1,
           np.linalg.norm(psi.norm_test()), time.time() - t0), flush=True)


if mode == "idmrg":
    E = ims.idmrg(psi, M, chi, int(kw.get("nmix", 20)))
    show("mix", E)
    for b in range(nblk):
        par = {"trunc_params": {"chi_max": chi, "svd_min": 1e-14}, "mixer": False, "max_sweeps": blk,
               "min_sweeps": blk, "N_sweeps_check": blk, "max_E_err": 1e-16, "max_S_err": 1e-16, "norm_tol": 1e9,
               "lanczos_params": lz}
        eng = dmrg.TwoSiteDMRGEngine(psi, M, par)
        E, _ = eng.run()
        show("blk%d" % b, E)
else:
    E = ims.idmrg(psi, M, chi, int(kw.get("nmix", 20)))
    show("mix", E)
    for b in range(nblk):
        opts = {"trunc_params": {"chi_max": chi, "svd_min": 1e-14}, "mixer": False, "max_sweeps": blk,
                "min_sweeps": blk, "N_sweeps_check": blk, "max_E_err": 1e-16, "max_S_err": 1e-16, "max_split_err": 1e-16,
                "check_overlap": False, "norm_tol": 1e9, "lanczos_params": lz}
        eng = TwoSiteVUMPSEngine(psi, M, opts)
        _, psi = eng.run()
        psi.canonical_form()
        E = float(np.mean(np.real(M.H_MPO.expectation_value(psi))))
        show("blk%d split=%.1e" % (b, eng.sweep_stats["max_split_err"][-1]), E)
