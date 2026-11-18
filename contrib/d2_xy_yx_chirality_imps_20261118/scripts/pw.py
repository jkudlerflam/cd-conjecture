"""Plane-wave (quasiparticle ansatz) dispersions on a converged iMPS from ims.py.
usage: python pw.py file=out3/....h5 model=target h=.. [d=..] [frame=rot] ks=0.01,0.02,... nev=4 q=0,1
ks are SITE momenta in units of pi; the engine gets p = 2*k*pi (2-site cell).  Prints the lowest nev excitation
energies in each parity sector q."""
import sys, time
import numpy as np
import h5py
from tenpy.tools import hdf5_io
from tenpy.networks.uniform_mps import UniformMPS
from tenpy.algorithms.plane_wave_excitation import PlaneWaveExcitationEngine
import ims

kw = dict(a.split("=", 1) for a in sys.argv[1:])
P = ims.params(kw)
M = ims.model(P)
with h5py.File(kw["file"], "r") as f:
    psi = hdf5_io.load_from_hdf5(f)["psi"]
psi.canonical_form()
u = UniformMPS.from_MPS(psi)
ks = [float(x) for x in kw.get("ks", "0,0.01,0.02,0.04").split(",")]
nev = int(kw.get("nev", 4)); qs = [int(x) for x in kw.get("q", "0,1").split(",")]
eng = PlaneWaveExcitationEngine(u, M, {"lanczos_params": {"N_max": int(kw.get("nmax", 60)), "N_min": 10,
                                                          "E_tol": 1e-10, "reortho": True}})
t0 = time.time()
print("# file", kw["file"], "chi", max(psi.chi), flush=True)
for q in qs:
    for k in ks:
        Es, _, N = eng.run(2 * np.pi * k, qtotal_change=[q], num_ev=nev)
        Es = np.sort(np.real(Es))
        print("PW q=%d k=%.5f  " % (q, k) + " ".join("%.6f" % e for e in Es) + "  N=%d t=%.0f" % (N, time.time() - t0),
              flush=True)
