"""Inspect saved real-space correlators (connected) in out3 npz files.
usage: python corr.py file.npz ...   prints staggered chirality correlator G(r) = (-1)^r C_kk(r) (connected),
its local log-slope, and the dominant wavevector of C_kk near pi from a windowed Fourier transform."""
import sys
import numpy as np

for fn in sys.argv[1:]:
    d = np.load(fn)
    if "c_r" not in d.files:
        print(fn, "no correlators"); continue
    r = d["c_r"]; kap = d["kap"]
    for i0 in (0,):
        c = np.real(d["c_kk%d" % i0]) - kap[i0] * kap[(i0 + r) % 2]
        G = c * (-1.0) ** r
        print("==", fn.split("/")[-1], "i0=%d kstag=%+.5f xi=?" % (i0, float(d["kstag"])))
        for rr in (4, 8, 16, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024):
            j = np.searchsorted(r, rr)
            if j < len(r):
                print("   r=%5d C=%+.3e  (-1)^r C=%+.3e" % (r[j], c[j], G[j]))
        # dominant wavevector: uniform grid part r < 64 plus fine sampling; use r in [3, 63]
        m = r < 64
        rs, cs = r[m], c[m]
        qs = np.linspace(0.8, 1.0, 401) * np.pi
        F = np.array([np.abs(np.sum(cs * np.exp(1j * q * rs))) for q in qs])
        print("   peak |FT| (r<64) at q/pi = %.4f" % (qs[np.argmax(F)] / np.pi))
        # sign pattern of G for r in 3..63: count sign changes of G
        sg = np.sign(G[m]); nch = np.sum(sg[1:] != sg[:-1])
        print("   sign changes of (-1)^r C for 3<=r<64: %d" % nch)
