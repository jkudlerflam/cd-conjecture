"""C188 (Claude Opus 5.5): driver for ringk (matrix-free momentum Lanczos) on the cyclic stratum
h = a XZ + b YX + g ZY.  Usage:
  python3 ringk_run.py L a,b,g m1,m2,... niter out.jsonl [exe]
For each sector m, runs ringk, removes Lanczos ghosts (Cullum-Willoughby) and records the lowest
converged levels."""
import sys, os, json, subprocess, time, tempfile
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]], complex); SZ = np.diag([1., -1.]).astype(complex)


def bond_abg(a, b, g):
    return a * np.kron(SX, SZ) + b * np.kron(SY, SX) + g * np.kron(SZ, SY)


def levels(al, be, nlev=4, tol=1e-9):
    """Eigenvalues of the Lanczos tridiagonal with Cullum-Willoughby ghost removal and a convergence check."""
    n = len(al)
    def ev(al, be):
        T = np.diag(al) + np.diag(be[:len(al) - 1], 1) + np.diag(be[:len(al) - 1], -1)
        return np.linalg.eigvalsh(T)
    w = ev(al, be)
    wh = ev(al[1:], be[1:]) if n > 2 else np.array([])
    out = []
    i = 0
    while i < len(w) and len(out) < nlev + 2:
        x = w[i]
        mult = np.sum(np.abs(w - x) < tol * max(1, abs(x)) * 1e3)
        spurious = mult == 1 and len(wh) and np.min(np.abs(wh - x)) < tol * max(1, abs(x)) * 10
        if not spurious and (not out or abs(x - out[-1]) > 1e-8 * max(1, abs(x))):
            out.append(float(x))
        i += max(1, mult)
    # convergence: compare with a shorter run
    w2 = ev(al[:n - 30], be[:n - 30]) if n > 60 else w
    conv = [float(np.min(np.abs(w2 - x))) for x in out]
    return out[:nlev], conv[:nlev]


if __name__ == "__main__":
    L = int(sys.argv[1]); a, b, g = [float(x) for x in sys.argv[2].split(",")]
    ms = [int(x) % L for x in sys.argv[3].split(",")]; niter = int(sys.argv[4]); fout = sys.argv[5]
    exe = sys.argv[6] if len(sys.argv) > 6 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "ringk")
    h = bond_abg(a, b, g)
    td = tempfile.mkdtemp(prefix="ringk_")
    hf = os.path.join(td, "h.txt")
    with open(hf, "w") as f:
        for o in range(4):
            for i in range(4):
                f.write(f"{h[o, i].real:.17g} {h[o, i].imag:.17g}\n")
    for m in ms:
        t0 = time.time()
        of = os.path.join(td, f"lz_{m}.txt")
        subprocess.run([exe, str(L), str(m), str(niter), "7", hf, of], check=True, stderr=subprocess.DEVNULL)
        dat = np.loadtxt(of, skiprows=1, ndmin=2)
        dim = int(open(of).readline().split()[2])
        al, be = dat[:, 0], dat[:, 1]
        lev, conv = levels(al, be)
        rec = dict(a=a, b=b, g=g, L=L, m=m, dim=dim, niter=len(al), levels=lev, conv=conv, sec=round(time.time() - t0, 1))
        with open(fout, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
