"""C188 follow-up at the corner point (coordinator, 2026-11-09). See README.md.
Usage: python3 run085.py L a,b,g
The sweep loop is copied from dmrg3.dmrg (same start, same local solve), with a log line per sweep,
a full record at the end of each chi stage, and an early stop in the final stage."""
import sys, os, json, time
import numpy as np
import scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from dmrg3 import _lupdate, _rupdate, _heff, entropy
from dmrg_tau import mpo_tau, cfit, site_expect, bond_expect, couplings, X, iY, Z

OLD = [16, 32, 64, 128, 160, 160]     # C188's schedule, one sweep each
MID = [240, 240]
FINAL, MIN_FINAL, MAX_FINAL = 320, 2, 6
TOL_E, TOL_S = 1e-8, 1e-4             # early stop, tighter than K1 and K2
CUTOFF, LTOL, d = 1e-10, 1e-9, 2
STAGE_MIN = 160                       # full records from this chi up


def local(M, j, Lenv, Renv, Ws, chi, dt, right):
    th = np.tensordot(M[j], M[j + 1], axes=([2], [0])); shape = th.shape
    op = _heff(Lenv[j], Ws[j], Ws[j + 1], Renv[j + 2], shape, dt)
    if th.size <= 16:
        Hd = np.array([op.matvec(e) for e in np.eye(th.size, dtype=dt)]).T
        w, v = np.linalg.eigh(Hd)
    else:
        w, v = spla.eigsh(op, k=1, which="SA", v0=th.reshape(-1), tol=LTOL, ncv=min(20, th.size - 1))
    E = float(np.real(w[0])); th = v[:, 0].reshape(shape[0] * d, d * shape[3])
    U, s, Vh = np.linalg.svd(th, full_matrices=False)
    keep = min(chi, max(1, int(np.sum(s > CUTOFF * s[0]))))
    tr = float(np.sum(s[keep:] ** 2) / np.sum(s ** 2))
    s = s[:keep] / np.linalg.norm(s[:keep])
    if right:
        M[j] = U[:, :keep].reshape(shape[0], d, keep)
        M[j + 1] = (s[:, None] * Vh[:keep]).reshape(keep, d, shape[3])
        Lenv[j + 1] = _lupdate(Lenv[j], M[j], Ws[j])
    else:
        M[j + 1] = Vh[:keep].reshape(keep, d, shape[3])
        M[j] = (U[:, :keep] * s[None, :]).reshape(shape[0], d, keep)
        Renv[j + 1] = _rupdate(Renv[j + 2], M[j + 1], Ws[j + 1])
    return E, tr


def entropies(M, chi):
    """Schmidt entropies of every cut, on a copy (the state M has its centre at site 0)."""
    M = list(M); L = len(M); S = np.zeros(L - 1)
    for j in range(L - 1):
        th = np.tensordot(M[j], M[j + 1], axes=([2], [0])); shape = th.shape
        U, s, Vh = np.linalg.svd(th.reshape(shape[0] * d, d * shape[3]), full_matrices=False)
        keep = min(chi, max(1, int(np.sum(s > CUTOFF * s[0]))))
        s = s[:keep] / np.linalg.norm(s[:keep]); S[j] = entropy(s)
        M[j] = U[:, :keep].reshape(shape[0], d, keep)
        M[j + 1] = (s[:, None] * Vh[:keep]).reshape(keep, d, shape[3])
    return S, M


def main():
    L = int(sys.argv[1]); a, b, g = [float(x) for x in sys.argv[2].split(",")]
    fs = os.path.join(HERE, f"sweeps_L{L}.jsonl"); ff = os.path.join(HERE, f"stages_L{L}.jsonl")
    Ws = mpo_tau(a, b, g, L); dt = np.result_type(*[W.dtype for W in Ws], float)
    rng = np.random.default_rng(0)
    D = [1]
    for j in range(L):
        D.append(min(4, d ** (j + 1), d ** (L - j - 1)))
    M = [rng.standard_normal((D[j], d, D[j + 1])).astype(dt) for j in range(L)]
    for j in range(L - 1, 0, -1):
        A = M[j].reshape(M[j].shape[0], -1)
        Q, R = np.linalg.qr(A.T.conj())
        M[j] = Q.T.conj().reshape(-1, d, M[j].shape[2])
        M[j - 1] = np.tensordot(M[j - 1], R.T.conj(), axes=([2], [0]))
    M[0] /= np.linalg.norm(M[0])
    Lenv = [None] * (L + 1); Renv = [None] * (L + 1)
    Lenv[0] = np.ones((1, 1, 1), dtype=dt); Renv[L] = np.ones((1, 1, 1), dtype=dt)
    for j in range(L - 1, 0, -1):
        Renv[j] = _rupdate(Renv[j + 1], M[j], Ws[j])
    sched = OLD + MID + [FINAL] * MAX_FINAL
    t0 = time.time(); Eprev = None; Sprev = None; nfinal = 0; chis_done = []
    for sweep, chi in enumerate(sched):
        tr = 0.0
        for j in range(L - 1):
            E, t = local(M, j, Lenv, Renv, Ws, chi, dt, True); tr = max(tr, t)
        for j in range(L - 2, -1, -1):
            E, t = local(M, j, Lenv, Renv, Ws, chi, dt, False); tr = max(tr, t)
        chis_done.append(chi)
        S, Mc = entropies(M, chi)
        S4 = float(np.mean(S[L // 2 - 3:L // 2 + 1]))
        dE = None if Eprev is None else E - Eprev
        dS = None if Sprev is None else S4 - Sprev
        rec = dict(L=L, sweep=sweep, chi=chi, E=E, dE=dE, trunc=tr, S_mid4=S4, dS=dS,
                   c_fit=cfit(S, L, 0.25, 0.75), sec=round(time.time() - t0, 1))
        with open(fs, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"L={L} sweep {sweep} chi={chi} E={E:.10f} dE={dE} trunc={tr:.1e} S4={S4:.5f} dS={dS} ({rec['sec']}s)", flush=True)
        nfinal += chi == FINAL
        stop = (chi == FINAL and nfinal >= MIN_FINAL and abs(dE) <= TOL_E and abs(dS) <= TOL_S)
        last_of_stage = stop or sweep == len(sched) - 1 or sched[sweep + 1] != chi
        if last_of_stage and chi >= STAGE_MIN:
            xx = bond_expect(Mc, X, X); yy = -bond_expect(Mc, iY, iY); zz = bond_expect(Mc, Z, Z)
            e = np.array([np.dot(couplings(a, b, g, j), (xx[j], yy[j], zz[j])) for j in range(L - 1)])
            mx = site_expect(Mc, X); mz = site_expect(Mc, Z)
            full = dict(a=a, b=b, g=g, L=L, chis=list(chis_done), E=E, E_per_site=E / L, trunc=tr, dE_last=dE, dS_last=dS,
                        c_fit=cfit(S, L, 0.25, 0.75), c_fit_wide=cfit(S, L, 0.15, 0.85), S=[round(float(s), 7) for s in S],
                        S_mid4=S4, e=[round(float(x), 8) for x in e], mx=[round(float(x), 5) for x in mx],
                        mz=[round(float(x), 5) for x in mz], mmax=float(max(np.max(np.abs(mx)), np.max(np.abs(mz)))),
                        E_bonds=float(np.sum(e)), sec=round(time.time() - t0, 1))
            with open(ff, "a") as f:
                f.write(json.dumps(full) + "\n")
            print(f"L={L} STAGE chi={chi} E={E:.10f} sum(e)={full['E_bonds']:.10f} S4={S4:.5f} c_fit={np.round(full['c_fit'], 3).tolist()} mmax={full['mmax']:.4f}", flush=True)
        if stop:
            break
        Eprev, Sprev = E, S4
    print(f"L={L} done", flush=True)


if __name__ == "__main__":
    main()
