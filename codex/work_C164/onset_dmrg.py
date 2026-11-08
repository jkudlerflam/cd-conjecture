"""C164 Task 3: direct test past a two-valley onset (Claude Opus agent, 2026-11-08).

Input: a tuned two-valley bond v_c (qtune endpoint, E_A ~ E_B small, ring L = 14 band).
Directions (from central differences of the tracked valley energies at L1 = 12, L2 = 14):
  u_mu  : minimal-norm direction with dE_A = dE_B = -1 (both valleys go down together),
  u_lam : minimal-norm direction with dE_A = +1/2, dE_B = -1/2 (valley splitting).
Version 2 (bisection). Coarse grid, then bisection in lam at depth m = 2 on the valley-weight difference,
then L = 64 and 128 at the switching point lam*, where an immiscible gas switches valleys abruptly (c = 1 on
both sides) while a miscible one has both valleys populated over a window of width of order the pseudospin
exchange (c = 2). Valley weights: w_v = max S_+(q) over |q - s k_v| < 0.3, with the helicity sign s fixed by
whichever choice aligns better with the peaks.
Grid (coarse): v(mu, lam) = v_c + mu u_mu + lam u_lam with mu = (E_A + E_B)/2 + m delta, m = 1, 2, 4, so that both
valleys sit at -m delta to first order (dilute gas, chemical potential m delta; the linear extrapolation from
the ring band is only a guide, hence three depths), and lam scanning the splitting E_A - E_B over
[-2 delta, 2 delta] (nl values).
Each grid point: open-chain DMRG (C155 dmrg.py, C144-validated), entanglement fits (dmrg_c.profile_fits),
connected correlations C(r) = sum_a <s^a_i s^a_{i+r}> - <s^a_i><s^a_{i+r}> from i0 = L/4, structure factor
S(q) = sum_r cos(q r) C(r) w(r) (Hann window), and its peaks.
Pre-fixed c = 2 criteria (spec Task 3): c fits > 1.5 not falling from L = 128 to 256, truncation < 1e-7,
two incommensurate S(q) peaks near k_A and k_B, FES slope > 1.5 at chi >= 256. This script does the L = 64
and 128 stages; points with c_prof(128) > 1.3 are written to promote.jsonl for L = 256 and FES.
usage: onset_dmrg.py tuned.jsonl index out.json [delta_frac] [L list]
"""
import sys, json, time
import numpy as np
import est, qopt, qtune, dmrg as D
from dmrg_c import chain_bond, profile_fits

PA = est.PAULI


def directions(v, kA, kB, h=1e-3):
    Jm = np.zeros((2, 9))
    for i in range(9):
        dv = np.zeros(9); dv[i] = h
        Ep, _, _ = qtune.track(v + dv, kA, kB); Em, _, _ = qtune.track(v - dv, kA, kB)
        Jm[:, i] = (np.array(Ep) - np.array(Em)) / (2 * h)
    P = np.linalg.pinv(Jm)
    return P @ np.array([-1.0, -1.0]), P @ np.array([0.5, -0.5]), Jm


def correlations(M, i0, rmax):
    """C^{ab}(r) = <s^a_i0 s^b_{i0+r}> - <s^a_i0><s^b_{i0+r}>, r = 0..rmax, and <s^a_j> on the window.
    M from dmrg.dmrg: left-canonical up to the last site (which carries the norm)."""
    L = len(M)
    R = [None] * (L + 1); R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        A = M[j]; R[j] = np.einsum('asc,cd,bsd->ab', np.conj(A), R[j + 1], A)
    Lf = np.ones((1, 1))
    for j in range(i0):
        A = M[j]; Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A)
    norm = np.einsum('ab,ab->', Lf, R[i0]).real
    jmax = min(L - 1, i0 + rmax)
    n = jmax - i0 + 1
    mag = np.zeros((n, 3)); C = np.zeros((n, 3, 3))
    TI = Lf; Ts = None
    for j in range(i0, jmax + 1):
        A = M[j]
        for a in range(3):
            mag[j - i0, a] = np.einsum('cd,cd->', np.einsum('ab,asc,st,btd->cd', TI, np.conj(A), PA[a], A), R[j + 1]).real / norm
        if j == i0:
            for a in range(3):
                for b in range(3):
                    C[0, a, b] = np.einsum('cd,cd->', np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), PA[a] @ PA[b], A), R[j + 1]).real / norm
            C[0] -= np.outer(mag[0], mag[0])
            Ts = [np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), PA[a], A) for a in range(3)]
        else:
            for a in range(3):
                for b in range(3):
                    C[j - i0, a, b] = np.einsum('cd,cd->', np.einsum('ab,asc,st,btd->cd', Ts[a], np.conj(A), PA[b], A), R[j + 1]).real / norm
                Ts[a] = np.einsum('ab,asc,bsd->cd', Ts[a], np.conj(A), A)
            C[j - i0] -= np.outer(mag[0], mag[j - i0])
        TI = np.einsum('ab,asc,bsd->cd', TI, np.conj(A), A)
    return C, mag


def sq_peaks(C, nbg, nq=720):
    """Helicity-resolved transverse structure factor S_+(q) = e_+^dag S(q) e_+, e_+ = (u + i w)/sqrt 2 in the frame
    of the local moment nbg, S^{ab}(q) = sum_r e^{-iqr} C^{ab}(r) w(r) with C^{ab}(-r) = C^{ba}(r) (Hann window).
    For real C one has S_-(q) = S_+(-q), so S_+ on [0, 2 pi) carries the full momentum information and, unlike
    the eigenvalues of S(q), distinguishes q from -q."""
    n = len(C); r = np.arange(n); w = 0.5 * (1 + np.cos(np.pi * r / n))
    nb = np.asarray(nbg, float)
    nb = nb / np.linalg.norm(nb) if np.linalg.norm(nb) > 1e-6 else np.array([0, 0, 1.0])
    a0 = np.array([1.0, 0, 0]) if abs(nb[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(nb, a0); u /= np.linalg.norm(u); wv = np.cross(nb, u)
    ep = (u + 1j * wv) / np.sqrt(2)
    qs = np.linspace(0, 2 * np.pi, nq, endpoint=False)
    C0 = 0.5 * (C[0] + C[0].T)
    S = []
    for q in qs:
        ph = np.exp(-1j * q * r[1:])[:, None, None] * w[1:, None, None]
        Sq = C0 + np.sum(ph * C[1:], 0) + np.sum(np.conj(ph) * np.transpose(C[1:], (0, 2, 1)), 0)
        S.append(np.real(np.conj(ep) @ Sq @ ep))
    S = np.array(S)
    loc = [i for i in range(nq) if S[i] >= S[i - 1] and S[i] >= S[(i + 1) % nq]]
    loc.sort(key=lambda i: -S[i])
    return qs, S, [(float(qs[i]), float(S[i])) for i in loc[:4]]


def run_point(v, L, chis, seed=0):
    h, edge = chain_bond(v)
    t0 = time.time()
    r = D.dmrg(h, L, chis=chis, seed=seed, edge_field=edge)
    S = r['S']
    cp, cs, mid = profile_fits(S, L)
    i0 = L // 4; rmax = L // 2
    C, mag = correlations(r['M'], i0, rmax)
    qs, Sq, peaks = sq_peaks(C, mag[min(len(mag) - 1, L // 4)])
    return dict(L=L, E=r['E'], c_prof=cp, c_sm=cs, Smid=mid, trunc=r['hist'][-1][2], S=[float(x) for x in S],
                C=C.tolist(), Sq=Sq[::4].tolist(), Sq_full=Sq.tolist(), mag0=mag[0].tolist(), peaks=peaks, t=time.time() - t0)


def weights(run, kA, kB):
    qs = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    out = {}
    Sq = np.array(run['Sq_full'])
    for sgn in (1, -1):
        wv = []
        for k in (kA, kB):
            d = np.abs((qs - sgn * k + np.pi) % (2 * np.pi) - np.pi)
            wv.append(float(Sq[d < 0.3].max()))
        out[sgn] = wv
    return out


if __name__ == '__main__':
    import os
    fin, idx, fout = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    dfrac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    Ls = [int(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else [64, 128]
    nbis = int(sys.argv[6]) if len(sys.argv) > 6 else 7
    T = [json.loads(l) for l in open(fin)]
    rec = T[idx]
    last = rec['hist'][-1]
    v = np.array(last['v']); inf = last['info']
    kA, kB = inf['kA'], inf['kB']
    umu, ulam, Jm = directions(v, kA, kB)
    delta = dfrac * inf['bw']
    g = 0.5 * (inf['EA'] + inf['EB'])
    out = dict(version=2, src=rec['src'], v_c=v.tolist(), info=inf, Jm=Jm.tolist(), delta=delta, points=[])
    L0 = Ls[0]
    chis0 = (16, 32, 64, 96, 128, 128)

    def point(mfac, lam, L, chis):
        vv = v + (g + mfac * delta) * umu + lam * ulam
        vv = vv / np.linalg.norm(vv)
        rr = run_point(vv, L, chis)
        w = weights(rr, kA, kB)
        pt = dict(mu=float(mfac * delta), lam=float(lam), L=L, v=vv.tolist(), w=w,
                  **{k: rr[k] for k in ('E', 'c_prof', 'c_sm', 'Smid', 'trunc', 'peaks', 'mag0', 't', 'S', 'Sq', 'C')})
        out['points'].append(pt)
        json.dump(out, open(fout + '.tmp', 'w')); os.replace(fout + '.tmp', fout)
        print(f"mu={mfac * delta:.4f} lam={lam:+.4f} L={L} c_prof={rr['c_prof']:.3f} c_sm={rr['c_sm']:.3f} trunc={rr['trunc']:.1e} "
              f"w(+)={np.round(w[1], 2).tolist()} w(-)={np.round(w[-1], 2).tolist()} peaks={[(round(q, 3), round(s_, 2)) for q, s_ in rr['peaks'][:3]]} {rr['t']:.0f}s", flush=True)
        return pt

    off = -(inf['EA'] - inf['EB'])
    lo, hi = off - 2 * delta, off + 2 * delta                       # lo: valley A lower; hi: valley B lower
    coarse = [point(2.0, lam, L0, chis0) for lam in (lo, off, hi)]
    # helicity convention: the sign whose valley windows hold the larger peaks at the two ends
    sc = {s: coarse[0]['w'][s][0] + coarse[2]['w'][s][1] for s in (1, -1)}
    sgn = max(sc, key=sc.get); out['sign'] = sgn
    f = lambda p: p['w'][sgn][0] - p['w'][sgn][1]
    fl, fh = f(coarse[0]), f(coarse[2])
    out['bracket'] = [fl, fh]
    if fl * fh < 0:
        a, b, fa = lo, hi, fl
        for it in range(nbis):
            m = 0.5 * (a + b); pm = point(2.0, m, L0, chis0); fm = f(pm)
            if fm * fa > 0:
                a, fa = m, fm
            else:
                b = m
        lstar = 0.5 * (a + b)
    else:
        lstar = off
    out['lam_star'] = float(lstar)
    for L in Ls:
        chis = chis0 if L <= 64 else (16, 32, 64, 128, 192, 256, 256)
        point(2.0, lstar, L, chis)
    point(4.0, lstar, L0, chis0)
