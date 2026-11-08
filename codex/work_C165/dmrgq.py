#!/usr/bin/env python3
"""C165 open-chain DMRG driver (Claude Opus agent, 2026-11-08), built on the validated dmrg.py (C144 + C155
edge_field, made d-generic in dmrg_d.py by inferring d from the bond) and dmrg_c.profile_fits (C155).

usage: dmrgq.py MODEL L CHI out.json
  MODEL = 'fam:a=0.9,Dz=1.2,g=0.3,...'   bond of classical.fam (field split over the two sites, f/2 at the ends)
        = 'tfi'  (-XX - Z, c = 1/2), 'xx' (XX + YY, c = 1), 'isxx' (rung of a critical TFI qubit and an XX qubit,
                  d = 4, c = 3/2; the calibration outside B_2)
Records: energy, truncation, S(l), c_prof / c_sm fits (Calabrese-Cardy with parity term, window L/8..7L/8),
the bulk staggered moment and uniform moment (sites L/4..3L/4), and connected correlators C_aa(r) from site L/4,
with the structure-factor peak positions.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import dmrg_d as D
from dmrg_c import profile_fits
from est import bond_from_JB, PAULI
from classical import fam, nfloat

I2 = np.eye(2, dtype=complex)


def model(name):
    if name.startswith('fam:'):
        kw = eval('dict(' + name[4:] + ')')
        J, B = fam(**kw) if 'nfloat' not in kw else nfloat(**kw['nfloat'])
        return bond_from_JB(J, B), sum(B[a] * PAULI[a] for a in range(3)) / 2, 2
    if name == 'tfi':
        J = np.diag([-1., 0, 0]); B = np.array([0, 0, -1.])
        return bond_from_JB(J, B), sum(B[a] * PAULI[a] for a in range(3)) / 2, 2
    if name == 'xx':
        return bond_from_JB(np.diag([1., 1, 0]), np.zeros(3)), None, 2
    if name == 'isxx':
        hi = bond_from_JB(np.diag([-1., 0, 0]), np.array([0, 0, -1.])).reshape(2, 2, 2, 2)
        hx = bond_from_JB(np.diag([1., 1, 0]), np.zeros(3)).reshape(2, 2, 2, 2)
        one = np.eye(4).reshape(2, 2, 2, 2)
        T = np.einsum('abcd,efgh->aebfcgdh', hi, one) + np.einsum('abcd,efgh->aebfcgdh', one, hx)
        edge = np.kron(-PAULI[2] / 2, I2)
        return T.reshape(16, 16), edge, 4
    raise SystemExit('unknown model')


def renv(M):
    L = len(M); R = [None] * (L + 1); R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        A = M[j]; R[j] = np.einsum('asc,cd,bsd->ab', np.conj(A), R[j + 1], A)
    return R


def local(M, R, ops):
    """<O_j> for every site (left-canonical M: left environment is the identity)."""
    out = np.zeros((len(ops), len(M)), complex)
    for j, A in enumerate(M):
        for k, O in enumerate(ops):
            out[k, j] = np.einsum('asc,st,atd,cd->', np.conj(A), O, A, R[j + 1])
    return out


def corr(M, R, O, i0, rmax):
    """<O_i0 O_{i0+r}>, r = 1..rmax. The left environment of site i0 is the identity for left-canonical M."""
    A = M[i0]; T = np.einsum('asc,st,atd->cd', np.conj(A), O, A)
    out = []
    for r in range(1, rmax + 1):
        B = M[i0 + r]
        out.append(np.einsum('ab,asc,st,btd,cd->', T, np.conj(B), O, B, R[i0 + r + 1]))
        T = np.einsum('ab,asc,bsd->cd', T, np.conj(B), B)
    return np.array(out)


def main():
    name, L, CHI, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    h, edge, d = model(name)
    chis = tuple(c for c in (16, 32, 64, 96, 128, 160, 192, 256, 320) if c < CHI) + (CHI, CHI)
    t0 = time.time()
    r = D.dmrg(h, L, chis=chis, seed=0, edge_field=edge)
    S = r['S']; cp, cs, mid = profile_fits(S, L)
    res = dict(model=name, L=L, chis=chis, E=r['E'], trunc=r['hist'][-1][2], dE_last=r['hist'][-1][1] - r['hist'][-2][1],
               S=[float(x) for x in S], c_prof=cp, c_sm=cs, Smid=mid)
    if d == 2:
        M = r['M']; R = renv(M)
        nrm = np.einsum('asc,asd,cd->', np.conj(M[0]), M[0], R[1]).real
        loc = (local(M, R, PAULI) / nrm).real
        j = np.arange(L); w = (j >= L // 4) & (j < 3 * L // 4)
        res['stag'] = [float(np.mean(((-1.0) ** j * loc[a])[w])) for a in range(3)]
        res['unif'] = [float(np.mean(loc[a][w])) for a in range(3)]
        res['loc'] = loc.tolist()
        i0 = L // 4; rmax = L // 2
        res['corr'] = {}
        for a, nm in enumerate('xyz'):
            c = (corr(M, R, PAULI[a], i0, rmax) / nrm).real - loc[a][i0] * loc[a][i0 + 1:i0 + rmax + 1]
            res['corr'][nm] = c.tolist()
            k = np.linspace(0, np.pi, 401)
            win = np.hanning(2 * rmax + 1)[rmax + 1:]
            Sk = 1 + 2 * (np.cos(np.outer(k, np.arange(1, rmax + 1))) @ (c * win))
            res.setdefault('Sk_peak', {})[nm] = [float(k[np.argmax(Sk)]), float(Sk.max()), float(Sk[-1])]
    res['time_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), default=float)
    print(f"{name} L={L} chi={CHI} E={res['E']:.10f} trunc={res['trunc']:.1e} c_prof={cp:.3f} c_sm={cs:.3f} "
          + (f"stag={np.round(res['stag'],4).tolist()} Sk={res['Sk_peak']}" if d == 2 else '') + f" t={res['time_s']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
