#!/usr/bin/env python3
# GPT-6 Codex 2026-11-09: C183 fc.py wrapper; numerical solver unchanged.
# Additions: fixed source imports, physical X/Y and hopping profiles, local MPS checkpoint, corrected V/Delta metadata.
"""[C183 copy of C177 fc.py: adds --cp (full couplings as a JSON string, incl. Kxz) and --kxz; otherwise unchanged.]
C177 (Claude Opus 5.5 agent, 2026-11-09): open-chain DMRG of the full qubit chain at the Fibonacci cusp.

Uses C166's validated dmrg.py unchanged (two-site, complex, edge_field slot) and c166lib's observables.
Bond: band.bond(Jx, Jy, D, Jz, hz, Jxz, hx), i.e. h = Jx XX + Jy YY + D(XY - YX) + Jz ZZ + (hz/2)(Z1 + 1Z)
      + Jxz (XZ + ZX) + (hx/2)(X1 + 1X).
The open chain carries the full on-site field on every site (the two end sites get the missing half through
edge_field), so it is the plain open chain of the bond.  A constant e_shift = Jz + hz per bond (the classical
energy of the vacuum bond) is subtracted from h so that ARPACK's relative tolerance acts on the small
low-energy scale; it changes no state.  Energies reported are with the shift removed (E, ebond, emid).
Couplings come from tune.json (band-tuned D and h_z at V/Delta = r) unless given explicitly.

Observables: S(l), c_prof, c_sm, S_mid (dmrg_c.profile_fits), Sbar = mean of S over the three middle cuts,
<Z_j>, density n_j = (1 - <Z_j>)/2, bond energies and emid (mean over the central L/4 bonds), the hopping
profiles <s+_j s-_{j+1}> and <s+_j s-_{j+2}> (complex), and the row G(r) = <s+_{i0} s-_{i0+r}>, i0 = L/4,
whose Fourier weights at the two valley momenta measure the valley content.
usage: fc.py TAG L --r 8 --chis 16 32 64 ... [--hz X] [--D X] [--jxpy X] [--jxz X] [--hx X] [--outdir oc]
"""
import os, sys, json, time, argparse
from pathlib import Path
SOURCE=Path(__file__).resolve().parents[2]/"codex/work_C183"
sys.path.insert(0, str(SOURCE))
import numpy as np
import dmrg as D
import dmrg2 as D2
import dmrg_c as DC
import c166lib as X
import band as B

SP = np.array([[0, 1], [0, 0]], complex)   # s+ : |down> -> |up>, removes a particle (particle = down spin)
SM = SP.T.copy()                            # s- : creates a particle
SZ = B.SZ


def load_couplings(r, here):
    T = json.load(open(os.path.join(here, 'tune.json')))
    for d in T:
        if abs(d['r'] - r) < 1e-9:
            return dict(d['couplings'])
    raise KeyError(r)


def correlator_row(M, i, O1, O2):
    """<O1_i O2_j> for all j > i (complex), M left-canonical with the centre on the last site."""
    L = len(M)
    R = X._right_ids(M)
    Lf = np.ones((1, 1), complex)
    for k in range(i):
        A = M[k]
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    A = M[i]
    t = np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), O1, A, optimize=True)
    out = []
    for j in range(i + 1, L):
        A = M[j]
        Y = np.einsum('ab,asc,st,btd->cd', t, np.conj(A), O2, A, optimize=True)
        out.append(complex(np.einsum('cd,cd->', Y, R[j + 1], optimize=True)))
        t = np.einsum('ab,asc,bsd->cd', t, np.conj(A), A, optimize=True)
    return np.array(out)


def bond_corr_row(M, i, A2, B2):
    """<A2_{i,i+1} B2_{j,j+1}> for j = i+2..L-2, A2 and B2 given as (op_left, op_right) product pairs."""
    L = len(M)
    R = X._right_ids(M)
    Lf = np.ones((1, 1), complex)
    for k in range(i):
        A = M[k]
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    t = Lf
    for s, O in ((i, A2[0]), (i + 1, A2[1])):
        A = M[s]
        t = np.einsum('ab,asc,st,btd->cd', t, np.conj(A), O, A, optimize=True)
    out = []
    for j in range(i + 2, L - 1):
        y = t
        for s, O in ((j, B2[0]), (j + 1, B2[1])):
            A = M[s]
            y = np.einsum('ab,asc,st,btd->cd', y, np.conj(A), O, A, optimize=True)
        out.append(complex(np.einsum('cd,cd->', y, R[j + 2], optimize=True)))
        A = M[j]
        t = np.einsum('ab,asc,bsd->cd', t, np.conj(A), A, optimize=True)
    return np.array(out)


def pp_row(M, i, th):
    """<P_i P_j>, P_b = X_b X_{b+1} + Y_b Y_{b+1} = 2 (s+_b s-_{b+1} + h.c.) in the physical frame (valley-polarization
    order parameter, odd under the J_x + J_y = 0 symmetry); M is in the gauge frame with angles th."""
    L = len(M)
    def terms(b):
        ph = np.exp(-1j * (th[b] - th[b + 1]))
        return [(2 * ph, (SP, SM)), (2 * np.conj(ph), (SM, SP))]
    tot = 0
    for ca, A2 in terms(i):
        for j0 in [None]:
            pass
        rows = []
        for cb_idx in range(2):
            # phases of the second bond depend on j through th[j] - th[j+1]; build per j below
            Bop = (SP, SM) if cb_idx == 0 else (SM, SP)
            rows.append(bond_corr_row(M, i, A2, Bop))
        js = np.arange(i + 2, L - 1)
        ph = np.exp(-1j * (th[js] - th[js + 1]))
        tot = tot + ca * (2 * ph * rows[0] + 2 * np.conj(ph) * rows[1])
    return tot.real


def hop_profiles(M, rng):
    """<s+_j s-_{j+r}> for r in rng and all j (complex)."""
    L = len(M)
    R = X._right_ids(M)
    out = {r: np.zeros(L - r, complex) for r in rng}
    Lf = np.ones((1, 1), complex)
    for j in range(L - 1):
        A = M[j]
        t = np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), SP, A, optimize=True)
        for r in range(1, max(rng) + 1):
            if j + r >= L:
                break
            B_ = M[j + r]
            if r in rng:
                Y = np.einsum('ab,asc,st,btd->cd', t, np.conj(B_), SM, B_, optimize=True)
                out[r][j] = np.einsum('cd,cd->', Y, R[j + r + 1], optimize=True)
            t = np.einsum('ab,asc,bsd->cd', t, np.conj(B_), B_, optimize=True)
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tag'); ap.add_argument('L', type=int)
    ap.add_argument('--r', type=float, default=8.0)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 32, 64, 64])
    ap.add_argument('--hz', type=float, default=None)
    ap.add_argument('--mu', type=float, default=None, help='dressed one-particle mu_eff (tune.json linear map); overrides the tuned -2.2')
    ap.add_argument('--D', type=float, default=None)
    ap.add_argument('--jxpy', type=float, default=0.0)
    ap.add_argument('--jxz', type=float, default=0.0)
    ap.add_argument('--hx', type=float, default=0.0)
    ap.add_argument('--kxz', type=float, default=0.0)
    ap.add_argument('--cp', default=None, help='C183: JSON dict of all couplings (Jx Jy D Jz hz Jxz Kxz hx); overrides the rest')
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--res', type=float, default=1e-8)
    ap.add_argument('--old', action='store_true', help='use dmrg.py (ARPACK) instead of dmrg2.py')
    ap.add_argument('--cutoff', type=float, default=1e-12)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--noshift', action='store_true')
    ap.add_argument('--outdir', default='oc')
    ap.add_argument('--nocorr', action='store_true')
    ap.add_argument('--parity', default=None, choices=['even', 'odd'], help='select a fermion-parity sector with lam * prod Z, lam = 1')
    a = ap.parse_args()
    here = str(SOURCE)
    os.makedirs(a.outdir, exist_ok=True)
    cp = load_couplings(a.r, here) if a.cp is None else {}
    if a.mu is not None and a.cp is None:
        Tn = [d for d in json.load(open(os.path.join(here, 'tune.json'))) if abs(d['r'] - a.r) < 1e-9][0]
        cp['hz'] = cp['hz'] - (a.mu - Tn['target_mu']) * Tn['t2'] / Tn['dc0_dhz']
    if a.hz is not None: cp['hz'] = a.hz
    if a.D is not None: cp['D'] = a.D
    cp['Jx'] = (1.0 + a.jxpy) / 2; cp['Jy'] = (a.jxpy - 1.0) / 2
    cp['Jxz'] = a.jxz; cp['hx'] = a.hx; cp['Kxz'] = a.kxz
    if a.cp is not None:
        cp = {k: float(v) for k, v in (json.load(open(a.cp[1:])) if a.cp.startswith('@') else json.loads(a.cp)).items()}
    h = B.bond(**cp)
    shift = 0.0 if a.noshift else cp['Jz'] + cp['hz']
    hs = h - shift * np.eye(4)
    L = a.L
    edge = 0.5 * (cp['hz'] * SZ + cp.get('hx', 0.0) * B.SX)
    plam = None if a.parity is None else (1.0 if a.parity == 'odd' else -1.0)
    fn = os.path.join(a.outdir, f"{a.tag}_{L}.json")
    t0 = time.time()
    # staggered real gauge exp(-i theta_j Z_j/2), theta_j = (-1)^j pi/4, used when it makes every bond real
    th = np.array([(np.pi / 4) * (-1) ** j for j in range(L)])
    u = lambda t: np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])
    hg = [np.kron(u(th[j]), u(th[j + 1])) @ hs @ np.kron(u(th[j]), u(th[j + 1])).conj().T for j in range(2)]
    gauge = (not a.old) and max(np.abs(x.imag).max() for x in hg) < 1e-14
    if not gauge:
        th = np.zeros(L)
    if a.old:
        r = D.dmrg(hs, L, chis=tuple(a.chis), cutoff=a.cutoff, seed=a.seed, verbose=True, lanczos_tol=a.tol, edge_field=edge)
        hl = [h] * (L - 1)
    elif gauge:
        hgl = [hg[j % 2].real.copy() for j in range(L - 1)]
        r = D2.dmrg(hgl, L, chis=tuple(a.chis), cutoff=a.cutoff, seed=a.seed, verbose=True, res_tol=a.res, edge_field=edge.real.copy(), parity_lam=plam)
        hl = [x + shift * np.eye(4) for x in hgl]
    else:
        r = D2.dmrg(hs, L, chis=tuple(a.chis), cutoff=a.cutoff, seed=a.seed, verbose=True, res_tol=a.res, edge_field=edge, parity_lam=plam)
        hl = [h] * (L - 1)
    S = np.asarray(r['S'])
    cpf, csm, mid = DC.profile_fits(S, L)
    c = L // 2
    Sbar = float(np.mean(S[c - 2:c + 1]))
    M = r['M']
    mag = X.onesite(M, [B.SX, B.SY, SZ])
    Xprof,Yprof,Z = mag
    # Undo any on-site Z rotation used for a real solver gauge.
    transverse=(Xprof+1j*Yprof)*np.exp(-1j*th)
    Xprof=transverse.real;Yprof=transverse.imag
    eb = [X.bond_energies(M, hl[p]) for p in range(min(2, L - 1))]
    e = np.array([eb[j % 2][j] for j in range(L - 1)])
    w = L // 8
    emid = float(np.mean(e[c - w:c + w]))
    Pval=None
    if plam is not None:
        pe=np.ones((1,1),complex)
        for A in M:
            pe=np.einsum('ab,asc,st,btd->cd',pe,np.conj(A),SZ,A,optimize=True)
        Pval=float(np.real(pe[0,0]))
    pshift = 0.0 if plam is None else -plam*Pval
    out = dict(tag=a.tag, L=L, r=(4*cp['Jz']/abs(cp['Jx']-cp['Jy']) if cp['Jx']!=cp['Jy'] else None), couplings=cp, shift=shift, chis=a.chis, E=r['E'] + shift * (L - 1) + pshift, parity=a.parity,
               P_total=Pval,
               hist=[(x[0], x[1] + shift * (L - 1), x[2]) for x in r['hist']], trunc=r['hist'][-1][2],
               c_prof=cpf, c_sm=csm, Smid=mid, Sbar=Sbar, S=S.tolist(), X=Xprof.tolist(), Y=Yprof.tolist(), Z=Z.tolist(),
               Xmid=float(np.mean(Xprof[c-w:c+w])),Ymid=float(np.mean(Yprof[c-w:c+w])),
               n_mid=float(np.mean((1 - Z[c - w:c + w]) / 2)), N=float(np.sum((1 - Z) / 2)),
               ebond=e.tolist(), emid=emid, t_dmrg=time.time() - t0, gauge=bool(gauge), stats=r.get('stats'))
    # Save a local, git-ignored checkpoint and a cheap direct hopping-response profile.
    checkpoint_dir=os.environ.get("CD_DMRG_CHECKPOINT_DIR","/private/tmp/cd-conjecture-mps")
    os.makedirs(checkpoint_dir,exist_ok=True)
    statefile=os.path.join(checkpoint_dir,f"{a.tag}_{L}_mps.npz")
    np.savez_compressed(statefile,**{f"M{j}":A for j,A in enumerate(M)},gauge_angles=th)
    Pbond=np.kron(B.SX,B.SX)+np.kron(B.SY,B.SY)
    Pop=[]
    for j in range(2 if gauge else 1):
        ug=np.kron(u(th[j]),u(th[j+1]));Pop.append(ug@Pbond@ug.conj().T)
    Pvalues=[X.bond_energies(M,p) for p in Pop]
    Pprofile=np.array([Pvalues[j%len(Pop)][j] for j in range(L-1)])
    out['P_hop']=Pprofile.tolist();out['P_hop_mid']=float(np.mean(Pprofile[c-w:c+w]))
    out['mps_checkpoint']=statefile
    if not a.nocorr:
        hp = hop_profiles(M, (1, 2))
        for rr in (1, 2):  # undo the gauge: physical <s+_i s-_j> = exp(-i th_i + i th_j) <...>'
            hp[rr] = hp[rr] * np.exp(-1j * th[:L - rr] + 1j * th[rr:])
        out['hop1'] = [[z.real, z.imag] for z in hp[1]]
        out['hop2'] = [[z.real, z.imag] for z in hp[2]]
        i0 = L // 4
        G = correlator_row(M, i0, SP, SM) * np.exp(-1j * th[i0] + 1j * th[i0 + 1:])
        out['G_i0'] = i0; out['G'] = [[z.real, z.imag] for z in G]
        out['P'] = (4 * hp[1].real).tolist()          # <X_j X_{j+1} + Y_j Y_{j+1}> profile
        out['PP'] = pp_row(M, i0, th).tolist()         # <P_{i0} P_j>, j = i0 + 2 .. L - 2
    out['time_s'] = time.time() - t0
    json.dump(out, open(fn + '.tmp', 'w'), default=float); os.replace(fn + '.tmp', fn)
    print(f"done {fn} E={out['E']:.10f} Sbar={Sbar:.4f} c_prof={cpf:.3f} trunc={out['trunc']:.1e} n_mid={out['n_mid']:.4f} "
          f"N={out['N']:.3f} emid={emid:.10f} t={out['time_s']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
