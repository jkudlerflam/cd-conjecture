#!/usr/bin/env python3
"""Full microscopic OBC screen using the unchanged, validated C190 Z2 solver.

No projected Hilbert space. The optional staggered rotation is computational;
all couplings and reported observables refer to the uniform original chain.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / 'codex/work_C190'))
import z2dmrg as Z


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n')
    os.replace(tmp, path)


def build_model(L, p, distance=1):
    h2 = Z.pauli_bond(Jx=p.get('Jx', 0), Jy=p.get('Jy', 0),
                      Jz=p.get('Jz', 0), D=p.get('Dz', 0))
    # Field is physical on every site, including both ends.
    hz = p.get('hz', 0)
    shift = -(p.get('Jz', 0) + hz)
    h1 = hz * Z.ZZ_ + shift * Z.I2
    gauge = 'stag' if distance == 1 and abs(p.get('Jx', 0) + p.get('Jy', 0)) < 1e-14 else None
    gph = Z.gauge_phases(L, gauge)
    if distance == 2:
        assert p.get('Jz', 0) == p.get('Dz', 0) == 0 and p['Jx'] == p['Jy']
        mpo = Z.MPO(L, [np.zeros((4, 4))] * (L - 1), [h1] * L,
                    nnn=[(2*p['Jx'], 2*p['Jy'])] * (L - 2))
    else:
        h2s = [Z.rotate_bond(h2, gph[j], gph[j+1]) for j in range(L-1)] if gauge else [h2]*(L-1)
        mpo = Z.MPO(L, h2s, [h1]*L)
    return mpo, gph, shift


def entropy_stats(S, L):
    S = np.asarray(S)
    x = np.arange(1, L)
    central3 = float(np.mean(S[(x >= L//2-1) & (x <= L//2+1)]))
    lo, hi = 3*L/8, 5*L/8
    xx = np.r_[lo, x[(x > lo) & (x < hi)], hi]
    fractional = float(np.trapezoid(np.interp(xx, x, S), xx)/(hi-lo))
    profiles = {}
    for fraction in (0.25, 0.375):
        mask = (x >= fraction*L) & (x <= (1-fraction)*L)
        chord = np.log((2*L/np.pi)*np.sin(np.pi*x[mask]/L))
        coef = np.polyfit(chord, S[mask], 1)
        profiles[str(fraction)] = dict(c=6*float(coef[0]),
            rms=float(np.sqrt(np.mean((np.polyval(coef, chord)-S[mask])**2))))
    return dict(S_central3=central3, S_fractional=fractional, profiles=profiles)


def long_correlations(psi, gph, density):
    """Physical XX,YY,ZZ from L/4 through 3L/4; parity makes <X>=<Y>=0."""
    state=psi.copy()
    origin=psi.L//4
    for j in range(origin):
        Z._qr_move_right(state,j)
    ops={'X':[g*Z.SP+g.conjugate()*Z.SM for g in gph],
         'Y':[-1j*g*Z.SP+1j*g.conjugate()*Z.SM for g in gph],
         'Z':[Z.ZZ_]*psi.L}
    out={'origin':origin}
    for axis, oo in ops.items():
        C=state.M[origin]
        env=Z._apply_site(np.eye(C.shape[0]),C,oo[origin])
        values=[]
        for j in range(origin+1, min(psi.L,3*psi.L//4+1)):
            value=float(np.trace(Z._apply_site(env,state.M[j],oo[j])).real)
            if axis=='Z':
                value-=(1-2*density[origin])*(1-2*density[j])
            values.append([j-origin,value])
            env=Z._apply_site(env,state.M[j],None)
        out[axis]=values
    return out


def run_point(point, L, parity, chi=64, seed=17, resume=False, stages=None, max_sweeps=6,
              output_dir=None, checkpoint_dir=None):
    outdir = Path(output_dir or HERE/'dmrg')
    ckdir = Path(checkpoint_dir or '/private/tmp/d2-order-one-checkpoints')
    ckdir.mkdir(parents=True, exist_ok=True)
    label = f"{point['name']}_L{L}_P{parity}_s{seed}"
    out = outdir / (label + '.json')
    ckpt = ckdir / (label + '_mps.npz')
    p = point['couplings']
    distance = point.get('distance', 1)
    mpo, gph, shift = build_model(L, p, distance)
    history = []
    if resume:
        old = json.loads(out.read_text())
        assert old['couplings'] == p and old['L'] == L and old['parity'] == parity and old.get('distance',1) == distance
        with np.load(ckpt,allow_pickle=False) as saved:
            saved_meta=json.loads(str(saved['meta']))
        if 'history' in saved_meta:
            assert len(saved_meta['history'])==len(old['history']), 'Checkpoint/JSON sweep mismatch; recover from checkpoint metadata explicitly'
        psi = Z.MPS.load(ckpt)
        history = old['history']
    else:
        assert not out.exists(), f'Output exists; use --resume: {out}'
        psi = Z.random_mps(L, parity, D0=4, seed=seed, dtype=mpo.dtype)
    record = dict(name=point['name'], L=L, parity=parity, seed=seed, couplings=p, distance=distance,
        solver='codex/work_C190/z2dmrg.py',
        solver_sha256=hashlib.sha256((ROOT/'codex/work_C190/z2dmrg.py').read_bytes()).hexdigest(),
        checkpoint=str(ckpt), status='running', history=history)
    def checkpoint(psi, _):
        temp = ckpt.with_name(ckpt.stem + '.tmp.npz')
        psi.save(temp, meta=record)
        os.replace(temp, ckpt)
        atomic_json(out, record)
    def log(row):
        row['sweep_global'] = len(history)
        row['E_physical'] = row['E'] + mpo.const - L*shift
        row.update(entropy_stats(row['Sprof'], L) if 'Sprof' in row else {})
        history.append(row)
        print(json.dumps({k:row[k] for k in ['sweep_global','chi','E_physical','res_max','trunc','t','S_fractional'] if k in row}), flush=True)
    if stages is None:
        stages = [chi] if resume else sorted(set([min(16,chi), min(32,chi), chi]))
    for stage_idx, D in enumerate(stages):
        initial = not resume and stage_idx == 0
        sweeps = []
        if initial:
            sweeps += [dict(chi=D, alpha=1e-3, res_tol=1e-8, maxit=500),
                       dict(chi=D, alpha=1e-5, res_tol=1e-9, maxit=500)]
        sweeps += [dict(chi=D, alpha=0., res_tol=1e-10, maxit=600)] * max_sweeps
        def stop(h):
            clean = [r for r in h if r['alpha'] == 0]
            return len(clean) >= 3 and abs(clean[-1]['E']-clean[-2]['E']) < 1e-8 and \
                abs(clean[-1]['S_fractional']-clean[-2]['S_fractional']) < 1e-3 and \
                clean[-1]['res_max'] <= 1e-8
        psi, _ = Z.dmrg(psi, mpo, sweeps, cutoff=1e-14, log=log, ckpt=checkpoint, stop=stop)
    assert psi.check()
    S = Z.entropies(psi)
    m = Z.measure(psi, mpo, gph=gph, gwin=min(12,L), cp_rmax=min(12,L//2))
    # C190 leaves this placeholder zero unless explicit energy operators are supplied.
    m.pop('e_bond', None)
    density = np.asarray(m['dens'])
    m['long_connected_correlations']=long_correlations(psi,gph,density)
    bulk = density[L//4:3*L//4]
    record.update(status='done', chi=chi, E_physical=Z.energy(psi, mpo)+mpo.const-L*shift,
        Sprof=S.tolist(), **entropy_stats(S,L), measurements=m,
        density_bulk=float(np.mean(bulk)), density_bulk_range=float(np.ptp(bulk)),
        parity_expectation=Z.parity_expectation(psi))
    clean = [r for r in history if r['alpha'] == 0 and r['chi'] == chi]
    record['convergence'] = dict(local_residual=clean[-1]['res_max'], discarded_weight=clean[-1]['trunc'],
        delta_E=abs(clean[-1]['E']-clean[-2]['E']) if len(clean)>1 else None,
        delta_S=abs(clean[-1]['S_fractional']-clean[-2]['S_fractional']) if len(clean)>1 else None,
        canonical_vs_sweep_S=record['S_fractional']-clean[-1]['S_fractional'])
    cc=record['convergence']
    cc['sweep_converged']=bool(cc['local_residual']<=1e-8 and cc['delta_E'] is not None and cc['delta_E']<=1e-8
                               and cc['delta_S'] is not None and cc['delta_S']<=1e-3)
    cc['bond_dimension_converged']='not established by a completed schedule alone'
    checkpoint(psi, [])
    return record


def validate():
    from ed_screen import build_sector
    from scipy.sparse.linalg import eigsh
    configs = [dict(name='validation_symmetric',couplings=dict(Jx=.5,Jy=-.5,Jz=.5,Dz=.1875,hz=-.5)),
               dict(name='validation_attractive',couplings=dict(Jx=.5,Jy=-.25,Jz=-.625,Dz=.5,hz=0)),
               dict(name='validation_xx',couplings=dict(Jx=.5,Jy=.5,Jz=0,Dz=0,hz=0)),
               dict(name='validation_xx2',distance=2,couplings=dict(Jx=.5,Jy=.5,Jz=0,Dz=0,hz=0))]
    rows=[]
    for point in configs:
        for P in (0,1):
            H,basis=build_sector(8,point['couplings'],parity=P,periodic=False,
                                 distance=point.get('distance',1))
            ev,v=eigsh(H,k=2,which='SA',tol=1e-13,
                       v0=np.random.default_rng(17).normal(size=H.shape[0]))
            order=np.argsort(ev);ev=ev[order];v=v[:,order]
            full=np.zeros(256,complex);full[basis]=v[:,0]
            p=np.linalg.svd(full.reshape(16,16),compute_uv=False)**2
            exactS=Z.entropy(p)
            exists=(HERE/'validation_runs'/f"{point['name']}_L8_P{P}_s17.json").exists()
            result=run_point(point,8,P,chi=16,stages=[16],max_sweeps=5,resume=exists,
                             output_dir=HERE/'validation_runs')
            row=dict(name=point['name'],parity=P,energy_error=abs(result['E_physical']-ev[0]),
                sector_gap=float(ev[1]-ev[0]),entropy_comparison_valid=bool(ev[1]-ev[0]>1e-9),
                entropy_error=abs(result['Sprof'][3]-exactS) if ev[1]-ev[0]>1e-9 else None)
            assert row['energy_error']<1e-8 and (row['entropy_error'] is None or row['entropy_error']<1e-7), row
            rows.append(row)
    atomic_json(HERE/'dmrg_validation.json',rows)
    print(json.dumps(rows,indent=2))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--validate',action='store_true')
    ap.add_argument('--config',type=Path)
    ap.add_argument('--name')
    ap.add_argument('--L',type=int,default=24)
    ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--chi',type=int,default=64)
    ap.add_argument('--seed',type=int,default=17)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--max-sweeps',type=int,default=6)
    args=ap.parse_args()
    if args.validate:
        validate();return
    cfg=json.loads(args.config.read_text())
    point=next(p for p in cfg['points'] if p['name']==args.name)
    run_point(point,args.L,args.parity,args.chi,args.seed,args.resume,max_sweeps=args.max_sweeps)


if __name__=='__main__':
    main()
