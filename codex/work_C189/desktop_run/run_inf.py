"""C189 production driver for idmrg.py. Usage: python3 run_inf.py <cfg.json> <outdir>
cfg keys: tag, p, q, wa, phi, [t2, V2, V3, blockade, fermionic], chis (list), min_it, max_it, tolS, tolE, chi0.
gauge is 1 when |phi - pi/2| < 0.2 (C176's rule), else 0.
At every chi the window is iterated until |S_t - S_{t-2}| < tolS and |e_t - e_{t-2}| < tolE (two-iteration
comparison, because a 2q-periodic state alternates between the two mid bonds), with at least min_it and at most
max_it iterations. Then the transfer spectra are measured and one jsonl record is written. Iteration lines go to
<outdir>/<tag>.iter.jsonl, chi records to <outdir>/<tag>.jsonl, a checkpoint to <outdir>/ckpt/<tag>.pkl.
"""
import numpy as np, json, sys, os, time, pickle
from idmrg import Model, IDMRG


def make_model(c):
    phi = c.get('phi', 0.0)
    w = c['wa'] * np.exp(1j * phi) if abs(phi) > 1e-12 else c['wa']
    gauge = c.get('gauge', 1 if abs(phi - np.pi / 2) < 0.2 else 0)
    return Model(c['p'], c['q'], w, t2=c.get('t2', -1.0), V2=c.get('V2', -2.0), V3=c.get('V3', -2.0 / 3.0),
                 blockade=c.get('blockade', True), gauge=gauge, fermionic=c.get('fermionic', False))


def spec_rec(ev, ev0, n2):
    ev = np.asarray(ev)
    lam0 = abs(ev0[0])
    out = []
    for x in ev:
        r = abs(x) / lam0
        xi = float(-n2 / np.log(r)) if 0 < r < 1 else (float('inf') if r >= 1 else 0.0)
        out.append([float(np.real(x)), float(np.imag(x)), xi, float(np.angle(x) / n2)])
    return out


def main():
    c = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    os.makedirs(os.path.join(outdir, 'ckpt'), exist_ok=True)
    tag = c['tag']
    M = make_model(c)
    ck = os.path.join(outdir, 'ckpt', f'{tag}.pkl')
    done = set()
    fout = os.path.join(outdir, f'{tag}.jsonl')
    if os.path.exists(fout):
        for l in open(fout):
            if l.strip():
                done.add(json.loads(l)['chi'])
    if os.path.exists(ck):
        with open(ck, 'rb') as f:
            dm = pickle.load(f)
        dm.M = M; dm.W.M = M
        if not hasattr(dm.W, 'noise'):
            dm.W.noise = c.get('v0noise', 0.0)
    else:
        dm = IDMRG(M, chi0=c.get('chi0', 16))
        dm.W.noise = c.get('v0noise', 0.0)
        dm.sweep_window(2 * c.get('chi0', 16), nsw=4)
    n2 = 2 * M.q
    for chi in c['chis']:
        if chi in done:
            continue
        hist = []
        t0 = time.time()
        for it in range(c.get('max_it', 400)):
            t1 = time.time()
            r = dm.iterate(chi, tol=c.get('tol_eig', 1e-10))
            r['secs'] = time.time() - t1
            hist.append(r)
            with open(os.path.join(outdir, f'{tag}.iter.jsonl'), 'a') as f:
                f.write(json.dumps(r) + '\n')
            if it + 1 >= c.get('min_it', 20) and len(hist) >= 3:
                dS = abs(hist[-1]['S'] - hist[-3]['S']); dE = abs(hist[-1]['e_site'] - hist[-3]['e_site'])
                if dS < c.get('tolS', 1e-5) and dE < c.get('tolE', 1e-9):
                    break
            if it % 10 == 0:
                with open(ck, 'wb') as f:
                    pickle.dump(dm, f)
        conv = dict(iters=len(hist), dS=abs(hist[-1]['S'] - hist[-3]['S']) if len(hist) >= 3 else None,
                    dE=abs(hist[-1]['e_site'] - hist[-3]['e_site']) if len(hist) >= 3 else None)
        t2 = time.time()
        meas = {}
        try:
            ev0 = dm.transfer_eigs(0, k=10)
            meas['T0'] = spec_rec(ev0, ev0, n2)
            for nm, dq, st in (('T1b', 1, False), ('T1f', 1, True), ('T2', 2, False), ('T0f', 0, True)):
                meas[nm] = spec_rec(dm.transfer_eigs(dq, string=st, k=6), ev0, n2)
        except Exception as e:
            meas['error'] = repr(e)
        try:
            meas['dens'] = dm.window_density()
        except Exception as e:
            meas['dens_error'] = repr(e)
        spec = dm.last['spec']
        rec = dict(tag=tag, cfg=c, chi=int(chi), e_site=hist[-1]['e_site'], e_prev=hist[-2]['e_site'],
                   S=hist[-1]['S'], S_prev=hist[-2]['S'], trunc=hist[-1]['trunc'], chi_mid=hist[-1]['chi_mid'],
                   conv=conv, secs=time.time() - t0, meas_secs=time.time() - t2, meas=meas,
                   schmidt=[float(x) for x in spec[:64]], nsv=len(spec))
        with open(fout, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        with open(ck, 'wb') as f:
            pickle.dump(dm, f)
        x = meas.get('T0', [[0, 0, 0]] * 2)
        print(f"{tag} chi={chi} it={len(hist)} e={rec['e_site']:.10f} S={rec['S']:.5f}/{rec['S_prev']:.5f} "
              f"xi0={x[1][2] if len(x) > 1 else None} trunc={rec['trunc']:.1e} t={rec['secs']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
