"""C194 finite-chain driver for the pair-breaking gap (C176's fibdmrg.py, unchanged).
Usage: python3 -E -s run194.py <tasks.json> <outdir> <nworkers>
tasks.json: list of dicts {tag, L, N, chis, [meas], [seed], [wa, phi, t2, V2, V3, blockade]}; defaults are the
symmetric line |w| = 1, phi = pi/2, bare couplings (C173 PROOF eq. 3.2), blockade on.
One jsonl line per sweep in <outdir>/<tag>.jsonl (E, trunc, S profile, maxchi, secs). After the last sweep a
final record with kind='final' holds the density profile (always) and, if meas is set, C176's full measurement
(analyze.reduce_meas with the one-body matrices stored). Tasks whose final record exists are skipped. Workers pull
tasks in the order given (put the expensive ones first). Each worker is single threaded."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import numpy as np, json, sys, time, multiprocessing as mp
from fibdmrg import Chain, DMRG, ID, NN
import analyze


def density(dm):
    ch = dm.ch; L = ch.L
    Tl = dm._site_tensors()
    E0 = {(('L', s), ('L', t)): np.ones((1, 1)) for s in (0, 1) for t in (0, 1)}
    Eid = [E0]
    for i in range(L - 1):
        Eid.append(dm._transfer(Tl, Eid[-1], i, ID))
    norm = dm._close(dm._transfer(Tl, Eid[L - 1], L - 1, ID))
    return [float(np.real(dm._close(dm._transfer(Tl, Eid[i], i, NN))) / np.real(norm)) for i in range(L)]


def run_task(c, outdir):
    tag = c['tag']; outf = os.path.join(outdir, tag + '.jsonl')
    if os.path.exists(outf):
        for l in open(outf):
            if '"kind": "final"' in l:
                return tag, 'skip'
        os.remove(outf)
    phi = c.get('phi', np.pi / 2)
    w = c.get('wa', 1.0) * np.exp(1j * phi) if abs(phi) > 1e-12 else c.get('wa', 1.0)
    gauge = c.get('gauge', 1 if abs(phi - np.pi / 2) < 0.2 else 0)
    kw = {k: c[k] for k in ('t2', 'V2', 'V3', 'blockade') if k in c}
    ch = Chain(c['L'], c['N'], w, gauge=gauge, **kw)
    dm = DMRG(ch, seed=c.get('seed', 1))
    t0 = time.time()
    for k, chi in enumerate(c['chis']):
        t1 = time.time()
        E, tr = dm.sweep(chi, tol=c.get('tol', 1e-11))
        rec = dict(tag=tag, kind='sweep', sweep=k, chi=int(chi), E=E, trunc=tr, maxchi=int(dm.chi_bond.max()),
                   Smid=float(dm.S[c['L'] // 2 - 1]), secs=time.time() - t1)
        with open(outf, 'a') as f:
            f.write(json.dumps(rec) + '\n')
    fin = dict(tag=tag, kind='final', cfg=c, E=E, S=dm.S.tolist(), trunc_bond=dm.trunc.tolist(), secs=time.time() - t0)
    fin['dens'] = density(dm)
    if c.get('meas'):
        t1 = time.time()
        m = dm.measure()
        fin['meas'] = analyze.reduce_meas(m, c['L'], c['N'], store_G=True)
        fin['meas_secs'] = time.time() - t1
    with open(outf, 'a') as f:
        f.write(json.dumps(fin) + '\n')
    return tag, f"E={E:.12f} t={time.time() - t0:.0f}s"


def worker(args):
    c, outdir = args
    try:
        tag, msg = run_task(c, outdir)
    except Exception as e:
        tag, msg = c['tag'], 'ERROR ' + repr(e)
    print(f"{time.strftime('%H:%M:%S')} {tag} {msg}", flush=True)
    return tag, msg


if __name__ == '__main__':
    tasks = json.load(open(sys.argv[1])); outdir = sys.argv[2]; nw = int(sys.argv[3])
    os.makedirs(outdir, exist_ok=True)
    with mp.Pool(nw, maxtasksperchild=1) as pool:
        for _ in pool.imap_unordered(worker, [(c, outdir) for c in tasks], chunksize=1):
            pass
