"""Independent literal-qubit ED; real gauge is used only for diagonalization.

Bit j denotes site j, bit=1 means Z=-1. All observable amplitudes below
are in the ORIGINAL complex Pauli frame. No Jordan-Wigner boundary shortcut.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parents[1]


class Chain:
    def __init__(self, L, parity, periodic=True):
        self.L, self.parity, self.periodic = L, parity, periodic
        self.states = np.array([s for s in range(1 << L)
                                if s.bit_count() % 2 == parity], dtype=np.int64)
        self.n = len(self.states)
        self.inv = np.full(1 << L, -1, dtype=np.int64)
        self.inv[self.states] = np.arange(self.n)
        self.bits = ((self.states[None, :] >> np.arange(L)[:, None]) & 1)
        self.z = 1 - 2 * self.bits
        self.mag = self.z.sum(axis=0)
        self.phase = (1j) ** ((np.arange(L) @ self.bits) % 4)
        self.bonds = [(j, (j+1) % L) for j in range(L if periodic else L-1)]
        row, col, val = [], [], []
        diag = np.zeros(self.n)
        for j, k in self.bonds:
            dest = self.inv[self.states ^ (1 << j) ^ (1 << k)]
            amp = 1j * (9/16 * self.z[k] + 7/16 * self.z[j])
            amp = amp * self.phase / self.phase[dest]
            assert np.max(abs(amp.imag)) < 1e-13
            row.append(dest); col.append(np.arange(self.n)); val.append(amp.real)
            diag += self.z[j] * self.z[k]
        self.H0 = coo_matrix((np.concatenate(val),
                             (np.concatenate(row), np.concatenate(col))),
                            shape=(self.n, self.n)).tocsr() + diags(diag)
        self.trans_dest = self.inv[((self.states << 1) & ((1 << L)-1))
                                   | (self.states >> (L-1))]

    def H(self, h):
        return self.H0 - diags(h * self.mag)

    def translate(self, v):
        w = np.empty_like(v); w[self.trans_dest] = v
        return w

    def momentum_states(self, V):
        """Project returned eigenvectors, including incomplete degenerate blocks.

        This can recover translation partners omitted by a real eigensolver.
        Reported low levels are still not claimed to exhaust all degeneracies.
        """
        found = {}
        for v in V.T:
            vs = []; w = v.copy()
            for _ in range(self.L):
                vs.append(w); w = self.translate(w)
            projections = np.fft.fft(np.array(vs), axis=0) / self.L
            for k, w in enumerate(projections):
                for old in found.get(k, []):
                    w = w - old * np.vdot(old, w)
                norm = np.linalg.norm(w)
                if norm > 1e-6:
                    found.setdefault(k, []).append(w / norm)
        return [(k, v) for k in sorted(found) for v in found[k]]

    def local(self, name, j, v):
        if name == 'Z':
            return self.z[j] * v
        k = (j + (1 if name == 'bond' else 2)) % self.L
        dest = self.inv[self.states ^ (1 << j) ^ (1 << k)]
        if name == 'bond':
            amp = 1j * (9/16 * self.z[k] + 7/16 * self.z[j])
        elif name == 'kappa':
            amp = 1j * (self.z[k] - self.z[j])
        elif name == 'B':
            amp = 1 + self.z[j] * self.z[k]
        elif name == 'hop2':
            amp = 1 - self.z[j] * self.z[k]
        else:
            raise ValueError(name)
        w = np.zeros_like(v, dtype=complex); w[dest] = amp * v
        if name == 'bond':
            w += self.z[j] * self.z[k] * v
        return w

    def structure(self, name, v):
        L = self.L
        W = np.array([self.local(name, j, v) for j in range(L)])
        means = W @ v.conj()
        raw = W.conj() @ W.T
        conn = raw - np.outer(means.conj(), means)
        C = np.array([np.mean([conn[j, (j+r) % L] for j in range(L)])
                      for r in range(L)])
        sf = np.fft.fft(C).real
        assert np.min(sf) > -1e-8
        peak = int(np.argmax(sf))
        op = np.sum((-1.) ** np.arange(L)[:, None] * W, axis=0)
        return dict(mean=float(means.real.mean()),
                    mean_nonuniformity=float(np.ptp(means.real)),
                    connected_re=C.real.tolist(), connected_im=C.imag.tolist(),
                    structure=sf.tolist(), peak_k=peak,
                    peak_q_over_pi=2*peak/L, peak_over_L=float(sf[peak]/L),
                    pi_over_L=float(sf[L//2]/L),
                    half_chain_staggered=float((-1)**(L//2)*C[L//2].real)), op


def run(L, h, nev=6):
    start = time.monotonic(); sectors = []; allstates = []
    for p in (0, 1):
        c = Chain(L, p); H = c.H(h)
        rng = np.random.default_rng(811 + L + p)
        vals, V = eigsh(H, k=nev, which='SA', tol=2e-12,
                        maxiter=12000, ncv=max(32, 4*nev),
                        v0=rng.normal(size=c.n))
        order = np.argsort(vals); vals, V = vals[order], V[:, order]
        residual = max(np.linalg.norm(H @ V[:, i] - vals[i]*V[:, i])
                       for i in range(nev))
        assert residual < 1e-8, residual
        states = []
        for E in vals:
            if states and abs(E - states[-1]['E']) < 1e-8:
                continue
            mask = abs(vals-E) < 1e-8
            for k, v in c.momentum_states(c.phase[:, None] * V[:, mask]):
                vr = v / c.phase
                res = np.linalg.norm(H @ vr-E*vr)
                if res > 1e-7:
                    raise AssertionError((L,h,p,k,res))
                states.append(dict(E=float(E), k=int(k), v=v))
        sectors.append(dict(parity=p, energies=vals.tolist(), residual=float(residual),
                            levels=[dict(E=s['E'], k=s['k']) for s in states]))
        allstates.append((c, states))
    p = int(sectors[1]['energies'][0] < sectors[0]['energies'][0] - 1e-9)
    c, states = allstates[p]; E = min(x['E'] for x in states)
    gs = [s for s in states if s['E'] < E + 1e-8]
    s = min(gs, key=lambda s: (min(s['k'], L-s['k']), s['k']))
    v, kg = s['v'], s['k']
    obs = {}; channels = {}
    for name in ('Z', 'bond', 'B', 'kappa', 'hop2'):
        obs[name], op = c.structure(name, v)
        channels[name] = [dict(gap=t['E']-E, k=t['k'],
                               weight=float(abs(np.vdot(t['v'], op))**2/L))
                          for t in states]
    levels = sorted((s['E']-E,p0,s['k']) for p0,(_,ss) in enumerate(allstates) for s in ss)
    full = np.zeros(1 << L, dtype=complex); full[c.states] = v
    ent = []
    for cut in range(1,L//2+1):
        sv = np.linalg.svd(full.reshape(1 << cut, -1), compute_uv=False)**2
        sv = sv[sv>1e-15]; ent.append(float(-sv @ np.log(sv)))
    return dict(L=L, h=h, E=E, e=E/L, parity=p, k=kg,
                ground_momenta_returned=[s['k'] for s in gs], sectors=sectors,
                low_levels=[dict(gap=x, parity=p0, k=k) for x,p0,k in levels],
                magnetization=float(c.mag @ abs(v)**2 / L),
                opposite_parity_gap=sectors[1-p]['energies'][0]-E,
                observables=obs, staggered_spectral_weights=channels,
                half_chain_entropy=ent[-1], entropies=ent,
                runtime_seconds=time.monotonic()-start)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sizes', default='12,16')
    ap.add_argument('--fields', default='0,.5,1,1.4,1.6,1.7,1.738,1.7466,1.755,1.8,1.9,2,2.2,2.5,3,4,6')
    ap.add_argument('--seconds', type=int, default=360)
    ap.add_argument('--nev', type=int, default=6)
    ap.add_argument('--output', default='coarse.json')
    args = ap.parse_args()
    output = HERE / 'results' / args.output
    data = dict(parameters=vars(args), records=[], complete=False,
                conventions='Pauli; full field; periodic even ring; one worker/thread')
    def save():
        output.write_text(json.dumps(data, indent=2)+'\n')
    def timeout(*_):
        raise TimeoutError('registered wall-clock cap reached')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(args.seconds)
    start = time.monotonic()
    try:
        with threadpool_limits(limits=1):
            for L in map(int,args.sizes.split(',')):
                assert 6 <= L <= 20 and L % 2 == 0
                for h in map(float,args.fields.split(',')):
                    record=run(L,h,args.nev); data['records'].append(record); save()
                    print(f'L={L} h={h:g} e={record["e"]:.10f} M={record["magnetization"]:.6f} '
                          f'gapP={record["opposite_parity_gap"]:.6g} '
                          f'N2={record["observables"]["Z"]["pi_over_L"]:.5f} '
                          f'k2={record["observables"]["kappa"]["pi_over_L"]:.5f} '
                          f't={record["runtime_seconds"]:.2f}s', flush=True)
        data['complete']=True
    finally:
        signal.alarm(0); data['runtime_seconds']=time.monotonic()-start; save()


if __name__ == '__main__':
    main()
