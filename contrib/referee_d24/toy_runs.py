"""Driver for the toy vacancy-completion runs (see toy_ed.py).

Usage
  python3 -B toy_runs.py sectors L1 L2 ...          # b=2, N in [L/2-2, L/2+2]
  python3 -B toy_runs.py sectors_b b lo hi L1 ...   # N in [L/2-lo, L/2+hi]
  python3 -B toy_runs.py dimer L J eta              # imprint of the colour period
Sector runs give the lowest levels, ground-state degeneracy and momenta in
each occupied-number sector (perfect necklace if b | N, one-defect necklace
otherwise).  Each call appends one JSON line per result to results.jsonl.
"""
import sys
import json
import time
import numpy as np
from toy_ed import Toy, lowest, momenta

BASE = dict(t=1.0, J=0.5, K=10.0, Jxy=1.0, Jz=0.5)
EPS2 = [0.0, 0.0]
ETA = {1: None, 2: [0.3, -0.3], 3: [0.3, -0.1, -0.2]}


def out(rec):
    rec['time'] = time.strftime('%Y-%m-%d %H:%M:%S')
    with open('results.jsonl', 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)


def sector(L, N, b, eps, k=4, theta=0.0, eta='default', base=None):
    w = [j % b for j in range(N)]
    n_up = N // 2
    toy = Toy(L, N, w, n_up, b)
    eta = ETA[b] if eta == 'default' else eta
    H = toy.hamiltonian(eps=eps, theta=theta, eta=eta, **(base or BASE))
    ev, vec = lowest(H, k=min(k, toy.dim - 1))
    return toy, H, ev, vec


def run_sectors(Ls, b=2, eps=EPS2, lo=2, hi=2):
    for L in Ls:
        for N in range(max(2, L // 2 - lo), min(L - 1, L // 2 + hi) + 1):
            t0 = time.time()
            toy, H, ev, vec = sector(L, N, b, eps)
            deg = ev[ev < ev[0] + 1e-8]
            P = momenta(toy.translation(), vec[:, :len(deg)], L)
            out(dict(kind='sector', L=L, N=N, b=b, eps=eps, eta=ETA[b], base=BASE, dim=toy.dim,
                     necklace_period=toy.P, defects=int(sum(
                         1 for j in range(N) if toy.w[(j + 1) % N] != (toy.w[j] + 1) % b)),
                     E=[float(e) for e in ev], gs_degeneracy=int(len(deg)),
                     gs_momenta_2pi_over_L=[float(p) for p in P],
                     seconds=round(time.time() - t0, 1)))


def run_dimer(L, J, eta):
    """Squeezed AB-minus-BA alternation of <Sz Sz> induced through the vacancy fluid."""
    N = L // 2
    base = dict(BASE)
    base['J'] = J
    toy, H, ev, vec = sector(L, N, 2, [0.0, 0.0], k=3, eta=[eta, -eta], base=base)
    psi = vec[:, 0]
    p2 = np.abs(psi) ** 2
    spin = toy.spins[toy.si]
    zz = {0: 0.0, 1: 0.0}
    for k in range(N):
        c = int(toy.w[k])
        ck = np.array(toy.w)[(k + toy.r) % N]
        s1 = ((spin >> k) & 1) - 0.5
        s2 = ((spin >> ((k + 1) % N)) & 1) - 0.5
        for col in (0, 1):
            m = ck == col
            zz[col] += float(np.sum(p2[m] * s1[m] * s2[m]))
    out(dict(kind='dimer', L=L, N=N, J=J, eta=eta, gap=float(ev[1] - ev[0]),
             zz_AB=zz[0] / (N / 2), zz_BA=zz[1] / (N / 2),
             alternation=(zz[0] - zz[1]) / (N / 2)))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'sectors':
        run_sectors([int(a) for a in sys.argv[2:]])
    elif mode == 'sectors_b':
        b = int(sys.argv[2]); lo = int(sys.argv[3]); hi = int(sys.argv[4])
        run_sectors([int(a) for a in sys.argv[5:]], b=b, eps=[0.0] * b, lo=lo, hi=hi)
    elif mode == 'dimer':
        run_dimer(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))
