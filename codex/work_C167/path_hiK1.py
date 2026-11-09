"""C167 (Claude Opus 5.5, 2026-11-08): two-site VUMPS along the harmonic-amplitude path of hi_K_1.

J(lam) = (S0 + lam (S - S0)) O,  f(lam) = (m.f) m + lam f_perp,  with J = S O the polar decomposition of hi_K_1
(C159 stage-B point, stratum K_I).  lam = 0 is the zero-field XXZ parent (Delta = -0.154), lam = 1 is hi_K_1.
Optional --scale1 and --scale2 rescale the charge-1 (field) and charge-2 parts separately (Proposition R test).
usage: python3 path_hiK1.py LAM [--s1 1] [--s2 1] --chis 8 16 24 32 --out loc/p_LAM.json [--npz state.npz]
(an existing --npz is used as the initial state, also across different LAM: hysteresis test)
"""
import sys, os, json, argparse, subprocess, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import est as E, tk

V_HIK1 = [-0.19943689208151885, 0.5346005217434073, -0.07328358510629586, 0.0, 0.0, -0.4444987569968905,
          -0.16400859080229255, 0.6667673968887314, 0.0]


def bond(lam, s1=1.0, s2=1.0):
    v = np.asarray(V_HIK1); v = v / np.linalg.norm(v)
    J, f = E.J_f(v)
    p = tk.spiral_parent(J, f)
    S, O, m = p['S'], p['O'], p['m']
    P = np.eye(3) - np.outer(m, m); mm = np.outer(m, m)
    S0 = p['s_m'] * mm + p['s_perp'] * P
    Sq1 = P @ S @ mm + mm @ S @ P
    Sq2 = P @ S @ P - p['s_perp'] * P
    Jl = (S0 + lam * s1 * Sq1 + lam * s2 * Sq2) @ O
    fl = mm @ f + lam * s1 * (P @ f)
    return Jl, fl


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('lam', type=float)
    ap.add_argument('--s1', type=float, default=1.0)
    ap.add_argument('--s2', type=float, default=1.0)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 16, 24, 32])
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--out', required=True)
    ap.add_argument('--npz', default=None)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    J, f = bond(a.lam, a.s1, a.s2)
    cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fesx.py'), 'gen',
           '--J'] + [repr(float(x)) for x in J.ravel()] + ['--B'] + [repr(float(x)) for x in f] + \
          ['--n', str(a.n), '--chis'] + [str(c) for c in a.chis] + ['--out', a.out, '--seed', str(a.seed)] + (['--npz', a.npz] if a.npz else [])
    subprocess.run(cmd, check=True)
