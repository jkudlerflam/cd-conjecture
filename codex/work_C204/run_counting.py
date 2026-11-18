"""Independent check of Theorem 4.11 / Corollaries 4.12-4.13 on true Bethe roots."""
import numpy as np, json, sys, time
from c204lib import *

def analyse(eta, L, n, cont, hf):
    K = Kern(eta); nu = n / L; Q = cont.Q
    I = np.arange(1, n + 1) - (n + 1) / 2
    # initial guess: Z0(x_j) = I_j / L, Z0 on a grid
    ug = np.concatenate([-Q - np.geomspace(12, 1e-3, 300), Q * np.linspace(-1, 1, 1201), Q + np.geomspace(1e-3, 12, 300)])
    ug = np.unique(ug)
    Zg = cont.Z0(ug) - nu * K.q / 2 * 0  # Z0 as defined by p/2pi + int Theta r  (equals PROOF's Z0)
    # PROOF: Z0 = p/(2pi) + T*N0 - q nu/2 = p/(2pi) + int Theta(u-v) r dv
    x0 = np.interp(I / L, Zg, ug)
    x, res, it = bethe_solve(eta, L, n, x0=x0)
    N0 = cont.N0(x)
    j = np.arange(1, n + 1)
    KL = float(np.max(np.maximum(np.abs(j - L * N0), np.abs(j - 1 - L * N0))))
    # signed: lower deficit (N_x below N_0) and upper excess
    low = float(np.max(L * N0 - j))      # >0 means N_x(v-) ... left-limit deficit
    up = float(np.max(j - 1 - L * N0 + 1))  # j - L N0
    rQ = float(cont.s(np.array([Q]))[0])
    lamQ = eta * L * rQ
    # sigma_L vs s on a set of points
    upts = np.unique(np.concatenate([x, Q * np.linspace(-1, 1, 801), Q + np.linspace(0, 3, 200), (x[:-1] + x[1:]) / 2]))
    sig = K.a1(upts) + K.T(upts[:, None] - x[None, :]).sum(1) / L
    sv = cont.s(upts)
    DL = float(eta * L * np.max(np.abs(sig - sv)))
    inside = np.abs(upts) <= Q
    lam = eta * L * sv
    shell = inside & (lam >= 2)
    minratio = float(np.min(sig[shell] / sv[shell])) if shell.any() else None
    # energy
    a = PI * np.sin(eta)
    eint = cont.h * np.sum(cont.w * (-a * K.a1(cont.v)) * cont.r)
    dE = float(np.sum(-a * K.a1(x)) - L * eint)
    out = dict(eta=eta, onepD=1 - np.cos(eta), L=L, n=n, nu=nu, Q=Q, hfac=hf, etaL=eta * L, lamQ=lamQ,
               KL=KL, DL=DL, minratio_shell=minratio, LdE=L * dE, res=res, it=it, xmax=float(x[-1]))
    return out

if __name__ == '__main__':
    onepDs = [float(s) for s in sys.argv[1].split(',')]
    hf = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    cases = []
    for L in [80, 320, 1280, 4000]:
        for nu in [0.05, 0.2, 0.35, 0.45]:
            cases.append((L, int(round(nu * L))))
    cases += [(4000, 3), (1280, 5), (320, 2), (1001, 500), (4001, 2000), (257, 128), (999, 37)]
    results = []
    for x in onepDs:
        eta = float(np.arccos(1 - x))
        Qg = 1.0
        for (L, n) in cases:
            t = time.time()
            cont = continuum_for_nu(eta, n / L, Qguess=Qg, hfac=hf)
            r = analyse(eta, L, n, cont, hf)
            r['time'] = time.time() - t
            results.append(r)
            print(json.dumps(r), flush=True)
    json.dump(results, open('counting_%s_h%d.json' % (sys.argv[1].replace(',', '_'), hf), 'w'), indent=1)
