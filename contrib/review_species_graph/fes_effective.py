#!/usr/bin/env python3
"""Finite-entanglement scaling of a scaled-down logical chain of the
species-graph candidate, in cell coordinates.

Each cell carries one length bit eta (long/short cell, the C8 length fluid)
and n_pay payload occupations n_a (logical fermion flavours).  All species
are hard-core particles with their own Jordan-Wigner string, as in C13's
numerical subclass; this retains every U(1) that the species torus protects
at quartic order (the graph forbids flavour conversion below six fields).

First-order reference (the proposal's MICROSCOPIC_HAMILTONIAN.md section 5):
  length: -J(eta+_j eta-_{j+1} + h.c.), J = 1/2, zero chemical potential (nu = 1/2)
  payload a: -t_a(n+_j n-_{j+1} + h.c.) + n_a,  t_a = 1/(2 cos k_a)
Interacting runs add, with random signs and overall strength g, every
two-site term the proposal lists as generated at finite coupling:
length-dependent payload hopping, payload-dependent marker hopping,
correlated marker+payload hops, and all density-density couplings
(on-site and nearest neighbour, within and between species).

Estimators (both calibrated on the free point, where c = 1 + n_pay exactly):
  c_xi  = 6 dS/d ln xi   over consecutive chi
  c_chi from dS/d ln chi = 1/(sqrt(12/c)+1)  (Pollmann et al. 2009)

Uses ../d2_followup_220_1282/vumps.py (numpy + scipy) without modifying it.
usage: fes_effective.py NAME --npay 1 --g 0.2 --seed 1 --chis 16 24 32 48 64 96 128
"""
import os, sys, json, time, argparse, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
sys.path.insert(0, HERE)
import numpy as np
import vumps as V

sp = np.array([[0., 0.], [1., 0.]])   # raising: |1><0|
sm = sp.T
nn = np.diag([0., 1.])
I2 = np.eye(2)


def site_op(ops, k, nsp):
    """Operator on the nsp-species site with ops[k] on species k."""
    out = np.array([[1.]])
    for s in range(nsp):
        out = np.kron(out, ops.get(s, I2))
    return out


def bond(npay, g, seed, kF, Jlen=0.5, mu_len=0.0):
    nsp = 1 + npay
    d = 2 ** nsp
    rng = np.random.default_rng(seed)
    r = lambda: g * rng.uniform(-1, 1)
    one = lambda s, O: site_op({s: O}, None, nsp)
    Id = np.eye(d)
    eta_n, eta_p, eta_m = one(0, nn), one(0, sp), one(0, sm)
    n_ = [one(1 + a, nn) for a in range(npay)]
    p_ = [one(1 + a, sp) for a in range(npay)]
    m_ = [one(1 + a, sm) for a in range(npay)]
    t = [1 / (2 * np.cos(k)) for k in kF]
    f = [k / np.pi for k in kF]
    h = np.zeros((d * d, d * d))
    K = lambda A, B: np.kron(A, B)
    # length fluid, payload-dependent amplitude
    beta = r()
    ntot = sum(n_) if npay else 0 * Id
    hopL = K(eta_p, eta_m) + K(eta_m, eta_p)
    h += -Jlen * hopL - Jlen * beta * 0.5 * ((K(ntot, Id) + K(Id, ntot) - 2 * sum(f) * K(Id, Id)) @ hopL
                                             + hopL @ (K(ntot, Id) + K(Id, ntot) - 2 * sum(f) * K(Id, Id)))
    for a in range(npay):
        alpha = r()
        hopP = K(p_[a], m_[a]) + K(m_[a], p_[a])
        lenfac = K(eta_n, Id) + K(Id, eta_n) - K(Id, Id)
        h += -t[a] * hopP - t[a] * alpha * 0.5 * (lenfac @ hopP + hopP @ lenfac)
        h += 0.5 * (K(n_[a], Id) + K(Id, n_[a]))                    # chemical potential +1 per particle
        h += r() * K(n_[a], n_[a])                                    # same-flavour nn density
        h += r() * (K(eta_n, n_[a]) + K(n_[a], eta_n))                # length-payload nn density
        h += r() * 0.5 * (K(eta_n @ n_[a], Id) + K(Id, eta_n @ n_[a]))  # on-site length-payload
        gam = r()                                                     # correlated marker+payload hop
        cor = K(eta_p @ p_[a], eta_m @ m_[a])
        h += gam * (cor + cor.T)
        for b in range(a + 1, npay):
            h += r() * (K(n_[a], n_[b]) + K(n_[b], n_[a]))           # inter-flavour nn density
            h += r() * 0.5 * (K(n_[a] @ n_[b], Id) + K(Id, n_[a] @ n_[b]))
    h += r() * K(eta_n, eta_n)
    h += -mu_len * 0.5 * (K(eta_n, Id) + K(Id, eta_n))
    assert np.allclose(h, h.T)
    dens = [eta_n] + n_
    return h, d, dens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--npay', type=int, default=1)
    ap.add_argument('--g', type=float, default=0.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--kF', type=float, nargs='+', default=None, help='payload Fermi momenta / pi')
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64])
    ap.add_argument('--maxit', type=int, default=600)
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    kF = [x * np.pi for x in (a.kF or [0.27, 0.29, 0.31][:a.npay])]
    h, d, dens = bond(a.npay, a.g, a.seed, kF)
    W = V.mpo_from_term(h, 2, d=d)
    out = a.out or os.path.join(HERE, 'data', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), c_free=1 + a.npay, results=[])
    t0 = time.time()

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['densities'] = V.local_expect(psi, dens).ravel().tolist()
        rec['S'] = [float(x) for x in rec['S']]
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} dens={np.round(rec['densities'], 4).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    V.run([W], 1, a.chis, d=d, tol=1e-9, maxit=a.maxit, seed=a.seed, mode=a.mode,
          log=lambda rec: None, on_chi=save)


if __name__ == '__main__':
    main()
