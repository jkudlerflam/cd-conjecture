#!/usr/bin/env python3
"""Rotating-frame (spiral gauge) decomposition of a nearest-neighbour qubit bond.

H = sum_j sigma_j^T J sigma_{j+1} + h . sigma_j  (Pauli matrices).
Define tau_j = R^j sigma_j with R = R_n(theta).  Then
sigma_j^T J sigma_{j+1} = tau_j^T R^j K R^{-j} tau_{j+1},  K = J R^T.
The part of K invariant under rotations about n is uniform in the new frame,
  P_n(K) = a (1 - n n^T) + b n n^T + c [n]_x ,
the rest oscillates with wave vectors theta and 2 theta.  We pick theta(n) so
that c = 0 (no residual DM along n) and then minimise the oscillating weight
||K - P_n(K)||_F / ||K||_F over the axis n.  The uniform model is an XXZ chain
a(xx+yy) + b zz in a field (h.n) along z; Delta_eff = b/|a| in the
antiferromagnetic convention.
"""
import json, sys
import numpy as np
from scipy.optimize import minimize

def rot(n, th):
    n = np.asarray(n, float); n = n/np.linalg.norm(n)
    X = np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
    return np.eye(3)+np.sin(th)*X+(1-np.cos(th))*X@X, X

def proj(M, n):
    n = np.asarray(n, float); n = n/np.linalg.norm(n)
    _, X = rot(n, 0)
    P = np.outer(n, n); Q = np.eye(3)-P
    a = np.trace(Q@M)/2; b = n@M@n; c = -np.trace(X@M)/2
    return a, b, c, a*Q+b*P+c*X

def best_theta(J, n):
    # c(theta) for K = J R^T; solve c = 0 on a grid, pick the root with smallest residual
    ths = np.linspace(-np.pi, np.pi, 721)
    cs = np.array([proj(J@rot(n,t)[0].T, n)[2] for t in ths])
    out = []
    for i in range(len(ths)-1):
        if cs[i] == 0 or cs[i]*cs[i+1] < 0:
            lo, hi = ths[i], ths[i+1]
            for _ in range(60):
                mid = (lo+hi)/2
                if proj(J@rot(n,lo)[0].T, n)[2]*proj(J@rot(n,mid)[0].T, n)[2] <= 0: hi = mid
                else: lo = mid
            out.append((lo+hi)/2)
    return out

def resid(J, n):
    best = (np.inf, None)
    for t in best_theta(J, n):
        K = J@rot(n,t)[0].T
        a, b, c, P = proj(K, n)
        r = np.linalg.norm(K-P)/np.linalg.norm(K)
        if r < best[0]: best = (r, t)
    return best

def analyse(J, h):
    J = np.asarray(J, float); h = np.asarray(h, float)
    best = (np.inf, None, None)
    for th0 in np.linspace(0, np.pi, 13)[:-1]:
        for ph0 in np.linspace(0, 2*np.pi, 13)[:-1]:
            n = np.array([np.sin(th0)*np.cos(ph0), np.sin(th0)*np.sin(ph0), np.cos(th0)])
            r, t = resid(J, n)
            if r < best[0]: best = (r, t, n)
    # polish
    def f(x):
        n = np.array([np.sin(x[0])*np.cos(x[1]), np.sin(x[0])*np.sin(x[1]), np.cos(x[0])])
        return resid(J, n)[0]
    n0 = best[2]; x0 = [np.arccos(np.clip(n0[2],-1,1)), np.arctan2(n0[1], n0[0])]
    res = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-10))
    x = res.x; n = np.array([np.sin(x[0])*np.cos(x[1]), np.sin(x[0])*np.sin(x[1]), np.cos(x[0])])
    r, t = resid(J, n)
    K = J@rot(n,t)[0].T
    a, b, c, P = proj(K, n)
    hn = h@n; hperp = np.linalg.norm(h-hn*n)
    return dict(axis=n.tolist(), theta=t, a=a, b=b, c=c, osc=r, Delta_eff=b/abs(a),
                h_par=hn/abs(a), h_perp=hperp/abs(a), scale=abs(a))

if __name__ == '__main__':
    pars = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'candidates_params.json'))
    out = {}
    print(f"{'id':10s} {'osc':>6s} {'theta/pi':>8s} {'Delta_eff':>9s} {'h_par':>7s} {'h_perp':>7s} axis")
    for k, p in pars.items():
        r = analyse(p['J'], p['field']); out[k] = r
        print(f"{k:10s} {r['osc']:6.3f} {r['theta']/np.pi:8.4f} {r['Delta_eff']:9.4f} {r['h_par']:7.3f} {r['h_perp']:7.3f} {np.round(r['axis'],3)}")
    json.dump(out, open('frame_analysis.json', 'w'), indent=1)
