"""Central-switch cell with a small left/right arm asymmetry (xL != 1) and a
phase phiT_L on carrier-T hops in the left arm: is T non-Hermitian, are the
Fermi-sea centres displaced, and is c = n kept?"""
import sys, json
import numpy as np, scipy.linalg as la, scipy.sparse.linalg as sla
from model import Chain, DEFAULT
from model2 import LaneChain, segment_parts, assemble_seg
from onebody import blocks, max_c
from switch_phase import fermi_points

def h_precise(ch, P, p):
    H0, V = assemble_seg(P, p)
    H00, _ = assemble_seg(P, dict(p, tauA=0.0, tauB=0.0))
    C = (H0 - H00).tocsr()
    dm = ch.data_modes
    vac = [c for c in ch.code if all(((ch.states[c] >> q) & 1) == 0 for q in dm)][0]
    cols, signs = [vac], [1.0]
    for q in dm:
        sg, new = Chain.apply([(q, 1)], [], ch.states[vac]); cols.append(ch.index[new]); signs.append(sg)
    mask = np.ones(ch.dim, bool); mask[ch.code] = False; Q = np.where(mask)[0]
    VQP = V[Q][:, cols].toarray() * np.array(signs)[None, :]
    X0 = sla.spsolve(H00[Q][:, Q].tocsc(), VQP)
    X1 = sla.spsolve(H0[Q][:, Q].tocsc(), C[Q][:, Q] @ X0)
    Dm = VQP.conj().T @ X1; Dm = 0.5 * (Dm + Dm.conj().T)
    return Dm[1:, 1:] - Dm[0, 0].real * np.eye(len(dm))

m = int(sys.argv[1]); tau = float(sys.argv[2])
ch = LaneChain('switch', m, 2, K=2, max_data=3, one_per_cell=True)
P = segment_parts(ch)
out = []
CASES = [(1.0, 0.0, 0.0), (1.0, 0.3, 0.0), (1.02, 0.0, 0.0), (1.02, 0.3, 0.0), (1.05, 0.3, 0.0), (1.0, 0.0, 1e-4), (1.05, 0.0, 1e-4), (1.05, 0.0, 1e-3)]
for (xL, ph, phU) in CASES:
    p = dict(DEFAULT, wD=-0.3, tauA=tau, tauB=tau, xL=xL, phiT_L=ph, phiT=phU)
    h = h_precise(ch, P, p)
    T, M = blocks(h, m)
    pts = fermi_points(T, M)
    cen = []
    for q in pts:
        if len(q) == 2:
            a, b = q; c = 0.5 * (a + b)
            if b - a > np.pi:
                c = c + np.pi if c < 0 else c - np.pi
            cen.append(c)
    r = dict(m=m, tau=tau, xL=xL, phiT_L=ph, phiT_all=phU, c_max=max_c(T, M, nk=4000)[0],
             T_herm_defect=float(abs(T - T.conj().T).max() / abs(T).max()),
             n_fermi_points_mu0=sum(len(q) for q in pts), centres_mu0=cen)
    print(json.dumps(r), flush=True); out.append(r)
json.dump(out, open(f'switch_break_m{m}.json', 'w'), indent=1)
