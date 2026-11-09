"""C166 task 1 (Claude Opus 5.5 agent): same Hamiltonian in dmrg.py (open chains) and vumps.py/fesx.py (cells).
(a) The 2-site operator encoded by fesx.mpo_cell equals est.bond_from_JB(J, f) plus the other half of the field.
(b) The J used in C163's fes/T40.json and fes/T50.json equals est.J_f(v) of hi_T_40 and dist_T_50.
(c) Open-chain DMRG (dmrg_c.chain_bond + dmrg.py) against sparse ED at L = 8, 12, 16.
(d) Periodic ED e0/L at even L = 8..20 against the VUMPS e0 (two-site, chi 192, C163).
(e) A 2-site-cell VUMPS at chi = 16 run here is compared with C163's chi = 16 value.
usage: python3 ed_check.py > ed_check.out"""
import json, sys, time
import numpy as np
import c166lib as X, est as E, fesx as F, vumps as V, dmrg_c as DC

names = ['hi_T_40', 'dist_T_50', 'cand_60', 'cand_61', 'KI22', 'KI26', 'N34']
out = {}
for nm in names:
    J, f = X.JB(nm)
    Ws = F.mpo_cell([J, J], [f, f])
    D = Ws[0].shape[0]
    H2 = sum(np.kron(Ws[0][0, w], Ws[1][w, D - 1]) for w in range(D))
    hb = E.bond_from_JB(J, f) + sum(f[a] * (np.kron(E.PAULI[a], X.I2) + np.kron(X.I2, E.PAULI[a])) / 2 for a in range(3))
    print(f'(a) {nm:10s} |H2_mpo - h_bond - f/2 field| = {np.max(np.abs(H2 - hb)):.1e}   v = {np.round(X.load_vec(nm), 4).tolist()}')
for nm, fn in (('hi_T_40', 'T40'), ('dist_T_50', 'T50')):
    d = json.load(open(f'../work_C163/typhon_out/fes/{fn}.json'))
    J, f = X.JB(nm)
    print(f'(b) {nm}: |J_fes - J_f(v)| = {np.max(np.abs(np.array(d["args"]["J"]).reshape(3, 3) - J)):.1e}, B_fes = {d["args"]["B"]}')
for nm in names:
    J, f = X.JB(nm)
    for L in (8, 12, 16):
        e_ed = X.ed_e0(J, f, L)[0]
        r, h, v = X.dmrg_edges(nm, L, chis=(16, 32, 64, 64))
        print(f'(c) {nm:10s} L={L:2d} E_ED={e_ed:.12f} E_DMRG={r["E"]:.12f} diff={r["E"] - e_ed:+.1e}', flush=True)
        out.setdefault(nm, {})[f'open{L}'] = [float(e_ed), float(r['E'])]
for nm in ('hi_T_40', 'dist_T_50'):
    J, f = X.JB(nm)
    ev = []
    for L in (8, 10, 12, 14, 16, 18, 20):
        t = time.time()
        e = X.ed_e0(J, f, L, periodic=True)[0] / L
        ev.append((L, float(e)))
        print(f'(d) {nm:10s} ring L={L:2d} e0/L={e:.10f} ({time.time() - t:.0f}s)', flush=True)
    out[nm]['ring'] = ev
    fn = 'T40' if nm == 'hi_T_40' else 'T50'
    d = json.load(open(f'../work_C163/typhon_out/fes/{fn}.json'))
    print(f'    VUMPS e0 (n=2, chi 192, C163) = {d["results"][-1]["e"]:.10f}')
    # (e) re-run chi 8..16 two-site VUMPS here
    W = F.mpo_cell([J, J], [f, f])
    psi, res = V.run(W, 2, [8, 12, 16], mode='damped', maxit=600, log=lambda r: None)
    print(f'(e) {nm}: VUMPS here chi=16 e={res[-1]["e"]:.10f}; C163 chi=16 e={d["results"][2]["e"]:.10f}', flush=True)
    out[nm]['vumps16_here'] = res[-1]['e']; out[nm]['vumps16_c163'] = d['results'][2]['e']
json.dump(out, open('ed_check.json', 'w'), indent=1)
