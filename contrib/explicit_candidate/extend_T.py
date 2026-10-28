"""Follow-up checks on the second-order hopping T of the 23-letter bond.
  (i)  Richardson test of the cancellation-free g^2 coefficient (bond23 method C)
       against the full-in-g one-particle T (method B) at g and g/2, n = 3, 4.
  (ii) g^2 coefficient at n = 5, 6 for several express energies D_x, with the
       express phase 0.4 and pi/(2n): eigenvalues of the Hermitian part.
  (iii) Reality of the Hilbert-Schmidt-projected T of ycell at t_d = 0, n = 2.
  (iv) Regenerate logical_chain.analyze at the README working points n = 1, 2.
Usage: python3 extend_T.py  ->  extend_T_results.json"""
import os, sys, json, time
import numpy as np
from bond23 import params, onebody_T, herm_eigs, HERE

out = {}
# (i)
rich = []
for n in (3, 4):
    g = 0.05
    TB1, _, _ = onebody_T(n, params(gp=g, gm=g), 'all')
    TB2, _, _ = onebody_T(n, params(gp=g / 2, gm=g / 2), 'all')
    TC, _, _ = onebody_T(n, params(gp=g, gm=g), 'g2')       # coefficient of g^2
    T2_rich = (4 * TB2 / (g / 2) ** 2 - TB1 / g ** 2) / 3
    rich.append(dict(n=n, rel_diff_g2_vs_richardson=float(np.abs(T2_rich - TC).max() / np.abs(TC).max()),
                     rel_diff_g2_vs_fullg_at_g=float(np.abs(TB1 / g ** 2 - TC).max() / np.abs(TC).max()),
                     rel_diff_g2_vs_fullg_at_g_half=float(np.abs(TB2 / (g / 2) ** 2 - TC).max() / np.abs(TC).max())))
    print(rich[-1], flush=True)
out['richardson'] = rich

# (ii)
ext = []
for n in (2, 3, 4, 5, 6):
    for Dx in (4.0, 10.0, 20.0, 40.0):
        for k0 in (0.4, np.pi / (2 * n)):
            t0 = time.time()
            P = params(Dx=Dx, DM=Dx, phx={'0': k0, '+': 0.0, '-': 0.0})
            TC, E0, res = onebody_T(n, P, 'g2')
            ev, im, sc = herm_eigs(TC)
            rec = dict(n=n, Dx=Dx, phx0=k0, T2_scale=sc, herm_eigs_over_T=list(ev),
                       imag_after_phase=im, min_over_max=float(ev.min() / ev.max()),
                       asym=float(np.linalg.norm((TC - TC.T)) / np.linalg.norm(TC)),
                       seconds=time.time() - t0)
            print(json.dumps(rec), flush=True)
            ext.append(rec)
out['g2_coefficient'] = ext

# (iii) and (iv)
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_spatial'))
import ycell
from logical_chain import analyse, IR_PARAMS
Py = ycell.params(td=0.0, **IR_PARAMS)
c = ycell.build_cell(2, Py); E0, G, gap = ycell.code(c); H2, k = ycell.second_order(c, G, E0)
L = ycell.logical(H2, 2)
T = L['T']; R = T * np.exp(-1j * np.angle(T[0, 0]))
T1 = ycell.one_particle_T(H2, 2)
out['hs_projection_n2_td0'] = dict(imag_after_phase=float(np.abs(R.imag).max() / np.abs(T).max()),
                                   mu_over_T=float(np.abs(L['muA'] + L['muB']).max() / np.abs(T).max()),
                                   hs_vs_onebody=float(np.abs(T - T1).max() / np.abs(T).max()),
                                   nonGaussian_over_T=float(L['res'] / np.abs(T).max()))
print(out['hs_projection_n2_td0'], flush=True)
regen = []
for (n, td) in ((1, 1e-3), (2, 1e-3)):
    o = analyse(n, td)
    regen.append({k_: o[k_] for k_ in ['n', 'td', 'cell_gap', 'poles', 'T_scale', 'mu_over_T',
                                        'nonGaussian_over_T', 'herm_part_eigs_over_T', 'roots_per_band',
                                        'centres', 'delta4', 'delta2', 'vmin_over_T']})
    print(json.dumps(regen[-1], default=float), flush=True)
out['logical_chain_regenerated'] = regen
with open(os.path.join(HERE, 'extend_T_results.json'), 'w') as f:
    json.dump(out, f, indent=1, default=float)
