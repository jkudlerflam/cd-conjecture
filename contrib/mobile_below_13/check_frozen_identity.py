"""Exact identity: at t=0 (vacancy positions frozen) the block of the d=12 ring
in which every vacancy sits in an inter-cell gap (right after a marker) equals
the frozen eleven-letter ring on the squeezed word, with the port deleted on
every bond that carries a vacancy.  Checked entry by entry, including all
illegal states of the closure and the K penalty."""
import os, json
import numpy as np
import bond12 as b
from closure import closure, squeeze, VI

HERE = os.path.dirname(os.path.abspath(__file__))
L0, L1, R0, R1, E = ('L', 0), ('L', 1), ('R', 0), ('R', 1), ('e', 0)
DA, XB = ('d', 'A'), ('x', 'B')
V = b.V
cases = [
    ('p3_two_cells_one_vacancy', [L1, R0, R0, E, V, L1, L0, R0, E], 3),
    ('p3_two_cells_two_vacancies_one_gap', [L1, R0, R0, E, V, V, L0, R0, R1, E], 3),
    ('p3_two_cells_both_gaps', [L1, R0, R0, E, V, L1, L0, R0, E, V], 3),
    ('p2_three_cells_one_vacancy', [L1, R0, E, V, R0, R1, E, L1, R0, E], 2),
    ('p3_token_next_to_vacancy', [L1, R0, R0, E, V, DA, R1, R0, E], 3),
    ('p3_xtoken_next_to_vacancy', [L1, R0, R0, E, V, XB, R1, R0, E], 3),
]
out = []
for name, seed, p_ in cases:
    p = b.params(p_, K=7.3, musel=.93, field=.021, flux=.37, theta=-.21, g=.11, td=.05, D=5.)
    lam = .17
    Lm = len(seed)
    hp = b.bond(p, lam=lam, t=0.); h0 = b.bond(p, lam=0., t=0.)
    # mobile ring at t=0: all bonds carry the port (it can only act on material pairs)
    pK = dict(p); pK['K'] = 0.
    Hm, penm, mbm, sm = closure(seed, [b.bond(pK, lam=lam, t=0.)] * Lm)
    # frozen squeezed ring: port off on bonds that were cut by a vacancy
    sq = tuple(x for x in seed if x != V); Lf = len(sq)
    pos = [i for i, x in enumerate(seed) if x != V]
    cut = [any(seed[k % Lm] == V for k in range(pos[i] + 1, pos[i + 1] if i + 1 < Lf else pos[0] + Lm))
           for i in range(Lf)]
    bonds_f = [b.bond(pK, lam=(0. if cut[i] else lam), t=0.) for i in range(Lf)]
    Hf, penf, mbf, sf = closure(sq, bonds_f)
    # map mobile states to squeezed words: vacancy positions must be frozen
    vac = [i for i, x in enumerate(seed) if x == V]
    assert all(all(s[i] == VI for i in vac) for s in sm)
    fidx = {s: i for i, s in enumerate(sf)}
    order = [fidx[tuple(x for x in s if x != VI)] for s in sm]
    assert len(order) == len(sf) == len(set(order))
    Hfp = Hf[order][:, order]
    err_off = float(abs(Hm - Hfp).max())
    err_pen = float(np.max(abs(penm - penf[order])))
    err_mb = float(np.max(abs(mbm - mbf[order])))
    rec = dict(case=name, dimension=len(sm), cut_bonds=int(sum(cut)), max_entry_error=err_off,
               penalty_count_error=err_pen, squeezed_bad_error=err_mb)
    assert err_off < 1e-14 and err_pen == 0 and err_mb == 0
    print(json.dumps(rec), flush=True)
    out.append(rec)
# control: without the vacancy-context switch the identity fails
p = b.params(3, K=0., musel=.93, field=.021, flux=.37, theta=-.21, g=.11, td=.05, D=5.)
seed = [L1, R0, R0, E, V, DA, R1, R0, E]
Hm, *_ , sm = closure(seed, [b.bond(p, lam=.17, t=0., vswitch=False)] * 9)
sq = tuple(x for x in seed if x != V)
Hf, *_ , sf = closure(sq, [b.bond(p, lam=(0. if i == 3 else .17), t=0.) for i in range(8)])
fidx = {x: i for i, x in enumerate(sf)}
order = [fidx[tuple(x for x in s_ if x != VI)] for s_ in sm]
out.append(dict(case='control_without_vswitch', mobile_dimension=len(sm), frozen_dimension=len(sf),
                max_entry_error=float(abs(Hm - Hf[order][:, order]).max())))
print(out[-1])
json.dump(out, open(os.path.join(HERE, 'frozen_identity.json'), 'w'), indent=1)
