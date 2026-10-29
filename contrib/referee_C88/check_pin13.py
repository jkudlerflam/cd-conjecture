"""(2) Pin conservation, colouring uniqueness and N_matbad <= N_bad, exhaustive.
 (a) Every off-diagonal entry of the coloured 169 x 169 bond (material moves, port,
     vacancy hops, at nonzero g, t_d, lambda, any K) applied on every bond of every
     ring word over the 13 letters, L = 2..6, legal or not, conserves F (PROOF eq. (2)).
     The same for the uncoloured control (12 letters), L = 2..6.
 (b) On the same rings: N_matbad <= N_bad for every word with N > 0; zero-penalty
     words with N > 0 have exactly one colouring per occupancy pattern, and their
     number equals sum_N binom(L, N) tr(A^N) with A the squeezed allowed matrix
     (bijection with legal squeezed words); the all-vacancy sector has 2 colourings.
 (c) Rule-table facts used in the general proof of [H, F] = 0."""
import sys, os, json
from math import comb
import numpy as np
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13

pm = my13.params(3, K=37., musel=.93, field=.021, g=.13, td=.21, phiA=.37, psiA=1.3, theta=-.21)
out = {'conservation': [], 'colouring': []}


def digits(L, d):
    idx = np.arange(d ** L, dtype=np.int64)
    return idx, np.array([(idx // d ** (L - 1 - k)) % d for k in range(L)])


def Fvec(digs, ls):
    L = digs.shape[0]
    vac = np.array([my13.is_vac(l) for l in ls])
    isL1 = np.array([l == 'L1' for l in ls]); ismk = np.array([l in my13.MARKERS for l in ls])
    Fv = np.zeros(digs.shape[1], np.int64)
    for i in range(L):
        # next material letter after site i (cyclic); if none other than itself, it is itself
        nxt = digs[i].copy(); found = np.zeros(digs.shape[1], bool)
        for k in range(1, L + 1):
            dk = digs[(i + k) % L]
            take = (~found) & (~vac[dk])
            nxt[take] = dk[take]; found |= take
        Fv += (isL1[digs[i]] & ismk[nxt]).astype(np.int64)
    return Fv


for coloured in (True, False):
    h, ls = my13.bond(pm, .17, J=.7, t=1.3, mu=.4, coloured=coloured)
    d = len(ls)
    hs = sp.coo_matrix(h - np.diag(np.diag(h)))
    for L in (2, 3, 4, 5, 6):
        idx, digs = digits(L, d)
        Fv = Fvec(digs, ls)
        ntr = 0; viol = 0
        for i in range(L):
            j = (i + 1) % L
            pair = digs[i] * d + digs[j]
            order = np.argsort(pair, kind='stable'); ps = pair[order]
            st = np.searchsorted(ps, np.arange(d * d)); en = np.searchsorted(ps, np.arange(d * d), side='right')
            for r, c in zip(hs.row, hs.col):
                ws = order[st[c]:en[c]]
                if len(ws) == 0:
                    continue
                ar, br = divmod(int(r), d)
                new = ws + (ar - digs[i][ws]) * d ** (L - 1 - i) + (br - digs[j][ws]) * d ** (L - 1 - j)
                ntr += len(ws); viol += int(np.sum(Fv[new] != Fv[ws]))
        rec = dict(coloured=coloured, L=L, words=int(d ** L), directed_transitions=ntr, F_changes=viol)
        out['conservation'].append(rec); print(json.dumps(rec), flush=True)
        if coloured:
            vac = np.array([my13.is_vac(l) for l in ls])
            okp = np.array([[my13.allowed(a, b) for b in ls] for a in ls])
            nb = sum((~okp[digs[i], digs[(i + 1) % L]]).astype(int) for i in range(L))
            # squeezed bad
            nmb = np.zeros(len(idx), int)
            for i in range(L):
                nxt = digs[i].copy(); found = np.zeros(len(idx), bool)
                for k in range(1, L + 1):
                    dk = digs[(i + k) % L]; take = (~found) & (~vac[dk]); nxt[take] = dk[take]; found |= take
                nmb += ((~vac[digs[i]]) & (~okp[digs[i], nxt])).astype(int)
            Nocc = sum((~vac[digs[i]]).astype(int) for i in range(L))
            mask = Nocc > 0
            # occupancy-pattern key with vacancy colour erased
            key = np.zeros(len(idx), np.int64)
            for i in range(L):
                dd = np.where(vac[digs[i]], 11, digs[i]); key = key * 12 + dd
            leg = (nb == 0) & mask
            _, counts = np.unique(key[leg], return_counts=True)
            A = np.array([[my13.allowed(a, b) for b in my13.MAT] for a in my13.MAT], dtype=np.int64)
            pred = sum(comb(L, N) * int(np.trace(np.linalg.matrix_power(A, N))) for N in range(1, L + 1))
            rec = dict(L=L, nonempty_words=int(mask.sum()), violations_Nmatbad_le_Nbad=int(np.sum(nmb[mask] > nb[mask])),
                       zero_penalty_nonempty=int(leg.sum()), colourings_per_pattern=sorted(set(counts.tolist())),
                       predicted_legal_squeezed_count=pred,
                       all_vacancy_zero_penalty=int(np.sum((nb == 0) & ~mask)))
            out['colouring'].append(rec); print(json.dumps(rec), flush=True)

# (c) rule-table facts
R = my13.material_rules(pm, .17)
mk = set(my13.MARKERS)
pairs = [(a, b, c, e) for a, b, c, e, _ in R] + [(c, e, a, b) for a, b, c, e, _ in R]
out['rule_facts'] = dict(
    pin_pair_never_input_or_output=all(not (x == 'L1' and y in mk) for a, b, c, e in pairs for x, y in ((a, b), (c, e))),
    L1_output_right_partners=sorted({e for a, b, c, e in pairs if c == 'L1'}),
    L1_input_right_partners=sorted({b for a, b, c, e in pairs if a == 'L1'}),
    right_letter_becomes_L1=any(e == 'L1' and b != 'L1' for a, b, c, e in pairs),
    left_letter_becomes_L1_from_other=sorted({(a, b, c, e) for a, b, c, e in pairs if c == 'L1' and a != 'L1'}),
    right_marker_type_changes_left_inputs=sorted({a for a, b, c, e in pairs if b in mk and e in mk and b != e}),
    marker_to_nonmarker=any((b in mk) != (e in mk) or (a in mk) != (c in mk) for a, b, c, e in pairs))
print(json.dumps(out['rule_facts'], default=str))
json.dump(out, open(os.path.join(HERE, 'check_pin13.json'), 'w'), indent=1, default=str)
