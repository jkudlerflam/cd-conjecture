"""Exhaustive checks on all 12^L ring words, L=2..LMAX (vectorized numpy).

  (a) squeezing: N_matbad(w) <= N_bad(w) for every word with material
      (legal physical word => legal squeezed cyclic material word);
  (b) the pin count F (squeezed L_1 -> marker junctions) is conserved by
      every literal transition, legal or not;
  (c) every vacancy hop maps legal words to legal words except the single
      context  o(x)=1, E, V -> x, V, E  (the 'E-context' leak);
  (d) the vacancy-context switch never leaves the legal space;
  (e) count of transitions legal -> illegal by term type.
"""
import sys, os, json
import numpy as np
import bond12 as b

HERE = os.path.dirname(os.path.abspath(__file__))
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = 12
VI = b.IX[b.V]
ALW = np.array([[b.allowed(a, c) for c in b.LET] for a in b.LET])
OUT1 = np.array([b.outgoing(a) for a in b.LET])
ISMARK = np.array([b.marker(a) for a in b.LET])
L1 = b.IX[('L', 1)]
E = b.IX[('e', 0)]

p = b.params(4, K=0., musel=.93, field=.021, flux=.37, theta=-.21, g=.11, td=.05)
h = b.bond(p, lam=.17, t=.4)
off = [(r, c, h[r, c]) for r in range(d * d) for c in range(d * d) if r != c and h[r, c] != 0]


def kind(r, c):
    a, bb = divmod(c, d); a2, b2 = divmod(r, d)
    if VI in (a, bb) and VI in (a2, b2):
        if {a, bb} == {a2, b2}:
            return 'hop'
        return 'vswitch'
    return 'material'


def words(L):
    w = np.indices((d,) * L).reshape(L, -1).T.astype(np.int8)
    return w


def next_material(w):
    """for each position j, the letter of the next material site cyclically
    (excluding j itself unless it is the only material letter), -1 if none."""
    N, L = w.shape
    nxt = np.full((N, L), -1, np.int16)
    for j in range(L):
        found = np.zeros(N, bool)
        for k in range(1, L + 1):
            col = w[:, (j + k) % L]
            new = (~found) & (col != VI)
            nxt[new, j] = col[new]
            found |= new
    return nxt


def F_and_matbad(w):
    nxt = next_material(w)
    ismat = w != VI
    has = nxt >= 0
    F = np.sum(ismat & (w == L1) & has & ISMARK[np.where(has, nxt, 0)], axis=1)
    mb = np.sum(ismat & has & ~ALW[w, np.where(has, nxt, 0)], axis=1)
    return F, mb


def nbad(w):
    L = w.shape[1]
    return np.sum(~ALW[w, np.roll(w, -1, axis=1)], axis=1)


res = {}
for L in range(2, LMAX + 1):
    w = words(L)
    F, mb = F_and_matbad(w)
    nb = nbad(w)
    anymat = np.any(w != VI, axis=1)
    sq = int(np.sum(anymat & (mb > nb)))
    rec = dict(L=L, words=int(len(w)), squeeze_violations=sq, legal_words=int(np.sum(nb == 0)),
               F_violations=0, leaks={'hop': 0, 'vswitch': 0, 'material': 0}, hop_leaks_not_E_context=0,
               transitions=0)
    for i in range(L):
        j = (i + 1) % L
        if L == 2 and i == 1:
            pass
        pair = w[:, i].astype(np.int32) * d + w[:, j]
        for r, c, val in off:
            m = np.flatnonzero(pair == c)
            if m.size == 0:
                continue
            z = w[m].copy()
            z[:, i], z[:, j] = divmod(r, d)
            Fz, _ = F_and_matbad(z)
            rec['F_violations'] += int(np.sum(Fz != F[m]))
            rec['transitions'] += int(m.size)
            src_legal = nb[m] == 0
            dst_illegal = nbad(z) > 0
            leak = src_legal & dst_illegal
            k = kind(r, c)
            rec['leaks'][k] += int(np.sum(leak))
            if k == 'hop' and np.any(leak):
                a, bb = divmod(c, d)
                left = w[m[leak], (i - 1) % L]
                ok = (a == E) & (bb == VI) & (OUT1[left] == 1)
                rec['hop_leaks_not_E_context'] += int(np.sum(~ok))
    assert rec['squeeze_violations'] == 0 and rec['F_violations'] == 0
    assert rec['leaks']['vswitch'] == 0 and rec['hop_leaks_not_E_context'] == 0
    print(json.dumps(rec), flush=True)
    res[L] = rec
json.dump(res, open(os.path.join(HERE, 'ring_checks.json'), 'w'), indent=1)
