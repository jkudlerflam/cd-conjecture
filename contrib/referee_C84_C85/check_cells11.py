"""(2) C84 legal cells, the boundary pin, the finite-K vacuum splitting, and the
zero-port selector bound, all from my11.py.
 (a) brute-force zero-penalty words after a marker, count (5p+3)2^(p-1), decoder
     bijection, agreement with the C80 legal words with last bit fixed to 0
     (previous referee's my12.py), marker-free zero-penalty rings;
 (b) pin audit: every bad pair that a single literal move creates from a legal word,
     classified by type; the pin type L1-marker must never appear;
 (c) the order-1/K vacuum-code operator (exact, from the intermediate states),
     Walsh decomposition, comparison with C85 eq. (10), and ED of the full
     zero-excitation cell space at finite K;
 (d) the zero-port selector form bound (C84 eq. (20)) on complete single-cell
     Hilbert spaces including illegal and excited sectors."""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'referee_C80_C82'))
import my11, my12
from my11 import allowed, MARKERS

out = {}
# (a)
cnt = {}
for p in range(1, 8):
    W = my11.legal_words(p)
    dec = [my11.decode(w) for w in W]
    inj = len(set(dec)) == len(dec)
    inv = all(my11.encode_word(*x) == w for x, w in zip(dec, W))
    one_exc = all(sum(my11.excited(l) for l in w) <= 1 for w in W)
    rec = dict(count=len(W), formula=(5 * p + 3) * 2 ** (p - 1), decoder_injective=inj, encode_inverts=inv,
               at_most_one_excited=one_exc)
    if p <= 6:
        old = my12.legal_words(p)
        img = set()
        for w in old:
            j, s, bits = my12.decode(w)
            if bits[-1] == 0:
                img.add(tuple('E' if l == 'E0' else l for l in w))
        rec['equals_C80_legal_with_last_bit_0'] = img == set(W)
    cnt[p] = rec
    print(p, rec, flush=True)
out['legal_counts'] = cnt
ring = {}
for L in range(2, 8):
    zs = [w for w in itertools.product(my11.DATA, repeat=L)
          if all(allowed(w[i], w[(i + 1) % L]) for i in range(L))]
    kinds = set('L' if all(l[0] == 'L' for l in w) else 'R' if all(l[0] == 'R' for l in w) else 'mixed' for w in zs)
    ring[L] = dict(count=len(zs), kinds=sorted(kinds), expected=2 * 2 ** L)
out['marker_free_zero_penalty_rings'] = ring

# (b) pin audit
pm = my11.default_params(4, g=.1, td=.3)
mv = my11.rule_table(my11.move_rules(pm))
types = {}
for p in range(1, 7):
    for w in my11.legal_words(p):
        for left in MARKERS:
            word = (left,) + w
            for i in range(len(word) - 1):
                for a2, b2, val in mv.get((word[i], word[i + 1]), []):
                    z = list(word); z[i] = a2; z[i + 1] = b2
                    if i == 0 and a2 != word[0]:
                        continue
                    bad = [(z[k], z[k + 1]) for k in range(len(z) - 1) if not allowed(z[k], z[k + 1])]
                    if bad:
                        key = ';'.join(f'{a}|{b}' for a, b in bad) + f' via {word[i]}{word[i+1]}->{a2}{b2}'
                        types[key] = types.get(key, 0) + 1
pin_created = sum(v for k, v in types.items() if any(f'L1|{m}' in k.split(' via')[0] for m in MARKERS))
out['single_move_bad_pair_types'] = types
out['pin_pairs_created_from_legal_words'] = pin_created
print('pin created', pin_created, 'types', len(types), flush=True)


# (c) splitting
def vac(p):
    j = np.arange(p + 1)
    return np.sqrt(2 / (p + 2)) * np.sin(np.pi * (j + 1) / (p + 2))


def walsh(f, n):
    """f: dict bits->value. Return max |coefficient| at Walsh weight >= 2."""
    keys = list(itertools.product((0, 1), repeat=n))
    mx = 0.
    for S in keys:
        if sum(S) < 2:
            continue
        c = np.mean([f[a] * (-1) ** sum(s * x for s, x in zip(S, a)) for a in keys])
        mx = max(mx, abs(c))
    return mx


split = {}
pm0 = dict(pm); pm0['K'] = np.inf
for p in range(1, 10):
    n = p - 1
    v = vac(p)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(p)]
    res = {}
    for bits in itertools.product((0, 1), repeat=n):
        code = {my11.encode_word(j, 'v', bits): v[j] for j in range(p + 1)}
        # apply moves (K-independent part), collect illegal outputs
        img = {}
        for w, amp in code.items():
            word = ('E',) + w
            for i in range(len(word) - 1):
                for a2, b2, val in mv.get((word[i], word[i + 1]), []):
                    if i == 0 and a2 != word[0]:
                        continue
                    z = list(word); z[i] = a2; z[i + 1] = b2; z = tuple(z[1:])
                    nb = my11.nbad(z)
                    if nb:
                        img[z] = img.get(z, 0) + val * amp
        res[bits] = -sum(abs(x) ** 2 / my11.nbad(z) for z, x in img.items())  # coefficient of 1/K
    pred = {bits: -sum((1 - bits[r]) * np.sum(v[:r + 1] ** 2) for r in range(n)) for bits in res}
    rec = dict(max_diff_vs_C85_eq10=float(max(abs(res[b] - pred[b]) for b in res)),
               spread_times_K=float(max(res.values()) - min(res.values())),
               walsh_weight_ge2_max=float(walsh(res, n)) if n >= 2 else 0.)
    split[p] = rec
    print('split', p, rec, flush=True)
out['order_1_over_K_splitting'] = split

# ED of zero-excitation cell space (letters L0 L1 R0 R1, marker E)
ed = {}
for p in range(2, 7):
    n = p - 1
    v = vac(p)
    seeds = [w + ('E',) for w in itertools.product(['L0', 'L1', 'R0', 'R1'], repeat=p)]
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(p)]
    Ev = -2 * np.cos(np.pi / (p + 2))
    for K in (1e3, 1e4, 1e5):
        q = dict(pm); q['K'] = K
        states, H, V = my11.closure(seeds, q, bonds, left='E')
        Q = np.array([sum(my11.charge(l) for l in w) for w in states])
        worst = 0.
        for Qv in sorted(set(Q)):
            idx = np.flatnonzero(Q == Qv)
            Hs = H[np.ix_(idx, idx)].toarray()
            e = np.linalg.eigvalsh(Hs)
            bits_here = [b for b in itertools.product((0, 1), repeat=n) if 2 * sum(b) == Qv]
            if not bits_here:
                continue
            got = np.sort(e)[:len(bits_here)]
            # recompute raw coefficients (not stored above for JSON brevity)
            raw = []
            for b in bits_here:
                code = {my11.encode_word(j, 'v', b): v[j] for j in range(p + 1)}
                img = {}
                for w, amp in code.items():
                    word = ('E',) + w
                    for i in range(len(word) - 1):
                        for a2, b2, val in mv.get((word[i], word[i + 1]), []):
                            if i == 0 and a2 != word[0]:
                                continue
                            z = list(word); z[i] = a2; z[i + 1] = b2; z = tuple(z[1:])
                            if my11.nbad(z):
                                img[z] = img.get(z, 0) + val * amp
                raw.append(Ev - sum(abs(x) ** 2 for x in img.values()) / K)
            worst = max(worst, float(np.max(abs(np.sort(raw) - got))) * K ** 2)
        ed.setdefault(p, {})[str(int(K))] = dict(dim=len(states), max_error_times_K2=worst)
    print('ed', p, ed[p], flush=True)
out['splitting_ED_check'] = ed

# (d) selector form bound on complete single cells, musel = 0
sel = []
for k in range(0, 5):
    words = [w + (m,) for w in itertools.product(my11.DATA, repeat=k) for m in MARKERS]
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(k)]
    for K in (6., 10., 100.):
        q = my11.default_params(max(k, 1), D=4., Delta=2.2, g=.1, td=.5, phiA=.3, phiB=0., psiA=1.1,
                                psiB=np.pi / 2, K=K, musel=0., field=0.)
        states, H, V = my11.closure(words, q, bonds, left='E')
        assert len(states) == len(words)
        if len(states) < 3000:
            E = float(np.linalg.eigvalsh(H.toarray())[0])
        else:
            E = float(eigsh(H, k=1, which='SA', tol=1e-12)[0][0])
        Ev = -2 * np.cos(np.pi / (k + 2))
        lo = Ev - max(k - 1, 0) / (K - 4)
        rec = dict(k=k, K=K, dim=len(states), E=E, lower=lo, upper=Ev, holds=bool(lo - 1e-10 <= E <= Ev + 1e-10))
        sel.append(rec)
        print(rec, flush=True)
out['selector_bound_complete_cells'] = sel
json.dump(out, open(os.path.join(HERE, 'check_cells11.json'), 'w'), indent=1)
