"""(4) C85 bond and colouring.
 (a) my14.py (from C85 PROOF) vs bond14.py entry by entry, coloured and uncoloured;
     Hermiticity; charges Q, occupied number N, marker number M.
 (b) the three-state automaton reproduces C84's 121-pair table.
 (c) minimality: boolean rank of the allowed relation (fewest rectangles A x B inside
     the allowed set covering a required set of pairs), for three required sets:
     all allowed pairs; pairs adjacent in legal rings of cells with p >= 1 data sites;
     pairs adjacent in rings of cells of one fixed length p = 6 (same as p = 14).
     Fooling sets of size three. Number of deterministic three-colour tables.
 (d) exhaustive rings of 3 to 5 sites over 14 letters: unique zero-penalty colouring
     for each legal squeezed material word and occupancy pattern, and
     N_material_bad <= N_bad.
Usage: python3 check_lift14.py <scratch copy of codex/work_C85>"""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my11, my14
W = sys.argv[1]
sys.path.insert(0, W)
import bond14  # object under test

out = {}
name = {('L', 0): 'L0', ('L', 1): 'L1', ('R', 0): 'R0', ('R', 1): 'R1', ('d', 'A'): 'DA', ('d', 'B'): 'DB',
        ('x', 'A'): 'XA', ('x', 'B'): 'XB', ('e', 0): 'E', ('m', 'A'): 'MA', ('m', 'B'): 'MB',
        ('V', 0): 'V0', ('V', 1): 'V1', ('V', 2): 'V2'}
tests = []
for p_, lam, J, t, mu, K, flux in [(3, .1, 1., 1., 0., 100., np.pi / 6), (14, .37, .3, 1.7, .23, 1e4, np.pi / 28 + .01),
                                   (2, .05, 2., .4, -.5, 37., .41)]:
    pc = bond14.params(p_, D=7., g=.071, td=.003, K=K, musel=.4, field=.17, flux=flux, theta=-np.pi / 4)
    pm = my11.default_params(p_, D=7., Delta=pc['Delta'], g=.071, td=.003, phiA=flux, phiB=0.,
                             psiA=p_ * flux + np.pi / 2, psiB=np.pi / 2, K=K, musel=.4, field=.17)
    for col in (True, False):
        theirs = bond14.bond(pc, lam=lam, J=J, t=t, mu=mu, colored=col)
        mine, ls = my14.bond(pm, lam, J, t, mu, col)
        dd = len(ls)
        tl = bond14.LET if col else bond14.MAT + [('V', 0)]
        perm = [ls.index(name[l]) if (col or l[0] != 'V') else ls.index('V') for l in tl]
        P2 = np.array([perm[i] * dd + perm[j] for i in range(dd) for j in range(dd)])
        mp = mine[np.ix_(P2, P2)]
        q = np.array([0 if l[0] == 'V' else my11.charge(l) for l in ls]); qq = np.add.outer(q, q).ravel()
        no = np.array([l[0] != 'V' for l in ls]).astype(int); nn = np.add.outer(no, no).ravel()
        mk = np.array([l in my11.MARKERS for l in ls]).astype(int); mm = np.add.outer(mk, mk).ravel()
        rec = dict(p=p_, coloured=col, dim=dd * dd, max_entry_difference=float(np.max(abs(mp - theirs))),
                   hermiticity=float(np.max(abs(mine - mine.conj().T))),
                   Q_comm=float(np.max(abs(mine * (qq[:, None] - qq[None, :])))),
                   N_comm=float(np.max(abs(mine * (nn[:, None] - nn[None, :])))),
                   M_comm=float(np.max(abs(mine * (mm[:, None] - mm[None, :])))))
        print(rec, flush=True)
        tests.append(rec)
out['bond_tests'] = tests
M = my11.LETTERS
out['automaton_reproduces_C84_table'] = all(my14.ext_allowed(a, b) == my11.allowed(a, b) for a in M for b in M)
out['bond14_allowed_matches'] = all(
    my14.ext_allowed(name[a], name[b]) == bond14.allowed(a, b) for a in bond14.LET for b in bond14.LET)

# (c) minimality
allowed = {(a, b) for a in M for b in M if my11.allowed(a, b)}


def ring_pairs(ps):
    S = set()
    for p in ps:
        for w in my11.legal_words(p):
            word = w
            for i in range(len(word) - 1):
                S.add((word[i], word[i + 1]))
            S.add((word[-1], word[0]))  # marker to first letter of the next cell (all markers equivalent)
            for m in my11.MARKERS:
                S.add((m, word[0]))
    return S


req = dict(all_allowed=allowed, rings_p1_to_6=ring_pairs(range(1, 7)), rings_fixed_p6=ring_pairs([6]))
rects = set()
for r in range(1, len(M) + 1):
    for A in itertools.combinations(M, r):
        B = frozenset(b for b in M if all((a, b) in allowed for a in A))
        if B:
            A2 = frozenset(a for a in M if all((a, b) in allowed for b in B))
            rects.add((A2, B))
rects = list(rects)
mini = {}
for key, S in req.items():
    best = None
    for k in (1, 2, 3):
        for combo in itertools.combinations(rects, k):
            cov = set()
            for A, B in combo:
                cov |= {(a, b) for a in A for b in B}
            if S <= cov:
                best = k; break
        if best:
            break
    # fooling set of size 3 inside S
    fool = None
    for trip in itertools.combinations(sorted(S), 3):
        ok = all(((x[0], y[1]) not in allowed) or ((y[0], x[1]) not in allowed)
                 for x, y in itertools.combinations(trip, 2))
        if ok:
            fool = trip; break
    mini[key] = dict(required_pairs=len(S), min_rectangles=best, fooling_set_example=fool)
    print(key, mini[key], flush=True)
out['maximal_rectangles'] = len(rects)
out['minimality'] = mini
rows = {}
for a in M:
    rows.setdefault(frozenset(b for b in M if (a, b) in allowed), []).append(a)
out['distinct_rows'] = [sorted(v) for v in rows.values()]
out['deterministic_three_colour_tables_up_to_relabelling'] = 1 if len(rows) == 3 else None

# (d) exhaustive rings
ex = {}
for L in (3, 4, 5):
    viol = 0; nonunique = 0; legal_patterns = 0; zero_pen = 0
    seen = {}
    for w in itertools.product(my14.LET, repeat=L):
        occ = [i for i, l in enumerate(w) if not my14.is_vac(l)]
        if not occ:
            continue
        mat = [w[i] for i in occ]
        nmb = sum(not my11.allowed(mat[i], mat[(i + 1) % len(mat)]) for i in range(len(mat)))
        nb = sum(not my14.ext_allowed(w[i], w[(i + 1) % L]) for i in range(L))
        if nmb > nb:
            viol += 1
        if nb == 0:
            zero_pen += 1
            key = tuple('V' if my14.is_vac(l) else l for l in w)
            seen[key] = seen.get(key, 0) + 1
    # every legal squeezed word with any occupancy pattern should appear exactly once
    for key, c in seen.items():
        if c != 1:
            nonunique += 1
    for pat in itertools.product((0, 1), repeat=L):
        k = sum(pat)
        if k == 0:
            continue
        for mat in itertools.product(M, repeat=k):
            if all(my11.allowed(mat[i], mat[(i + 1) % k]) for i in range(k)):
                legal_patterns += 1
                it = iter(mat)
                key = tuple(next(it) if s else 'V' for s in pat)
                if seen.get(key, 0) != 1:
                    nonunique += 1
    ex[L] = dict(zero_penalty_words=zero_pen, legal_squeezed_patterns=legal_patterns,
                 N_material_bad_gt_N_bad=viol, colouring_not_unique_or_missing=nonunique)
    print(L, ex[L], flush=True)
out['exhaustive_rings'] = ex
json.dump(out, open(os.path.join(HERE, 'check_lift14.json'), 'w'), indent=1)
