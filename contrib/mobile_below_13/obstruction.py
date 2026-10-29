"""Finite enumerations behind the obstruction statements (exact, combinatorial).

Material relation: the pinless eleven-letter grammar (C88 automaton), legal(a,b)
iff o(a) in I(b).
 (1) USED: ordered material pairs that occur in legal frozen rings (inside cells
     of every length p<=6, all clock/token/marker states and bit strings, and
     across cell boundaries).
 (2) every maximal rectangle A x B inside the legal relation (A = letters allowed
     left of an uncoloured vacancy, B = letters allowed right of it).
 (3) minimal number of rectangles inside LEGAL covering USED (= minimal number of
     vacancy colours if every bond of every legal word must admit a vacancy,
     the material-independent placement used by C82/C83/C88).
 (4) for each maximal rectangle, which zero-port cell states a single vacancy can
     cross with legality-preserving swaps, and whether crossing is bit blind.
"""
import os, json, itertools
import bond12 as b

HERE = os.path.dirname(os.path.abspath(__file__))
M = b.MAT
leg = {(x, y): b.allowed(x, y) for x in M for y in M}
L = lambda a: ('L', a)
R = lambda a: ('R', a)
TOK = [('R', 0), ('d', 'A'), ('d', 'B'), ('x', 'A'), ('x', 'B')]
MK = [('e', 0), ('m', 'A'), ('m', 'B')]


def cells(p, include_pinned=False):
    out = []
    for bits in itertools.product((0, 1), repeat=p - 1):
        for j in range(p):
            for T in TOK:
                out.append(tuple([L(a) for a in bits[:j]] + [T] + [R(a) for a in bits[j:]] + [('e', 0)]))
        for m in MK:
            out.append(tuple([L(a) for a in bits] + [L(0), m]))
            if include_pinned:
                out.append(tuple([L(a) for a in bits] + [L(1), m]))
    return out


used = set()
firsts = set(); lasts = set()
for p in range(1, 7):
    for c in cells(p):
        assert all(leg[(c[i], c[i + 1])] for i in range(len(c) - 1))
        used |= {(c[i], c[i + 1]) for i in range(len(c) - 1)}
        firsts.add(c[0]); lasts.add(c[-1])
used |= {(x, y) for x in lasts for y in firsts}
assert all(leg[u] for u in used)
res = dict(used_pairs=len(used), legal_pairs=sum(leg.values()))

# (2) maximal rectangles
rects = []
for r in range(1, 12):
    for A in itertools.combinations(M, r):
        B = tuple(y for y in M if all(leg[(x, y)] for x in A))
        if not B:
            continue
        A2 = tuple(x for x in M if all(leg[(x, y)] for y in B))
        if A2 == A and (A, B) not in rects:
            rects.append((A, B))
res['maximal_rectangles'] = [dict(A=[str(x) for x in A], B=[str(y) for y in B]) for A, B in rects]

# (3) minimal cover of USED
best = None
for k in range(1, len(rects) + 1):
    for combo in itertools.combinations(range(len(rects)), k):
        cov = set()
        for i in combo:
            A, B = rects[i]; cov |= {(x, y) for x in A for y in B}
        if used <= cov:
            best = (k, combo); break
    if best:
        break
res['min_colours_material_independent'] = best[0]
res['min_cover'] = [dict(A=[str(x) for x in rects[i][0]], B=[str(y) for y in rects[i][1]]) for i in best[1]]
res['single_rectangle_covering_used'] = any(used <= {(x, y) for x in A for y in B} for A, B in rects)


# (4) single-vacancy crossing analysis: vacancy enters from the left neighbour
# marker E, must reach the bond after the cell's E, with swaps V a <-> a V that
# keep every bond legal in the extended grammar (V allowed after A, before B).
def crossing(A, B, cell):
    """Can a vacancy placed after the previous E cross the frozen cell word
    (letters held fixed) by swaps?  Returns True iff every bond of
    E + cell + next-first-letter admits V and every crossed letter keeps legality."""
    word = (('e', 0),) + cell + (('R', 0),)   # next cell's walker at position 0
    for i in range(len(word) - 1):
        if not (word[i] in A and word[i + 1] in B):
            return False
    return True


cross = []
for A, B in rects:
    per = {}
    for p in (3, 4):
        for c in cells(p):
            if any(x[0] in 'dxm' for x in c):
                continue  # zero port: clock words only
            j = next((i for i, x in enumerate(c) if x == ('R', 0) and all(y[0] == 'L' for y in c[:i])), None)
            key = 'end' if j is None else 'walker_at_%d' % j
            bits = tuple(x[1] for x in c if x[0] in 'LR' and x != c[j if j is not None else -1])
            per.setdefault((p, key), []).append(crossing(A, B, c))
    summary = {}
    for (p, key), v in per.items():
        summary['p%d_%s' % (p, key)] = 'all' if all(v) else ('none' if not any(v) else 'bit-dependent')
    cross.append(dict(A=[str(x) for x in A], B=[str(y) for y in B], crossing=summary))
res['single_vacancy_crossing'] = cross
print(json.dumps({k: v for k, v in res.items() if k != 'single_vacancy_crossing'}, indent=1))
for c in cross:
    print(c['A'], '|', c['B'])
    print('   ', {k: v for k, v in c['crossing'].items() if k.startswith('p4')})
json.dump(res, open(os.path.join(HERE, 'obstruction.json'), 'w'), indent=1)
