"""Frustration-graph invariants of translation-invariant Pauli Hamiltonians (C169).

A family is a Pauli string given as a string over {I,X,Y,Z}, placed at sites
j, j+1, ... and translated by the period p.  Two Pauli strings anticommute iff
the number of sites where both act nontrivially with different letters is odd.

Invariants under Clifford maps (on-site Cliffords, cluster entanglers, KW and
JW restricted to the symmetric algebra, any Clifford QCA):
  D(F,G) = { k : F_0 anticommutes with G_k }   (preserved exactly by a map that
           commutes with translation, up to a common shift of G's origin);
  the number of infinite connected components of the frustration graph.
"""
import itertools, sys

def anti(a, sa, b, sb):
    """a, b strings; sa, sb start sites. True if anticommute."""
    n = 0
    for i, ca in enumerate(a):
        if ca == 'I':
            continue
        jb = sa + i - sb
        if 0 <= jb < len(b):
            cb = b[jb]
            if cb != 'I' and cb != ca:
                n += 1
    return n % 2 == 1

def offsets(F, G, p=1, R=8):
    return sorted(k for k in range(-R, R + 1) if anti(F, 0, G, k * p))

def diam(D):
    return (max(D) - min(D)) if D else -1

def components(fams, L, p=1, ring=False):
    """Frustration graph on an open chain (ring=False) or ring of L cells (period p).
    Open chains give the infinite-chain component count without wrap-around merging.
    Returns list of component sizes (vertices with at least one edge)."""
    N = L * p
    verts = []
    for fi, (F, s) in enumerate(fams):
        for c in range(L):
            if not ring and c * p + s + len(F) > N:
                continue
            verts.append((fi, (c * p + s) % N))
    def ring_anti(F, sa, G, sb):
        # place on ring of N sites
        A = {}
        for i, ch in enumerate(F):
            if ch != 'I':
                A[(sa + i) % N if ring else sa + i] = ch
        n = 0
        for i, ch in enumerate(G):
            if ch != 'I':
                q = (sb + i) % N if ring else sb + i
                if q in A and A[q] != ch:
                    n += 1
        return n % 2 == 1
    nv = len(verts)
    adj = [[] for _ in range(nv)]
    for u in range(nv):
        for v in range(u + 1, nv):
            Fu, su = fams[verts[u][0]][0], verts[u][1]
            Fv, sv = fams[verts[v][0]][0], verts[v][1]
            if ring_anti(Fu, su, Fv, sv):
                adj[u].append(v); adj[v].append(u)
    seen = [False] * nv; sizes = []
    for u in range(nv):
        if seen[u] or not adj[u]:
            continue
        st = [u]; seen[u] = True; n = 0
        while st:
            x = st.pop(); n += 1
            for y in adj[x]:
                if not seen[y]:
                    seen[y] = True; st.append(y)
        sizes.append(n)
    return sizes

MODELS = {
    # one-site-cell models (p = 1): list of (string, start offset)
    'cluster-Ising (Smacchia et al., c=3/2)': [('XZX', 0), ('YY', 0)],
    'RZFA spin form (c=3/2)': [('Z', 0), ('XX', 0), ('ZZ', 0), ('XIX', 0)],
    "O'Brien-Fendley (c=7/10 TCI point)": [('XX', 0), ('Z', 0), ('XXZ', 0), ('ZXX', 0)],
    'Titvinidze-Japaridze (c=2)': [('XX', 0), ('YY', 0), ('XZX', 0), ('YZY', 0)],
    'J1-J2 XXZ (c=2 in field)': [('XX', 0), ('YY', 0), ('ZZ', 0), ('XIX', 0), ('YIY', 0), ('ZIZ', 0), ('Z', 0)],
    'two XX legs, zigzag (c=2)': [('XIX', 0), ('YIY', 0)],
}
MODELS_P2 = {
    # two-site cell (p = 2)
    'Ising x XX, zigzag (c=3/2)': [('XIX', 0), ('Z', 0), ('XIX', 1), ('YIY', 1)],
    'Ising x XX, d=4 rungless ladder (c=3/2)': [('XIX', 0), ('Z', 0), ('XIX', 1), ('YIY', 1)],
}

B2 = [(a, 0) for a in 'XYZ'] + [(a + b, 0) for a in 'XYZ' for b in 'XYZ']

def model_report(name, fams, p=1):
    worst = (-1, None)
    for (F, s), (G, t) in itertools.product(fams, repeat=2):
        D = [k for k in range(-8, 9) if anti(F, s, G, t + k * p)]
        d = diam(D)
        if d > worst[0]:
            worst = (d, (F, s, G, t, D))
    comp = components(fams, 30 // p * 1, p)
    big = [c for c in comp if c >= 8]
    return worst, len(big), comp

if __name__ == '__main__':
    # B_2: every pair of families
    wd = -1
    for (F, s), (G, t) in itertools.product(B2, repeat=2):
        d = diam(offsets(F, G))
        wd = max(wd, d)
    print('B_2: max offset diameter over all 144 family pairs =', wd)
    # B_2: max number of infinite components over all subsets of the 12 families
    best = (0, None)
    for r in range(1, 13):
        for S in itertools.combinations(B2, r):
            comp = components(list(S), 20)
            nbig = sum(1 for c in comp if c >= 10)
            if nbig > best[0]:
                best = (nbig, [x[0] for x in S])
    print('B_2: max number of infinite frustration components over all 4095 subsets =', best[0], 'e.g.', best[1])
    print()
    for name, fams in MODELS.items():
        (d, ex), nb, comp = model_report(name, fams, 1)
        print(f'{name}: max diam {d} (pair {ex[0]}@{ex[1]} vs {ex[2]}@{ex[3]}, D={ex[4]}), infinite components {nb}')
    for name, fams in MODELS_P2.items():
        (d, ex), nb, comp = model_report(name, fams, 2)
        print(f'{name} [p=2]: max diam in cells {d}, infinite components {nb}')
