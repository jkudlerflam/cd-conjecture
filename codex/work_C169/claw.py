"""PROOF §8 (C169): which Pauli-aligned B_2 Hamiltonians have a claw-free frustration graph?

Chapman-Flammia: a Pauli Hamiltonian is solvable by a generalized Jordan-Wigner
map iff its frustration graph is a line graph; claw-freeness is necessary for
that and for the 'free fermions in disguise' class (Fendley; Chapman, Elman,
Mann).  A claw is a vertex with three pairwise commuting neighbours.  We test
every subset of the 12 B_2 Pauli families on an open chain of 10 sites,
then reduce modulo the uniform Clifford relabellings of axes (S_3 x signs act
trivially on families) and site reflection (ab <-> ba).
"""
import itertools
from frustration import anti, B2

def graph(fams, L=10):
    verts = []
    for F, s in fams:
        for c in range(L - len(F) + 1):
            verts.append((F, c))
    adj = {v: set() for v in range(len(verts))}
    for u in range(len(verts)):
        for v in range(u + 1, len(verts)):
            if anti(verts[u][0], verts[u][1], verts[v][0], verts[v][1]):
                adj[u].add(v); adj[v].add(u)
    return verts, adj

def has_claw(adj):
    for v, nb in adj.items():
        nb = sorted(nb)
        for a, b, c in itertools.combinations(nb, 3):
            if b not in adj[a] and c not in adj[a] and c not in adj[b]:
                return True
    return False

def canon(S):
    perms = list(itertools.permutations('XYZ'))
    best = None
    for P in perms:
        m = dict(zip('XYZ', P))
        for refl in (False, True):
            T = []
            for F in S:
                G = ''.join(m[ch] for ch in F)
                if refl:
                    G = G[::-1]
                T.append(G)
            key = tuple(sorted(T))
            if best is None or key < best:
                best = key
    return best

if __name__ == '__main__':
    found = {}
    for r in range(1, 13):
        for S in itertools.combinations([f for f, _ in B2], r):
            verts, adj = graph([(f, 0) for f in S])
            if not any(adj.values()):
                continue        # commuting (classical) Hamiltonian
            if not has_claw(adj):
                found.setdefault(canon(S), S)
    print(len(found), 'claw-free non-commuting classes (up to axis permutation and reflection):')
    for k in sorted(found, key=lambda t: (len(t), t)):
        print('  ', ' + '.join(k))
