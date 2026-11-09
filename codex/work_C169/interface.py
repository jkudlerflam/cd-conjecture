"""PROOF §3 (C169): the qubit-interface test for blocking.

If a chain on (C^{2^b})^{⊗L} is a blocked B_2 chain up to block unitaries
U(2^b)^{⊗L}, the part of the Hamiltonian acting nontrivially on both sides of
a block boundary is sigma_last^T J sigma_first, so its left operator-Schmidt
space V_L (a canonical, basis-free subspace) satisfies
  (i) dim V_L <= 3, and (ii) V_L is unitarily a subspace of span{sigma^a ⊗ 1}:
      every v in V_L squares to a multiple of 1 (equivalently the Hermitian
      orthonormal basis obeys v_a v_b + v_b v_a = 2 delta_ab 1/d-normalized),
      and if dim V_L = 3, the triple product of an orthonormal basis is +-i 1.
Both are invariant under block unitaries, so failing either excludes the model.
"""
import itertools, numpy as np

P = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]], complex),
     'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1., -1]).astype(complex)}

def kron(s):
    M = np.array([[1.]], complex)
    for ch in s:
        M = np.kron(M, P[ch])
    return M

def cross_terms(fams, b):
    """fams: list of (coef, string) for a one-site-translation-invariant qubit
    chain.  Returns dict {(Lstring, Rstring): coef} of terms crossing the cut
    between sites b-1 and b (block 0 = sites 0..b-1, block 1 = b..2b-1)."""
    M = {}
    for c, s in fams:
        r = len(s)
        for start in range(-r + 1, 2 * b):
            sites = [(start + i, ch) for i, ch in enumerate(s) if ch != 'I']
            if not sites:
                continue
            lo = min(x for x, _ in sites); hi = max(x for x, _ in sites)
            if lo < 0 or hi >= 2 * b:
                continue  # extends beyond the two blocks: only possible if range > b
            if not (lo < b <= hi):
                continue
            Ls = ['I'] * b; Rs = ['I'] * b
            for x, ch in sites:
                if x < b: Ls[x] = ch
                else: Rs[x - b] = ch
            key = (''.join(Ls), ''.join(Rs))
            M[key] = M.get(key, 0) + c
    return M

def test(fams, b, name):
    M = cross_terms(fams, b)
    Ls = sorted({k[0] for k in M}); Rs = sorted({k[1] for k in M})
    A = np.zeros((len(Ls), len(Rs)))
    for (l, r), c in M.items():
        A[Ls.index(l), Rs.index(r)] += c
    U, s, Vt = np.linalg.svd(A)
    k = int(np.sum(s > 1e-10))
    d = 2 ** b
    V = [sum(U[i, a] * kron(Ls[i]) for i in range(len(Ls))) for a in range(k)]
    cliff = max(np.abs(V[a] @ V[c] + V[c] @ V[a] - 2 * (a == c) * np.eye(d)).max()
                for a in range(k) for c in range(k))
    triple = None
    if k == 3:
        T = V[0] @ V[1] @ V[2]
        triple = np.abs(T - T[0, 0] * np.eye(d)).max()
    ok = (k <= 3) and cliff < 1e-9 and (triple is None or triple < 1e-9)
    print(f'{name:38s} b={b}: dim V_L = {k}, Clifford defect = {cliff:.3f}, '
          f'triple defect = {triple if triple is None else round(triple, 3)}  -> '
          f'{"passes" if ok else "fails"}')
    return ok

MODELS = {
    'B_2 XXZ + field (control)': [(1, 'XX'), (1, 'YY'), (0.7, 'ZZ'), (0.3, 'Z')],
    'B_2 generic bond (control)': [(c, a + b) for c, (a, b) in zip(np.random.default_rng(0).normal(size=9), itertools.product('XYZ', repeat=2))] + [(0.4, 'X'), (0.2, 'Z')],
    'cluster-Ising (c=3/2)': [(1, 'XZX'), (1, 'YY')],
    'RZFA spin form (c=3/2)': [(-1, 'Z'), (-1, 'XX'), (-0.5, 'ZZ'), (-0.5, 'XIX')],
    'Titvinidze-Japaridze (c=2)': [(1, 'XX'), (1, 'YY'), (0.8, 'XZX'), (0.8, 'YZY')],
    'J1-J2 XXZ (c=2 in field)': [(1, 'XX'), (1, 'YY'), (1, 'ZZ'), (0.6, 'XIX'), (0.6, 'YIY'), (0.6, 'ZIZ')],
    'two XX legs, zigzag (c=2)': [(1, 'XIX'), (1, 'YIY')],
    "O'Brien-Fendley": [(1, 'XX'), (1, 'Z'), (0.4, 'XXZ'), (0.4, 'ZXX')],
}

def ladder_test(name, cross):
    """d = 4 ladders written directly as cross terms between (A,B) blocks."""
    fams = []
    print_name = name
    M = {}
    for c, l, r in cross:
        M[(l, r)] = M.get((l, r), 0) + c
    Ls = sorted({k[0] for k in M}); Rs = sorted({k[1] for k in M})
    A = np.zeros((len(Ls), len(Rs)))
    for (l, r), c in M.items():
        A[Ls.index(l), Rs.index(r)] += c
    U, s, Vt = np.linalg.svd(A); k = int(np.sum(s > 1e-10))
    V = [sum(U[i, a] * kron(Ls[i]) for i in range(len(Ls))) for a in range(k)]
    cliff = max(np.abs(V[a] @ V[c] + V[c] @ V[a] - 2 * (a == c) * np.eye(4)).max() for a in range(k) for c in range(k))
    print(f'{print_name:38s} d=4: dim V_L = {k}, Clifford defect = {cliff:.3f} -> {"passes" if k <= 3 and cliff < 1e-9 else "fails"}')

if __name__ == '__main__':
    for name, f in MODELS.items():
        for b in (2, 3):
            test(f, b, name)
    # d = 4 ladders: legs A (first qubit) and B (second qubit) of each rung
    ladder_test('Ising x XX ladder (c=3/2)', [(1, 'XI', 'XI'), (1, 'IX', 'IX'), (1, 'IY', 'IY')])
    ladder_test('two XX ladder (c=2)', [(1, 'XI', 'XI'), (1, 'YI', 'YI'), (1, 'IX', 'IX'), (1, 'IY', 'IY')])
    ladder_test('Ising x Ising ladder (c=1, not B_2)', [(1, 'XI', 'XI'), (1, 'IX', 'IX')])
