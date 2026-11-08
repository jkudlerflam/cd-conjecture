"""C158: check of the A/B algebra at XZ on an L = 8 ring (Pauli strings as symplectic vectors).
a_j = X_j Z_{j+1}, b_j = Z_j X_{j+1}.  Two Pauli strings commute iff their symplectic product is 0."""
import numpy as np
L = 8
def P(ops):
    x = np.zeros(L, int); z = np.zeros(L, int)
    for site, o in ops:
        s = site % L
        if o in 'XY': x[s] ^= 1
        if o in 'ZY': z[s] ^= 1
    return x, z
def anti(p, q):
    return int((p[0] @ q[1] + p[1] @ q[0]) % 2)
a = [P([(j, 'X'), (j + 1, 'Z')]) for j in range(L)]
b = [P([(j, 'Z'), (j + 1, 'X')]) for j in range(L)]
print('max anticommutator [a_j, b_k]:', max(anti(a[j], b[k]) for j in range(L) for k in range(L)))
print('a_j anticommutes with a_k only for |j-k| = 1:', all(anti(a[j], a[k]) == (min((j - k) % L, (k - j) % L) == 1) for j in range(L) for k in range(L)))
print('b likewise:', all(anti(b[j], b[k]) == (min((j - k) % L, (k - j) % L) == 1) for j in range(L) for k in range(L)))
j = 3
for name, op in (('X', P([(j, 'X')])), ('Z', P([(j, 'Z')])), ('Y', P([(j, 'Y')]))):
    print(f'{name}_{j} anticommutes with a at', [k for k in range(L) if anti(op, a[k])], 'and b at', [k for k in range(L) if anti(op, b[k])])
