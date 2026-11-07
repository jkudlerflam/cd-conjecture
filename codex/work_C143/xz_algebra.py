"""Check the operator identities used at the XZ point on a ring of L = 6:
a_j = X_j Z_{j+1}, b_j = Z_j X_{j+1}: [a_i, b_j] = 0 for all i, j; {a_j, a_{j+1}} = 0 = {b_j, b_{j+1}};
Y_j Y_{j+1} = a_j b_j."""
import numpy as np
from ed import SX, SY, SZ, I2
L = 6
def op(d):
    M = np.array([[1.0 + 0j]])
    for i in range(L): M = np.kron(M, d.get(i, I2))
    return M
a = [op({j: SX, (j + 1) % L: SZ}) for j in range(L)]
b = [op({j: SZ, (j + 1) % L: SX}) for j in range(L)]
comm = max(np.linalg.norm(a[i] @ b[j] - b[j] @ a[i]) for i in range(L) for j in range(L))
anti_a = max(np.linalg.norm(a[j] @ a[(j + 1) % L] + a[(j + 1) % L] @ a[j]) for j in range(L))
anti_b = max(np.linalg.norm(b[j] @ b[(j + 1) % L] + b[(j + 1) % L] @ b[j]) for j in range(L))
far = max(np.linalg.norm(a[i] @ a[j] - a[j] @ a[i]) for i in range(L) for j in range(L) if (i - j) % L not in (1, L - 1))
yy = max(np.linalg.norm(op({j: SY, (j + 1) % L: SY}) - a[j] @ b[j]) for j in range(L))
print("max||[a_i,b_j]||", comm, " max||{a_j,a_j+1}||", anti_a, " max||{b_j,b_j+1}||", anti_b, " max||[a_i,a_j]|| (|i-j|>1)", far, " max||YY - a b||", yy)
