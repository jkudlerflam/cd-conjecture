"""Exact free-fermion open-chain entanglement for the estimator calibration (same S-bar and c_pair as summ.py):
XX chain (c = 1), critical TFI (c = 1/2), and the decoupled ladders XX x XX (c = 2) and TFI x XX (c = 3/2),
whose entropies are exact sums since the ladder ground state is a product over the two legs."""
import numpy as np, math
def S_xx(L):
    # H = sum (XX+YY) = 2 sum (c+c + h.c.): hopping matrix, half filling ground state
    T = np.zeros((L, L)); 
    for j in range(L - 1): T[j, j + 1] = T[j + 1, j] = 2.0
    w, v = np.linalg.eigh(T); occ = v[:, w < 0]; C = occ @ occ.T
    S = []
    for l in range(1, L):
        e = np.linalg.eigvalsh(C[:l, :l]); e = e[(e > 1e-14) & (e < 1 - 1e-14)]
        S.append(float(-np.sum(e * np.log(e) + (1 - e) * np.log(1 - e))))
    return np.array(S)
def S_tfi(L):
    # H = -sum X_j X_{j+1} - sum Z_j, open; Majorana form H = (i/4) sum A_ab g_a g_b
    n = 2 * L; A = np.zeros((n, n))
    for j in range(L):           # -Z_j = i g_{2j} g_{2j+1} (convention), coefficient h = 1
        A[2 * j, 2 * j + 1] += 2.0; A[2 * j + 1, 2 * j] -= 2.0
    for j in range(L - 1):       # -X_j X_{j+1} = i g_{2j+1} g_{2j+2}
        A[2 * j + 1, 2 * j + 2] += 2.0; A[2 * j + 2, 2 * j + 1] -= 2.0
    # ground-state Majorana correlation G = <i g_a g_b> = -sign(iA) structure: G = i * sign of the Hermitian iA
    w, v = np.linalg.eigh(1j * A)
    G = (v @ np.diag(-np.sign(w)) @ v.conj().T)   # Hermitian, eigenvalues +-1
    S = []
    for l in range(1, L):
        nu = np.linalg.eigvalsh(G[:2 * l, :2 * l]); nu = np.clip(np.abs(nu[nu > 0]), 0, 1 - 1e-15)
        p = (1 + nu) / 2; q = 1 - p
        S.append(float(-np.sum(p * np.log(p + 1e-300) + q * np.log(q + 1e-300))))
    return np.array(S)
def sbar(S, L): m = L // 2; return float(np.mean(S[m - 3:m + 1]))
if __name__ == '__main__':
    sizes = [48, 96, 192, 384]
    data = {L: (S_xx(L), S_tfi(L)) for L in sizes}
    for name, f in [('XX c=1', lambda a, b: a), ('TFI c=1/2', lambda a, b: b), ('XXxXX c=2', lambda a, b: 2 * a), ('TFIxXX c=3/2', lambda a, b: a + b)]:
        sb = {L: sbar(f(*data[L]), L) for L in sizes}
        print(name, '  '.join(f"c({a},{b})={6*(sb[b]-sb[a])/math.log(b/a):.3f}" for a, b in zip(sizes, sizes[1:])))
