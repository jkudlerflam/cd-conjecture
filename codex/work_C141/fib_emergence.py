"""Emergence at the d = 2 Fibonacci point (PROOF.md section 5.3).
h* = n (x) n has an extensively degenerate ground space (no two adjacent up spins),
gapped by Delta* = 1 from the rest.  Perturb by delta (s+ s- + s- s+).
Prediction of the emergence lemma: the emergent sector has v = delta * v_eff + O(delta^2),
with v_eff the velocity of the first-order effective (projected) model.
We measure v from the lowest momentum-2pi/L excitation in the ground-state particle sector.
"""
import numpy as np, itertools
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

def basis(L, N, constrained):
    out = []
    for occ in itertools.combinations(range(L), N):
        s = 0
        for i in occ: s |= 1 << i
        if constrained:
            if s & ((s << 1) | (s >> (L-1))) & ((1 << L)-1): continue
        out.append(s)
    return out

def ham(L, N, delta, constrained):
    B = basis(L, N, constrained); idx = {s: i for i, s in enumerate(B)}
    rows, cols, vals = [], [], []
    for a, s in enumerate(B):
        diag = 0.0
        for i in range(L):
            j = (i+1) % L
            ni, nj = (s >> i) & 1, (s >> j) & 1
            diag += ni*nj
            if ni != nj:
                t = s ^ ((1 << i) | (1 << j))
                if t in idx:
                    rows.append(idx[t]); cols.append(a); vals.append(delta)
        if not constrained:
            rows.append(a); cols.append(a); vals.append(diag)
    n = len(B)
    return csr_matrix((vals, (rows, cols)), shape=(n, n)), B

def translate(B, L):
    idx = {s: i for i, s in enumerate(B)}
    perm = np.empty(len(B), dtype=int)
    for a, s in enumerate(B):
        t = ((s << 1) | (s >> (L-1))) & ((1 << L)-1)
        perm[a] = idx[t]
    return perm

def spectrum_with_k(L, N, delta, constrained, nev=8):
    H, B = ham(L, N, delta, constrained)
    if H.shape[0] <= 600:
        w, V = np.linalg.eigh(H.toarray()); w, V = w[:nev], V[:, :nev]
    else:
        w, V = eigsh(H, k=nev, which='SA')
    o = np.argsort(w); w, V = w[o], V[:, o]
    perm = translate(B, L)
    TV = np.zeros_like(V); TV[perm, :] = V
    out = []; i = 0
    while i < len(w):
        j = i
        while j+1 < len(w) and abs(w[j+1]-w[i]) < 1e-9*max(1, abs(w[i])): j += 1
        sub = V[:, i:j+1]; Tm = sub.conj().T @ TV[:, i:j+1]
        ev = np.linalg.eigvals(Tm)
        for e in ev:
            out.append((w[i], np.angle(e)))
        i = j+1
    return out

for L in [12, 14, 16, 18, 20]:
    # first-order model: find best N
    best = None
    for N in range(1, L//2+1):
        E = spectrum_with_k(L, N, 1.0, True, nev=4)[0][0]
        if best is None or E < best[0]: best = (E, N)
    N = best[1]
    sp = spectrum_with_k(L, N, 1.0, True, nev=24)
    E0 = sp[0][0]
    k1 = 2*np.pi/L
    exl = [e for e, kk in sp if abs(abs(kk)-k1) < 1e-6]
    if not exl: print(L, N, "no k1 state", sp[:6]); continue
    ex = min(exl)
    veff = (ex-E0)*L/(2*np.pi)
    line = f"L={L} N*={N} projected v_eff={veff:.5f} |"
    for delta in [0.01, 0.03, 0.1]:
        spf = spectrum_with_k(L, N, delta, False, nev=24)
        E0f = spf[0][0]
        exf = [e for e, kk in spf if abs(abs(kk)-k1) < 1e-6]
        vf = (min(exf)-E0f)*L/(2*np.pi) if exf else float('nan')
        line += f" d={delta}: v/d={vf/delta:.5f}"
    print(line)
