"""Item 2: exhaustive check that every literal transition of the referee's h9 / h10
conserves F (adjacent (L,R) pairs of the cyclic word with vacancies removed), the
charge Q, the head number and the occupied number, on all d^L ring words.

Also checks a stronger statement that the README does not spell out:
  B_i = [x_i is not an L] * [x_{i+1} is not an R]   ("cell-boundary bond")
is conserved bond by bond by h9 on every word, so the cell pattern is frozen.
For h10 the analogous statement is that the cyclic pattern of boundary bonds of
the squeezed word is conserved (checked up to rotation).

usage: ring_invariants.py Lmax9 Lmax10
Independent of the author's ring_checks.py (not read before writing this).
"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_bonds as rb

IS_L = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0], bool)
IS_R = np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0], bool)
IS_H = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 0], bool)
IS_V = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1], bool)
Q = np.array([0, 2, 0, 2, 0, -1, 1, -1, 1, 0])


def transitions(H, d, tol=1e-13):
    out = []
    for row in range(d * d):
        for col in range(d * d):
            if row != col and abs(H[row, col]) > tol:
                out.append((col // d, col % d, row // d, row % d))   # (r,s) -> (p,q)
    return out


def all_words(d, L):
    N = d ** L
    idx = np.arange(N, dtype=np.int64)
    w = np.empty((N, L), dtype=np.int8)
    for i in range(L):
        w[:, i] = (idx // d ** (L - 1 - i)) % d
    return w


def invariants(w, d):
    N, L = w.shape
    isL, isR, isH, isV = IS_L[w], IS_R[w], IS_H[w], IS_V[w]
    F = np.zeros(N, dtype=np.int16)
    for i in range(L):
        found = np.zeros(N, bool)
        nxtR = np.zeros(N, bool)
        for k in range(1, L + 1):
            j = (i + k) % L
            new = (~found) & (~isV[:, j])
            nxtR |= new & isR[:, j]
            found |= new
        F += (isL[:, i] & nxtR).astype(np.int16)
    q = Q[w].sum(1).astype(np.int16)
    nh = isH.sum(1).astype(np.int16)
    nocc = (~isV).sum(1).astype(np.int16)
    # boundary-bond bitmask for the unsqueezed word (meaningful for d = 9)
    bmask = np.zeros(N, dtype=np.int32)
    for i in range(L):
        j = (i + 1) % L
        bmask |= ((~isL[:, i]) & (~isR[:, j])).astype(np.int32) << i
    return dict(F=F, Q=q, heads=nh, occ=nocc, bmask=bmask)


def squeezed_boundary_canon(word):
    """canonical (rotation-minimal) boundary pattern of the squeezed cyclic word"""
    s = [x for x in word if x != 9]
    m = len(s)
    if m == 0:
        return ()
    pat = tuple(int((not IS_L[s[i]]) and (not IS_R[s[(i + 1) % m]])) for i in range(m))
    return min(pat[k:] + pat[:k] for k in range(m))


def run(d, L, H, names):
    w = all_words(d, L)
    inv = invariants(w, d)
    tr = transitions(H, d)
    res = dict(d=d, L=L, words=int(d ** L), bond_transitions=len(tr), checked=0,
               viol={k: 0 for k in names})
    pw = [d ** (L - 1 - i) for i in range(L)]
    canon_cache = {}
    for i in range(L):
        j = (i + 1) % L
        if L == 1:
            continue
        for (r, s, p, q) in tr:
            src = np.nonzero((w[:, i] == r) & (w[:, j] == s))[0]
            if len(src) == 0:
                continue
            dst = src + (p - r) * pw[i] + (q - s) * pw[j]
            res["checked"] += int(len(src))
            for k in names:
                if k == "sqz":
                    # python-level check of the squeezed boundary pattern (small L only)
                    bad = 0
                    for a_, b_ in zip(src, dst):
                        ca = canon_cache.get(a_)
                        if ca is None:
                            ca = canon_cache[a_] = squeezed_boundary_canon(w[a_])
                        cb = canon_cache.get(b_)
                        if cb is None:
                            cb = canon_cache[b_] = squeezed_boundary_canon(w[b_])
                        bad += ca != cb
                    res["viol"][k] += int(bad)
                else:
                    res["viol"][k] += int((inv[k][src] != inv[k][dst]).sum())
    return res


if __name__ == "__main__":
    L9 = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    L10 = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = np.random.default_rng(7)
    par = rb.random_par(rng, mobile=True)
    H9, H10 = rb.h9(par), rb.h10(par)
    results = []
    for L in range(2, L9 + 1):
        t0 = time.time()
        r = run(9, L, H9, ["F", "Q", "heads", "occ", "bmask"])
        r["seconds"] = round(time.time() - t0, 1)
        print(json.dumps(r), flush=True)
        results.append(r)
    for L in range(2, L10 + 1):
        t0 = time.time()
        names = ["F", "Q", "heads", "occ"] + (["sqz"] if L <= 5 else [])
        r = run(10, L, H10, names)
        r["seconds"] = round(time.time() - t0, 1)
        print(json.dumps(r), flush=True)
        results.append(r)
    json.dump(results, open(os.path.join(HERE, "ring_invariants.json"), "w"), indent=1)
