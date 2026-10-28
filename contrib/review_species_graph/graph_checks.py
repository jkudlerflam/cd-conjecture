"""Independent re-implementation of the girth-six selection rule.

Written from the statement in contrib/torus_bottleneck/README.md, without
reusing its checks.py.  For the point-line incidence graph of PG(2,p) it
(i) builds the graph from F_p^3 directly, (ii) checks degrees, unique
incidences, girth and connectivity, (iii) enumerates EVERY multiset of at
most `kmax` fermion fields (flavour, creation/annihilation) and lists the
species-neutral ones, and (iv) checks that each neutral monomial with at
most four fields conserves every flavour.  For p=2 it also enumerates all
six-field monomials and counts the flavour-changing neutral ones, which
must be exactly the two orientations of each hexagon.

numpy only.  Run: python3 graph_checks.py   (about a minute)
"""
import itertools, json, sys, time
from pathlib import Path
import numpy as np


def projective_plane(p):
    pts = []
    for v in itertools.product(range(p), repeat=3):
        if any(v):
            first = next(x for x in v if x)
            inv = pow(first, p - 2, p)
            w = tuple((x * inv) % p for x in v)
            if w not in pts:
                pts.append(w)
    lines = list(pts)  # covectors, same normalisation
    edges = [(i, j) for i, P in enumerate(lines) for j, X in enumerate(pts)
             if sum(a * b for a, b in zip(P, X)) % p == 0]
    return len(lines), len(pts), edges  # upper = lines, lower = points


def incidence_matrix(u, v, edges):
    W = np.zeros((u + v, len(edges)), dtype=np.int64)
    for a, (l, r) in enumerate(edges):
        W[l, a] = 1          # e_(upper, l)
        W[u + r, a] = -1     # -e_(lower, r)
    return W


def girth_bipartite(u, v, edges):
    n = u + v
    adj = [[] for _ in range(n)]
    for l, r in edges:
        adj[l].append(u + r); adj[u + r].append(l)
    best = 10 ** 9
    for s in range(n):
        dist = [-1] * n; par = [-1] * n; dist[s] = 0; q = [s]
        for x in q:
            for y in adj[x]:
                if dist[y] < 0:
                    dist[y] = dist[x] + 1; par[y] = x; q.append(y)
                elif par[x] != y:
                    best = min(best, dist[x] + dist[y] + 1)
    connected = all(d >= 0 for d in dist)
    return best, connected


def neutral_monomials(W, kmax, chunk=400000):
    """Enumerate all multisets of <= kmax fields.  A field is (flavour a, s),
    s=+1 creation, -1 annihilation; its species charge is s*W[:,a].
    Returns counts and the list of neutral flavour-changing ones."""
    q, N = W.shape
    fields = [(a, s) for a in range(N) for s in (+1, -1)]
    F = np.array([s * W[:, a] for a, s in fields], dtype=np.int64)       # (2N, q)
    Z = np.array([s * np.eye(N, dtype=np.int64)[a] for a, s in fields])  # (2N, N)
    out = {}
    bad = []
    for k in range(1, kmax + 1):
        total = neutral = changing = 0
        it = itertools.combinations_with_replacement(range(len(fields)), k)
        while True:
            block = list(itertools.islice(it, chunk))
            if not block:
                break
            idx = np.array(block, dtype=np.int64)
            charge = F[idx].sum(axis=1)
            mask = ~charge.any(axis=1)
            total += len(idx); neutral += int(mask.sum())
            if mask.any():
                z = Z[idx[mask]].sum(axis=1)
                chg = z.any(axis=1)
                changing += int(chg.sum())
                if chg.any() and k <= 4:
                    bad.extend(idx[mask][chg].tolist())
        out[k] = dict(multisets=total, neutral=neutral, neutral_flavour_changing=changing)
    return out, bad


def kernel_rank(W):
    return W.shape[1] - np.linalg.matrix_rank(W.astype(float))


def main():
    res = {}
    for p, kmax in [(2, 6), (3, 4)]:
        t0 = time.time()
        u, v, edges = projective_plane(p)
        N = len(edges); q = u + v
        W = incidence_matrix(u, v, edges)
        deg_u = np.bincount([l for l, r in edges], minlength=u)
        deg_v = np.bincount([r for l, r in edges], minlength=v)
        # any two points lie on exactly one common line, and dually
        A = np.zeros((u, v), int)
        for l, r in edges:
            A[l, r] = 1
        pp = A.T @ A; ll = A @ A.T
        g, conn = girth_bipartite(u, v, edges)
        mono, bad = neutral_monomials(W, kmax)
        s = p * p + p + 1
        triangles = s * (s - 1) * (s - p - 1) // 6
        r = dict(p=p, u=u, v=v, q=q, N=N, expected_N=(p + 1) * s,
                 degrees=sorted(set(deg_u.tolist()) | set(deg_v.tolist())),
                 two_points_one_line=bool((pp[~np.eye(v, dtype=bool)] == 1).all()),
                 two_lines_one_point=bool((ll[~np.eye(u, dtype=bool)] == 1).all()),
                 girth=g, connected=conn, rank_W=int(np.linalg.matrix_rank(W.astype(float))),
                 kernel_rank=int(kernel_rank(W)), expected_kernel_rank=N - q + 1,
                 monomials=mono,
                 flavour_changing_neutral_with_at_most_4_fields=len(bad),
                 hexagons=triangles,
                 seconds=round(time.time() - t0, 1))
        if kmax >= 6:
            r['six_field_conversions_expected'] = 2 * triangles
        res[f'p={p}'] = r
        print(json.dumps(r, indent=1)); sys.stdout.flush()
        assert r['N'] == r['expected_N'] and r['girth'] == 6 and r['connected']
        assert r['two_points_one_line'] and r['two_lines_one_point']
        assert r['kernel_rank'] == r['expected_kernel_rank']
        assert len(bad) == 0
        if kmax >= 6:
            assert mono[6]['neutral_flavour_changing'] == 2 * triangles
            assert mono[5]['neutral'] == 0  # odd field numbers are never neutral
    Path(__file__).with_name('graph_results.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
