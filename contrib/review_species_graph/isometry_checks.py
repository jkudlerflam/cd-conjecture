"""Independent re-implementation of the sparse two-port isometry, its
compressions, its species covariance, the charge-preserving word transducer,
the three-gate controlled swap, and the exact uniform-amplitude preparation.

Written from contrib/exact_route/{README,CIRCUIT,EXACT_AMPLITUDES}.md only;
none of the repository's verify scripts is imported.

Graphs: the 6-cycle (smallest girth-six bipartite graph), a graph with a
4-cycle and an isolated upper vertex, and (channel identities only, with
sparse Jordan-Wigner operators on 2^21 states) the Fano incidence graph.

numpy + scipy.  Run: python3 isometry_checks.py   (about a minute)
"""
import itertools, json, math
from pathlib import Path
import numpy as np
import scipy.sparse as sp

from graph_checks import projective_plane

RES = {}


# --------------------------------------------------------------- Fock space
def jw_ops(N, sparse=False):
    """f_e (annihilation) with Jordan-Wigner order e=0..N-1, occupation basis
    index sum_e n_e 2^(N-1-e).  Returns list of annihilators and parity."""
    dim = 2 ** N
    idx = np.arange(dim)
    occ = (idx[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1
    P = (-1.0) ** occ.sum(1)
    fs = []
    for e in range(N):
        src = idx[occ[:, e] == 1]
        dst = src - (1 << (N - 1 - e))
        sign = (-1.0) ** occ[src, :e].sum(1)
        M = sp.csr_matrix((sign, (dst, src)), shape=(dim, dim))
        fs.append(M if sparse else M.toarray())
    Pm = sp.diags(P).tocsr() if sparse else np.diag(P)
    return fs, Pm, occ


# --------------------------------------------------------------- one port map
def one_port(u, v, edges, fs, P=None, a=None, b=None):
    """V: F_N -> B (x) A (x) F_N.  B basis: (1,alpha) alpha<u then (0,beta).
    A basis: g_alpha (alpha<u) then h_e (e<N).  If P is given the stage uses
    g_e = P f_e (right stage).  Returns dense V of shape (q*(u+N)*D, D)."""
    N = len(edges); q = u + v; nA = u + N; D = fs[0].shape[0]
    a = a if a is not None else 1 / math.sqrt(2 * u)
    b = b if b is not None else 1 / math.sqrt(2 * N)
    ann = [P @ f for f in fs] if P is not None else fs   # g_e
    V = np.zeros((q, nA, D, D))
    for al in range(u):
        V[al, al] = a * np.eye(D)
    for e, (al, be) in enumerate(edges):
        V[u + be, al] += b * ann[e].T.conj()   # |g_alpha> g_e^dagger
        V[u + be, u + e] += b * ann[e]         # |h_e> g_e
    return V.reshape(q * nA * D, D)


def two_port(u, v, edges):
    N = len(edges); q = u + v; nA = u + N
    fs, P, occ = jw_ops(N)
    D = 2 ** N
    VL = one_port(u, v, edges, fs)
    VR = one_port(u, v, edges, fs, P=P)
    # W = (I_{B_L A_L} (x) V_R) V_L ; output order B_L, A_L, B_R, A_R, F
    W = (VR @ VL.reshape(q * nA, D, D).transpose(1, 0, 2).reshape(D, -1)).reshape(q * nA * D, q * nA, D)
    W = W.transpose(1, 0, 2).reshape(q * nA * q * nA * D, D)
    return W, VL, VR, fs, P, occ


def port_op(q, nA, D, which, X):
    """X acting on B_L (which='L') or B_R (which='R') of the two-port output."""
    if which == 'L':
        return lambda W: np.einsum('pk,kxy->pxy', X, W.reshape(q, -1, D)).reshape(W.shape)
    return lambda W: np.einsum('pk,akxy->apxy', X, W.reshape(q * nA, q, -1, D)).reshape(W.shape)


def check_graph(name, u, v, edges):
    N = len(edges); q = u + v; nA = u + N; D = 2 ** N
    W, VL, VR, fs, P, occ = two_port(u, v, edges)
    a2, b2 = 1 / (2 * u), 1 / (2 * N)
    eps = math.sqrt(a2 * b2)
    r = dict(u=u, v=v, N=N, q=q, out_dim=W.shape[0])
    r['isometry_err'] = float(np.abs(W.T @ W - np.eye(D)).max())
    errL = errR = miss = 0.0
    E = set(edges)
    for al in range(u):
        for be in range(v):
            F = np.zeros((q, q)); F[al, u + be] = 1.0
            cL = W.T @ port_op(q, nA, D, 'L', F)(W)
            cR = W.T @ port_op(q, nA, D, 'R', F)(W)
            if (al, be) in E:
                e = edges.index((al, be))
                errL = max(errL, np.abs(cL - eps * fs[e].T).max())
                errR = max(errR, np.abs(cR - eps * (1 - b2) * fs[e].T @ P).max())
            else:
                miss = max(miss, np.abs(cL).max(), np.abs(cR).max())
    r.update(left_compression_err=float(errL), right_compression_err=float(errR),
             missing_edge_compression=float(miss))
    # channel identities (8),(9) of exact_route/README via V_L only
    VLr = VL.reshape(q * nA, D, D)
    chP = max(np.abs(np.einsum('kxy,xz,kzw->yw', VLr, fs[j].T @ P, VLr) - (1 - b2) * fs[j].T @ P).max() for j in range(N))
    ch = max(np.abs(np.einsum('kxy,xz,kzw->yw', VLr, fs[j].T, VLr) - (1 - (2 * N - 1) * b2) * fs[j].T).max() for j in range(N))
    r.update(channel_dressed_err=float(chP), channel_undressed_err=float(ch))
    # species covariance of V_L with the flag charges of eq. (17)
    S = np.ones(q, int); C0 = 2 * S
    eu = np.eye(q, dtype=int)
    port_ch = [eu[al] for al in range(u)] + [eu[u + be] for be in range(v)]
    flag_ch = [C0 - eu[al] for al in range(u)] + [C0 + eu[al] - 2 * eu[u + be] for (al, be) in edges]
    w = np.array([eu[al] - eu[u + be] for (al, be) in edges])
    logq = occ @ w
    viol = 0
    V4 = VL.reshape(q, nA, D, D)
    for pb, fa, x, y in zip(*np.nonzero(np.abs(V4) > 1e-14)):
        if not np.array_equal(port_ch[pb] + flag_ch[fa] + logq[x], logq[y] + C0):
            viol += 1
    r['charge_violations_VL'] = viol
    r['flags_nonneg_sum_2q-1'] = bool(all((c >= 0).all() and c.sum() == 2 * q - 1 for c in flag_ch))
    r['flags_distinct'] = len({tuple(c) for c in flag_ch}) == len(flag_ch)
    # upper-species parity of flags is odd, of C0 even
    r['flag_upper_parity_odd'] = bool(all(c[:u].sum() % 2 == 1 for c in flag_ch))
    # products of two port operators on the SAME port compress to at most quadratic operators
    RES[name] = r
    print(name, json.dumps(r))
    assert r['isometry_err'] < 1e-12 and errL < 1e-12 and errR < 1e-12 and miss < 1e-12
    assert chP < 1e-12 and ch < 1e-12 and viol == 0
    return W, fs, P


# --------------------------------------------------------------- transducer on species words
def transducer_stage(u, v, edges, occ_in, Q_sign_right=False):
    """Basis map of CIRCUIT section 1 on actual species words.
    Species labels: 0..u-1 are (1,alpha); u..u+v-1 are (0,beta).
    Input word: N occupation symbols followed by the reference word (2 of each species).
    Returns list over branches of (amplitude, output_word, residual_occ)."""
    N = len(edges); q = u + v
    occ_sym = [edges[e][0] if occ_in[e] else u + edges[e][1] for e in range(N)]
    ref = [s for s in range(q) for _ in range(2)]
    inp = occ_sym + ref
    a = 1 / math.sqrt(2 * u); b = 1 / math.sqrt(2 * N)
    Q = sum(occ_in)
    out = []
    for al in range(u):  # upper branch
        port = al
        pool = list(ref); pool.remove(al)
        flag = sorted(pool)
        out.append((a, ('U', al), [port] + flag + occ_sym, list(occ_in)))
    for e, (al, be) in enumerate(edges):
        Pe = sum(occ_in[:e]) % 2
        new = list(occ_in)
        pool = list(ref)
        if occ_in[e] == 0:   # creation f_e^dagger, flag g_alpha
            new[e] = 1
            pool.remove(al)                     # |1,alpha> moves into data site e
            port = u + be                       # old |0,beta> of data site e goes to the port
            sign = (-1) ** Pe * ((-1) ** Q if Q_sign_right else 1)
        else:                # annihilation f_e, flag h_e
            new[e] = 0
            pool.remove(u + be); pool.remove(u + be)   # data site and port both take |0,beta>
            pool.append(al)                            # old |1,alpha> of the data site joins the flag
            port = u + be
            sign = (-1) ** Pe * ((-1) ** (Q - 1) if Q_sign_right else 1)
        flag = sorted(pool)
        new_sym = [edges[x][0] if new[x] else u + edges[x][1] for x in range(N)]
        out.append((b * sign, ('E', e), [port] + flag + new_sym, new))
    return inp, out


def check_transducer(name, u, v, edges, inputs, compare_V=None, fs=None, P=None):
    N = len(edges); q = u + v
    seen = {}; mism = 0; ninj = 0; nflag = 0
    flags_allowed = {tuple(sorted([s for s in range(q) for _ in range(2) if s != al] + [al])) for al in range(u)}
    D = 2 ** N
    Vt = np.zeros((q, u + N, D, D)) if compare_V is not None else None
    for occ_in in inputs:
        for right in ([False, True] if compare_V is not None else [False]):
            inp, out = transducer_stage(u, v, edges, occ_in, right)
            for amp, br, word, new in out:
                if sorted(word) != sorted(inp):
                    mism += 1
                key = (right, tuple(word))
                if key in seen and seen[key] != (tuple(occ_in), br):
                    ninj += 1
                seen[key] = (tuple(occ_in), br)
                if len(word) != 1 + (2 * q - 1) + N:
                    nflag += 1
                if Vt is not None and not right:
                    x = int(''.join(map(str, occ_in)), 2); y = int(''.join(map(str, new)), 2)
                    fl = br[1] if br[0] == 'U' else (edges[br[1]][0] if occ_in[br[1]] == 0 else u + br[1])
                    Vt[word[0], fl, y, x] += amp
    r = dict(inputs=len(inputs), branches_checked=len(seen), species_multiset_violations=mism,
             injectivity_violations=ninj, word_length_violations=nflag)
    if Vt is not None:
        r['transducer_vs_algebraic_VL'] = float(np.abs(Vt.reshape(-1, D) - compare_V).max())
    RES[name + '_transducer'] = r
    print(name, 'transducer', json.dumps(r))
    assert mism == 0 and ninj == 0 and nflag == 0
    if Vt is not None:
        assert r['transducer_vs_algebraic_VL'] < 1e-12


# --------------------------------------------------------------- local gates
def check_local_gates(q, ell):
    types = [('D', s) for s in range(q)] + [('blank',), ('bit', 0), ('bit', 1)] + \
            [('H', j, c) for j in range(ell) for c in (0, 1)] + \
            [('T', j, c, s) for j in range(ell) for c in (0, 1) for s in range(q)]
    assert len(types) == q + 3 + 2 * (q + 1) * ell
    ch = lambda t: (np.eye(q, dtype=int)[t[-1]] if t[0] in ('D', 'T') else np.zeros(q, int))

    def G_load(x, y):   # |D_s,H_jc> <-> |blank,T_jcs>
        if x[0] == 'D' and y[0] == 'H':
            return ('blank',), ('T', y[1], y[2], x[1])
        if x == ('blank',) and y[0] == 'T':
            return ('D', y[3]), ('H', y[1], y[2])
        return x, y

    def G_swap(x, y):   # |T_j1s, D_t> -> |T_j1t, D_s>
        if x[0] == 'T' and x[2] == 1 and y[0] == 'D':
            return ('T', x[1], 1, y[1]), ('D', x[3])
        return x, y

    def G_query(x, y):  # |D_s,H_jc> -> |D_s,H_j,c xor bit_j(s)>
        if x[0] == 'D' and y[0] == 'H':
            return x, ('H', y[1], y[2] ^ ((x[1] >> y[1]) & 1))
        return x, y

    res = {}
    for name, G in [('load', G_load), ('controlled_exchange', G_swap), ('query', G_query)]:
        imgs = [G(x, y) for x in types for y in types]
        perm = len(set(imgs)) == len(types) ** 2
        invol = all(G(*G(x, y)) == (x, y) for x in types for y in types)
        cons = all(np.array_equal(ch(x) + ch(y), ch(G(x, y)[0]) + ch(G(x, y)[1])) for x in types for y in types)
        res[name] = dict(permutation=perm, involution=invol, charge_conserving=cons)
        assert perm and invol and cons
    # three-gate controlled swap on all three-site words
    ok = True; nwords = 0
    for x in types:
        for y in types:
            for z in types:
                a, b = G_load(x, y); b, c = G_swap(b, z); a, b = G_load(a, b)
                nwords += 1
                if not np.array_equal(ch(x) + ch(y) + ch(z), ch(a) + ch(b) + ch(c)):
                    ok = False
                if x[0] == 'D' and y[0] == 'H' and z[0] == 'D':
                    want = (z, y, x) if y[2] == 1 else (x, y, z)
                    if (a, b, c) != want:
                        ok = False
    res['controlled_swap_all_words'] = dict(words=nwords, correct_and_conserving=ok)
    assert ok
    RES[f'local_gates_q{q}_l{ell}'] = res
    print('local gates', q, ell, json.dumps(res))


# --------------------------------------------------------------- exact uniform amplitudes
def uniform_prefix_state(M, ell):
    """Prepare M^{-1/2} sum_{x<M}|x> on ell bits by the prefix rule of
    EXACT_AMPLITUDES (Hadamard on full prefixes, R_{M,r} on the partial one)."""
    amps = {(): 1.0}
    exceptional = 0
    for r in range(ell):
        new = {}
        width = 2 ** (ell - r)
        nonH = 0
        for pre, amp in amps.items():
            lo = int(''.join(map(str, pre)), 2) * width if pre else 0
            if lo + width <= M:            # full
                c0 = c1 = amp / math.sqrt(2)
            elif lo >= M:                  # empty
                c0, c1 = amp, 0.0
            else:                          # partial
                t = M - lo
                A = min(t, width // 2); B = t - A
                c0, c1 = amp * math.sqrt(A / t), amp * math.sqrt(B / t)
                if B > 0:
                    nonH += 1
            new[pre + (0,)] = new.get(pre + (0,), 0) + c0
            new[pre + (1,)] = new.get(pre + (1,), 0) + c1
        exceptional += nonH
        amps = new
    vec = np.zeros(2 ** ell)
    for pre, amp in amps.items():
        vec[int(''.join(map(str, pre)), 2)] = amp
    return vec, exceptional


def check_amplitudes():
    worst = 0.0; maxexc_ok = True
    for ell in range(1, 9):
        for M in range(1, 2 ** ell + 1):
            vec, exc = uniform_prefix_state(M, ell)
            want = np.zeros(2 ** ell); want[:M] = 1 / math.sqrt(M)
            worst = max(worst, np.abs(vec - want).max())
            maxexc_ok &= exc <= ell
    RES['uniform_amplitudes'] = dict(max_err=float(worst), at_most_ell_exceptional=bool(maxexc_ok))
    print('amplitudes', RES['uniform_amplitudes'])
    assert worst < 1e-12 and maxexc_ok


def fano_channel():
    u, v, edges = projective_plane(2)
    N = len(edges)
    fs, P, occ = jw_ops(N, sparse=True)
    b2 = 1 / (2 * N); ua2 = 0.5
    errs = []
    for j in [0, 7, N - 1]:
        O = fs[j].T @ P
        D = ua2 * O
        for e in range(N):
            D = D + b2 * (fs[e] @ O @ fs[e].T + fs[e].T @ O @ fs[e])
        errs.append(float(abs(D - (1 - b2) * O).max()))
        O2 = fs[j].T
        D2 = ua2 * O2
        for e in range(N):
            D2 = D2 + b2 * (fs[e] @ O2 @ fs[e].T + fs[e].T @ O2 @ fs[e])
        errs.append(float(abs(D2 - (1 - (2 * N - 1) * b2) * O2).max()))
    RES['fano_channel_identities'] = dict(N=N, max_err=max(errs))
    print('fano channel', RES['fano_channel_identities'])
    assert max(errs) < 1e-12


def main():
    hexagon = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 0)]
    W, fs, P = check_graph('six_cycle', 3, 3, hexagon)
    four = [(0, 0), (0, 1), (1, 0), (1, 1)]          # 4-cycle plus isolated upper vertex 2
    check_graph('four_cycle_isolated', 3, 2, four)
    path = [(0, 0), (1, 0), (1, 1)]
    check_graph('path', 2, 2, path)
    # transducer: exhaustive on the hexagon (compare with algebraic V_L), sampled on Fano
    VL = one_port(3, 3, hexagon, fs)
    allocc = [list(map(int, bits)) for bits in itertools.product('01', repeat=6)]
    check_transducer('six_cycle', 3, 3, hexagon, allocc, compare_V=VL)
    u, v, edges = projective_plane(2)
    rng = np.random.default_rng(1)
    occs = [list(rng.integers(0, 2, len(edges))) for _ in range(300)]
    check_transducer('fano_sampled', u, v, edges, occs)
    check_local_gates(3, 2)
    check_local_gates(5, 3)
    check_amplitudes()
    fano_channel()
    Path(__file__).with_name('isometry_results.json').write_text(json.dumps(RES, indent=1) + '\n')


if __name__ == '__main__':
    main()
