"""Periodic check of the Bloch assembly: M = 3 cells of the 9-letter bond on a genuine ring (no padding).
The one-particle block of the lambda^2 g^2 effective Hamiltonian is compared with
delta_{mm'} (mu_L + mu_R) + T on (m, m+1) + T^dag on (m+1, m), with T, mu_L, mu_R taken from the padded
two-cell calculation.  This tests the padding, the ring geometry and the assembly at once."""
import numpy as np, itertools, scipy.sparse as sp, scipy.sparse.linalg as spla, json, sys
from collections import defaultdict
import bond9
from bond9 import IDX, d, HEADS
from twocell9 import cell_words, analyse, one_particle

def build_ring(n, M, p, cw):
    states = []
    for combo in itertools.product(range(len(cw)), repeat=M):
        cs = [cw[i] for i in combo]
        nexc = sum(c[1] != 'W' for c in cs)
        if nexc == 0 or (nexc == 2 and any(cs[m][1] != 'W' and cs[(m+1) % M][1] != 'W' for m in range(M))):
            states.append((sum((c[0] for c in cs), ()), tuple(c[1] for c in cs), tuple(c[2] for c in cs), sum((c[3] for c in cs), ())))
    index = {s[0]: i for i, s in enumerate(states)}
    w, o, bad = bond9.bond(p)
    by_in = defaultdict(list)
    for r, c in zip(*np.nonzero(w)): by_in[int(c)].append((int(r), w[r, c]))
    rows, cols, vals = [], [], []
    Ls = M*(n+1); dropped = 0
    for i, s in enumerate(states):
        word = s[0]
        rows.append(i); cols.append(i); vals.append(sum(o[x] for x in word))
        for site in range(Ls):
            j = (site+1) % Ls
            for r, v in by_in.get(word[site]*d + word[j], ()):
                nw = list(word); nw[site] = r//d; nw[j] = r % d
                k = index.get(tuple(nw))
                if k is None: dropped += 1          # leaves the truncated space (a third excited cell): higher order in lambda
                else: rows.append(k); cols.append(i); vals.append(v)
    return states, sp.csr_matrix((vals, (rows, cols)), shape=(len(states),)*2), dropped

def ring_blocks(n, M, p):
    cw = cell_words(n)
    states, Hfull, _ = build_ring(n, M, p, cw)
    _, Hnog, _ = build_ring(n, M, dict(p, g=0.0), cw)
    _, Hnol, _ = build_ring(n, M, dict(p, lam=0.0), cw)
    V = Hfull - Hnol; Wg = Hfull - Hnog; H0 = Hnog - V
    kind = np.array([0 if all(h == 'W' for h in s[1]) else 1 for s in states])
    vac = np.where(kind == 0)[0]; exc = np.where(kind == 1)[0]
    vv = np.sqrt(2/(n+2))*np.sin(np.pi*(np.arange(n+1)+1)/(n+2))
    E0 = M*(-2*np.cos(np.pi/(n+2)) + p.get('nu', 0.0))
    # code states with zero or one unit
    codes = [tuple([0]*(M*n))] + [tuple(1 if i == a else 0 for i in range(M*n)) for a in range(M*n)]
    cidx = {b: i for i, b in enumerate(codes)}
    P = sp.lil_matrix((len(vac), len(codes)))
    for a, i in enumerate(vac):
        s = states[i]
        if s[3] in cidx: P[a, cidx[s[3]]] = np.prod([vv[j] for j in s[2]])
    P = P.tocsc()
    VP = (V[exc][:, vac] @ P).toarray().astype(complex)
    lu = spla.splu((H0[exc][:, exc] - E0*sp.identity(len(exc))).tocsc())
    X0 = lu.solve(VP); Wee = Wg[exc][:, exc]
    X2 = lu.solve(Wee @ lu.solve(Wee @ X0))
    h2 = -(VP.conj().T @ X2)
    return h2[1:, 1:] - h2[0, 0]*np.eye(M*n), len(states)

if __name__ == '__main__':
    out = []
    for n in (1, 2):
        M = 3
        for a, D, td in [(2.6, 3.0, 0.07), (2.2, 0.4, 0.3)]:
            p = dict(lam=1.0, g=1.0, Delta=a-2*np.cos(np.pi/(n+2)), D=D, nu=0.0, f=0.0, phiA=np.pi/(2*max(n,1)), td=td)
            hr, nst = ring_blocks(n, M, p)
            r = analyse(n, p); h, vacE = one_particle(r['h2'], n, r['bitlist'])
            muL, muR, T = h[:n, :n], h[n:, n:], h[:n, n:]
            pred = np.zeros((M*n, M*n), complex)
            for m in range(M):
                pred[m*n:(m+1)*n, m*n:(m+1)*n] += muL + muR
                m2 = (m+1) % M
                pred[m*n:(m+1)*n, m2*n:(m2+1)*n] += T
                pred[m2*n:(m2+1)*n, m*n:(m+1)*n] += T.conj().T
            sc = np.abs(T).max()
            ev1 = np.linalg.eigvalsh((hr+hr.conj().T)/2); ev2 = np.linalg.eigvalsh(pred)
            rec = dict(n=n, cells=M, a=a, D=D, td=td, ring_states=nst, max_entry_difference_over_T=float(np.abs(hr-pred).max()/sc),
                       max_level_difference_over_T=float(np.abs(ev1-ev2).max()/sc), mu_over_T=float(np.abs(muL+muR).max()/sc))
            out.append(rec); print(rec, flush=True)
    json.dump(out, open('ring_three_cells.json', 'w'), indent=1)
