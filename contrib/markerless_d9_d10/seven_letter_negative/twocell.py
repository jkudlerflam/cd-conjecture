"""Two open cells of the 7-letter bond: exact second-order (in the port) effective Hamiltonian
on the code space, g^2 coefficient, and a free-fermion (Gaussianity) test."""
import numpy as np, itertools, scipy.sparse as sp, scipy.sparse.linalg as spla
from bond7 import bond, LET, IDX, d

def cell_words(n, heads):
    out = []
    for h in heads:
        for j in range(n+1):
            for bits in itertools.product((0,1), repeat=n):
                w = [f'L{a}' for a in bits[:j]] + [h] + [f'R{a}' for a in bits[j:]]
                out.append((tuple(IDX[x] for x in w), h, j, bits))
    return out

def build(n, p):
    cw = cell_words(n, 'WAB')
    states = [(a[0]+b[0], (a[1],b[1]), (a[2],b[2]), a[3]+b[3]) for a in cw for b in cw]
    index = {s[0]:i for i,s in enumerate(states)}
    w, o, bad = bond(p)
    N = len(states); Ls = 2*(n+1)
    rows, cols, vals = [], [], []
    leaks = 0
    nz = [(int(r), int(c), w[r,c]) for r,c in zip(*np.nonzero(w))]
    from collections import defaultdict
    by_in = defaultdict(list)
    for r,c,v in nz: by_in[c].append((r,v))
    for i,s in enumerate(states):
        word = s[0]
        rows.append(i); cols.append(i); vals.append(sum(o[x] for x in word))
        for site in range(Ls-1):
            c = word[site]*d + word[site+1]
            for r,v in by_in.get(c, ()):
                nw = word[:site] + (r//d, r%d) + word[site+2:]
                k = index.get(nw)
                if k is None: leaks += 1
                else: rows.append(k); cols.append(i); vals.append(v)
    H = sp.csr_matrix((vals,(rows,cols)), shape=(N,N))
    return states, H, leaks

def analyse(n, p, verbose=True):
    p0 = dict(p); p0['gR'] = 0; p0['gL'] = 0
    pl = dict(p); pl['lam'] = 0
    states, Hfull, leaks = build(n, p)          # everything
    _, Hnog, _ = build(n, p0)                   # g = 0
    _, Hnol, _ = build(n, pl)                   # lam = 0
    V = (Hfull - Hnol)                          # port only
    Wg = (Hfull - Hnog)                         # conversions only
    H0 = Hnog - V                               # no port, no conversion
    heads = np.array([0 if s[1]==('W','W') else (1 if 'W' not in s[1] else 2) for s in states])
    vac = np.where(heads==0)[0]; exc = np.where(heads==1)[0]
    # vacuum sector ground space
    Hv = H0[vac][:,vac].toarray()
    ev, U = np.linalg.eigh(Hv)
    ncode = 2**(2*n)
    deg = np.sum(ev < ev[0]+1e-9)
    assert deg == ncode, (deg, ncode)
    E0 = ev[0]
    # code basis labelled by bits: P = |v>|v> x |bits>
    vv = np.sqrt(2/(n+2))*np.sin(np.pi*(np.arange(n+1)+1)/(n+2))
    bitlist = list(itertools.product((0,1), repeat=2*n))
    bidx = {b:i for i,b in enumerate(bitlist)}
    P = np.zeros((len(vac), ncode))
    for a,i in enumerate(vac):
        s = states[i]; P[a, bidx[s[3]]] = vv[s[2][0]]*vv[s[2][1]]
    assert np.allclose(Hv@P, E0*P, atol=1e-10), 'code is not the clock ground space'
    first = P.T @ (V[vac][:,vac] @ P)
    VP = V[exc][:,vac] @ P                      # excited <- code
    A0 = (H0[exc][:,exc] - E0*sp.identity(len(exc))).tocsc()
    lu = spla.splu(A0)
    X0 = lu.solve(VP.astype(complex))
    Wee = Wg[exc][:,exc]
    X1 = lu.solve(Wee @ X0); X2 = lu.solve(Wee @ X1)
    h0 = -(VP.conj().T @ X0)                    # order g^0
    h2 = -(VP.conj().T @ X2)                    # coefficient of g^2 (in units where gR,gL as given)
    return dict(states=states, leaks=leaks, E0=E0, first=np.abs(first).max(), h0=h0, h2=h2,
                bitlist=bitlist, gap=ev[deg]-ev[0], exc_gap=None)

def one_particle(h2, n, bitlist):
    m = 2*n
    idx1 = [bitlist.index(tuple(1 if i==a else 0 for i in range(m))) for a in range(m)]
    vacE = h2[bitlist.index(tuple([0]*m)), bitlist.index(tuple([0]*m))]
    h = h2[np.ix_(idx1, idx1)] - vacE*np.eye(m)
    return h, vacE

def jw_many_body(h, vacE, bitlist):
    """second-quantised sum_ab h_ab c_a^dag c_b + vacE in the bit basis (JW order = site order)"""
    m = h.shape[0]; N = len(bitlist); bidx = {b:i for i,b in enumerate(bitlist)}
    M = np.zeros((N,N), complex)
    for i,b in enumerate(bitlist):
        M[i,i] += vacE
        for bb in range(m):
            if not b[bb]: continue
            s1 = (-1)**sum(b[:bb]); t = list(b); t[bb] = 0
            for a in range(m):
                if t[a]: continue
                s2 = (-1)**sum(t[:a]); u = list(t); u[a] = 1
                M[bidx[tuple(u)], i] += h[a,bb]*s1*s2
    return M

if __name__ == '__main__':
    import sys
    n = int(sys.argv[1]) if len(sys.argv)>1 else 2
    np.set_printoptions(precision=5, suppress=True, linewidth=200)
    for label, gL in [('right rung only', 0.0), ('two rungs, real', 0.6), ('two rungs, flux', 0.6*np.exp(1.1j))]:
        p = dict(lam=1.0, gR=1.0, gL=gL, DA=0.7, DB=0.7, nu=0.0, f=0.0, chi=0.0)
        r = analyse(n, p)
        h, vacE = one_particle(r['h2'], n, r['bitlist'])
        M = jw_many_body(h, vacE, r['bitlist'])
        print(f"--- n={n} {label}: states={len(r['states'])} leaks={r['leaks']} first-order={r['first']:.1e} "
              f"h0 spread={np.ptp(np.diag(r['h0']).real):.1e} offdiag h0={np.abs(r['h0']-np.diag(np.diag(r['h0']))).max():.1e}")
        print("   hermiticity of h2:", np.abs(r['h2']-r['h2'].conj().T).max(), " |h2|:", np.abs(r['h2']).max(),
              " Gaussian residual:", np.abs(r['h2']-M).max())
        if n <= 3:
            print("   mu_L =\n", h[:n,:n]); print("   mu_R =\n", h[n:,n:]); print("   T (left<-right) =\n", h[:n,n:])
