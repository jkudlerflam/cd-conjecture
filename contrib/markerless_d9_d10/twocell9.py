"""Two open cells of the 9-letter bond: exact second-order (in the port) effective Hamiltonian on the
code space, its g^2 coefficient, free-fermion test, and the logical blocks mu_L, mu_R, T."""
import numpy as np, itertools, scipy.sparse as sp, scipy.sparse.linalg as spla
from collections import defaultdict
import bond9
from bond9 import IDX, d, HEADS

def cell_words(n):
    out = []
    for h in HEADS:
        for j in range(n+1):
            for bits in itertools.product((0,1), repeat=n):
                w = [f'L{a}' for a in bits[:j]] + [h] + [f'R{a}' for a in bits[j:]]
                out.append((tuple(IDX[x] for x in w), h, j, bits))
    return out

def build(n, p, cw=None):
    cw = cw or cell_words(n)
    padL, padR = (IDX['R0'],), (IDX['L0'],)
    states = [(padL+a[0]+b[0]+padR, (a[1],b[1]), (a[2],b[2]), a[3]+b[3]) for a in cw for b in cw]
    index = {s[0]:i for i,s in enumerate(states)}
    w, o, bad = bond9.bond(p)
    by_in = defaultdict(list)
    for r,c in zip(*np.nonzero(w)): by_in[int(c)].append((int(r), w[r,c]))
    rows, cols, vals, leaks = [], [], [], 0
    Ls = 2*(n+1)+2
    for i,s in enumerate(states):
        word = s[0]
        rows.append(i); cols.append(i); vals.append(sum(o[x] for x in word))
        for site in range(Ls-1):
            for r,v in by_in.get(word[site]*d + word[site+1], ()):
                k = index.get(word[:site] + (r//d, r%d) + word[site+2:])
                if k is None: leaks += 1
                else: rows.append(k); cols.append(i); vals.append(v)
    return states, sp.csr_matrix((vals,(rows,cols)), shape=(len(states),)*2), leaks

def analyse(n, p):
    cw = cell_words(n)
    states, Hfull, leaks = build(n, p, cw)
    _, Hnog, _ = build(n, dict(p, g=0.0), cw)
    _, Hnol, _ = build(n, dict(p, lam=0.0), cw)
    V = Hfull - Hnol; Wg = Hfull - Hnog; H0 = Hnog - V
    kind = np.array([0 if s[1]==('W','W') else (1 if 'W' not in s[1] else 2) for s in states])
    vac = np.where(kind==0)[0]; exc = np.where(kind==1)[0]
    ncode = 4**n
    vv = np.sqrt(2/(n+2))*np.sin(np.pi*(np.arange(n+1)+1)/(n+2))
    E0 = 2*(-2*np.cos(np.pi/(n+2)) + p.get('nu',0.0))
    bitlist = list(itertools.product((0,1), repeat=2*n)); bidx = {b:i for i,b in enumerate(bitlist)}
    P = sp.lil_matrix((len(vac), ncode))
    for a,i in enumerate(vac):
        s = states[i]; P[a, bidx[s[3]]] = vv[s[2][0]]*vv[s[2][1]]
    P = P.tocsc()
    Hv = H0[vac][:,vac]
    fbits = np.array([p.get('f',0.0)*sum(b) for b in bitlist])
    resid = abs((Hv@P - P@sp.diags(E0+fbits))).max()
    first = abs(P.T @ (V[vac][:,vac] @ P)).max() if V[vac][:,vac].nnz else 0.0
    VP = (V[exc][:,vac] @ P).toarray().astype(complex)
    A0 = (H0[exc][:,exc] - E0*sp.identity(len(exc))).tocsc()
    lu = spla.splu(A0)
    X0 = lu.solve(VP); Wee = Wg[exc][:,exc]
    X1 = lu.solve(Wee @ X0); X2 = lu.solve(Wee @ X1)
    h0 = -(VP.conj().T @ X0); h1 = -(VP.conj().T @ X1); h2 = -(VP.conj().T @ X2)
    return dict(nstates=len(states), leaks=leaks, code_resid=resid, first=first, h0=h0, h1=h1, h2=h2, bitlist=bitlist)

def one_particle(h2, n, bitlist):
    m = 2*n
    idx1 = [bitlist.index(tuple(1 if i==a else 0 for i in range(m))) for a in range(m)]
    i0 = bitlist.index(tuple([0]*m)); vacE = h2[i0,i0]
    return h2[np.ix_(idx1, idx1)] - vacE*np.eye(m), vacE

def jw_many_body(h, vacE, bitlist):
    m = h.shape[0]; N = len(bitlist); bidx = {b:i for i,b in enumerate(bitlist)}
    M = np.zeros((N,N), complex)
    for i,b in enumerate(bitlist):
        M[i,i] += vacE
        for bb in range(m):
            if not b[bb]: continue
            s1 = (-1)**sum(b[:bb]); t = list(b); t[bb] = 0
            for a in range(m):
                if t[a]: continue
                u = list(t); u[a] = 1
                M[bidx[tuple(u)], i] += h[a,bb]*s1*(-1)**sum(t[:a])
    return M

def fermi_points(mu, T, E=0.0):
    import scipy.linalg as sl
    n = mu.shape[0]
    A = np.block([[np.zeros((n,n)), np.eye(n)],[-T.conj().T, -(mu-E*np.eye(n))]])
    B = np.block([[np.eye(n), np.zeros((n,n))],[np.zeros((n,n)), T]])
    return sl.eigvals(A, B)

if __name__ == '__main__':
    import sys
    n = int(sys.argv[1]) if len(sys.argv)>1 else 2
    np.set_printoptions(precision=4, suppress=False, linewidth=220)
    for label, td, D in [('td=0, D=30', 0.0, 30.0), ('td=0.02, D=30', 0.02, 30.0), ('td=0.3, D=4', 0.3, 4.0)]:
        p = dict(lam=1.0, g=1.0, Delta=4-2*np.cos(np.pi/(n+2)), D=D, nu=0.0, f=0.0, phiA=np.pi/(2*(n+1)), td=td)
        r = analyse(n, p)
        h, vacE = one_particle(r['h2'], n, r['bitlist'])
        M = jw_many_body(h, vacE, r['bitlist'])
        muL, muR, T = h[:n,:n], h[n:,n:], h[:n,n:]; mu = muL+muR; sc = np.abs(T).max()
        h0 = r['h0']
        print(f"--- n={n} {label}: states={r['nstates']} leaks={r['leaks']} code_resid={r['code_resid']:.1e} first={r['first']:.1e}")
        print(f"   g^0: bit-dependence {np.ptp(np.diag(h0).real):.1e}, offdiag {np.abs(h0-np.diag(np.diag(h0))).max():.1e};  g^1: {np.abs(r['h1']).max():.1e}")
        print(f"   g^2: |h2|={np.abs(r['h2']).max():.3e}  herm={np.abs(r['h2']-r['h2'].conj().T).max():.1e}  Gaussian residual/|T|={np.abs(r['h2']-M).max()/sc:.1e}")
        ph = np.angle(np.trace(T)); Th = T*np.exp(-1j*ph)
        ev = np.linalg.eigvalsh((Th+Th.conj().T)/2)
        print(f"   |mu|/|T|={np.abs(mu).max()/sc:.2e}  phase(T)={ph/np.pi:.4f}pi  antiHerm(T e^-i ph)/|T|={np.abs(Th-Th.conj().T).max()/2/sc:.2e}  eig Herm/|T|={ev/sc}")
        z = fermi_points(mu, T); on = np.abs(np.abs(z)-1) < 1e-6
        ks = np.sort(np.angle(z[on]))
        print(f"   Fermi points on the circle: {on.sum()} of {2*n};  k/pi = {ks/np.pi};  reciprocal: {np.allclose(np.sort(ks), np.sort(-ks), atol=1e-7)}")
