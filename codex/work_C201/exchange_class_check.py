"""Exchange class (Astra A3 Sec.1): site = H0 (+) H1, bond h_W = -t(W+W^dag) + V n(x)n + (mu/2)(n(x)1+1(x)n),
W: H1(x)H0 -> H0(x)H1 arbitrary unitary. Checks: (i) E_L(W) >= E_L(scalar) = min_m [scalar hard-core ring with
hopping -t, no twist]; (ii) open-chain spectra equal the scalar ones with multiplicity; (iii) M(h_W) = M(h_sc)
numerically by gauge minimisation; (iv) ring spectrum contained in union of twisted scalar spectra (twist grid).
Claude Fable 5.2 checking Astra gpt-6.2-astra, work_C201."""
import numpy as np, itertools
from scipy.stats import unitary_group
from scipy.optimize import minimize
from potts_xxz_inclusion import xxz_sector
rng = np.random.default_rng(7)

def bond(d0, d1, W, t, V, mu):
    d = d0+d1; n = np.diag([0.]*d0+[1.]*d1); I = np.eye(d)
    h = V*np.kron(n, n)+0.5*mu*(np.kron(n, I)+np.kron(I, n))
    Wfull = np.zeros((d*d, d*d), dtype=complex)
    # input |a>|b> with a in H1 (index d0+i), b in H0 (index j); output |b'>|a'> b' in H0, a' in H1
    for i in range(d1):
        for j in range(d0):
            col = (d0+i)*d+j
            for jp in range(d0):
                for ip in range(d1):
                    row = jp*d+(d0+ip)
                    Wfull[row, col] = W[jp*d1+ip, i*d0+j]
    return h-t*(Wfull+Wfull.conj().T)

def chain_H(h, d, L, ring=True):
    dim = d**L; H = np.zeros((dim, dim), dtype=complex)
    h4 = h.reshape(d, d, d, d)
    for j in range(L if ring else L-1):
        a, b = j, (j+1) % L
        # apply h on sites (a,b) using tensor reshape
        T = np.eye(dim, dtype=complex).reshape([d]*L+[dim])
        T = np.moveaxis(T, [a, b], [0, 1])
        sh = T.shape
        T = (h @ T.reshape(d*d, -1)).reshape(sh)
        T = np.moveaxis(T, [0, 1], [a, b])
        H += T.reshape(dim, dim)
    return H

def scalar_ring_min(L, t, V, mu, phi=0.0, ring=True):
    """hard-core ring: -t(hop) + V n n + mu n ; via XXZ sector builder: xxz_sector gives hop coefficient +1 and
    0.5*Delta*(zz-1). Build directly instead."""
    out = []
    for m in range(L+1):
        confs = list(itertools.combinations(range(L), m)); idx = {c: i for i, c in enumerate(confs)}
        H = np.zeros((len(confs),)*2, dtype=complex)
        for c in confs:
            i = idx[c]; occ = set(c)
            for b in range(L if ring else L-1):
                a, a2 = b, (b+1) % L
                if a in occ and a2 in occ: H[i, i] += V
                H[i, i] += 0.5*mu*((a in occ)+(a2 in occ))
                if (a in occ) != (a2 in occ):
                    new = set(occ)
                    if a in occ: new.remove(a); new.add(a2); ph = np.exp(1j*phi) if (ring and b == L-1) else 1
                    else: new.remove(a2); new.add(a); ph = np.exp(-1j*phi) if (ring and b == L-1) else 1
                    H[idx[tuple(sorted(new))], i] += -t*ph
        out.append(np.linalg.eigvalsh(H))
    return out

def Mgauge(h, d):
    def spread(x):
        A = x[:d*d].reshape(d, d)+1j*x[d*d:].reshape(d, d); X = (A+A.conj().T)/2
        ev = np.linalg.eigvalsh(h+np.kron(X, np.eye(d))-np.kron(np.eye(d), X))
        return (ev[-1]-ev[0])/2
    best = spread(np.zeros(2*d*d))
    for trial in range(6):
        r = minimize(spread, 0.3*rng.standard_normal(2*d*d), method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-7, 'fatol': 1e-9})
        best = min(best, r.fun)
    return best, spread(np.zeros(2*d*d))

if __name__ == "__main__":
    for (d0, d1, Ls) in ((2, 2, (3, 4, 5)), (1, 2, (4, 5, 6)), (2, 1, (4, 5))):
        d = d0+d1
        for (t, V, mu) in ((1.0, 0.0, 0.0), (1.0, 0.8, -0.3), (1.0, -1.2, 0.5)):
            W = unitary_group.rvs(d0*d1, random_state=rng) if d0*d1 > 1 else np.array([[np.exp(0.4j)]])
            h = bond(d0, d1, W, t, V, mu)
            hs = bond(1, 1, np.array([[1.0]]), t, V, mu)
            Mg, M0 = Mgauge(h, d); Msg, Ms0 = Mgauge(hs, 2)
            for L in Ls:
                if d**L > 1100: continue
                ev = np.linalg.eigvalsh(chain_H(h, d, L, True))
                sc = scalar_ring_min(L, t, V, mu)
                Esc = min(e[0] for e in sc)
                evo = np.linalg.eigvalsh(chain_H(h, d, L, False))
                sco = np.unique(np.round(np.concatenate(scalar_ring_min(L, t, V, mu, ring=False)), 8))
                open_mis = max(np.min(np.abs(sco-e)) for e in np.unique(np.round(evo, 8)))
                # twisted inclusion on a grid of 720 twists
                grid = np.concatenate([np.concatenate(scalar_ring_min(L, t, V, mu, phi=2*np.pi*k/360)) for k in range(360)]) if d**L <= 300 else None
                inc = max(np.min(np.abs(grid-e)) for e in np.unique(np.round(ev, 8))) if grid is not None else float('nan')
                print(f"(d0,d1)=({d0},{d1}) t,V,mu={t},{V},{mu} L={L}: E_L(W)={ev[0].real:.8f} E_L(sc)={Esc:.8f} diff={ev[0]-Esc:+.2e} "
                      f"open-spectrum mismatch={open_mis:.1e} ring-in-twisted(grid) mismatch={inc:.1e}  M(W): X=0 {M0:.5f}, optimised {Mg:.5f}; M(sc) {Ms0:.5f}/{Msg:.5f}")
