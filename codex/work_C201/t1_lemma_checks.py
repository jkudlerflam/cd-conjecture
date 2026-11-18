"""Numerical checks of two lemmas inside Astra's uniform-decay proof (Theorem A).
(a) Marginal-correcting channel, eq. (10): for states r,s on B with z = ||r-s||_1/2 there is a channel D with D(r)=s
    and ||(id (x) D)(xi) - xi||_1 <= 2 sqrt(z) + z for every extension xi of r.
(b) Covariance step, eqs. (4)-(7), on ring ground states: with t-site windows, groups U,V,W of R windows,
    tau_ends = avg rho_{U_i W_j}, eta_ends = sum_v p_v gamma_U(v) (x) gamma_W(v) (v = IC-measurement record on V),
    check  E||Xbar_U - gamma_U(v)||_2^2 <= 4(a_R - a_2R)  and  ||tau - eta||_1 <= 4 q (a_R - a_2R),
    where a_R = q^4/R + R^-2 sum_{i != j} tr(rho_{ij} SWAP).
Claude Fable 5.2 checking Astra gpt-6.2-astra, work_C201."""
import numpy as np, itertools, scipy.linalg as sla
rng = np.random.default_rng(3)

def rand_state(n, rank=None):
    rank = rank or n
    A = rng.standard_normal((n, rank))+1j*rng.standard_normal((n, rank)); r = A@A.conj().T
    return r/np.trace(r).real
def tn(X): return np.sum(np.abs(np.linalg.eigvalsh((X+X.conj().T)/2)))
def psqrt(X):
    w, v = np.linalg.eigh((X+X.conj().T)/2); w = np.clip(w, 0, None); return (v*np.sqrt(w))@v.conj().T
def pinv_sqrt(X, tol=1e-12):
    w, v = np.linalg.eigh((X+X.conj().T)/2); wi = np.where(w > tol, 1/np.sqrt(np.clip(w, tol, None)), 0); return (v*wi)@v.conj().T

if __name__ == "__main__":
    print("(a) marginal-correcting channel")
    worst = 0
    for trial in range(400):
        dA, dB = rng.integers(2, 4), rng.integers(2, 5)
        xi = rand_state(dA*dB, rank=rng.integers(1, dA*dB+1))
        r = np.einsum('abac->bc', xi.reshape(dA, dB, dA, dB))
        eps = rng.choice([1e-3, 1e-2, 0.1, 0.5, 1.0])
        s = (1-eps)*r+eps*rand_state(dB, rank=rng.integers(1, dB+1)); s /= np.trace(s).real
        z = tn(r-s)/2
        w, v = np.linalg.eigh(r-s); pos = (v*np.clip(w, 0, None))@v.conj().T
        T = s+pos
        K = psqrt(s)@pinv_sqrt(T)
        res = s-K@r@K.conj().T
        out_r = K@r@K.conj().T+res          # = s by construction; check trace of weight
        wgt = (np.trace(r)-np.trace(K@r@K.conj().T)).real
        assert abs(np.trace(res).real-wgt) < 1e-9 and np.linalg.eigvalsh(res).min() > -1e-9
        omega = res/wgt if wgt > 1e-14 else np.eye(dB)/dB
        KK = np.kron(np.eye(dA), K)
        comp = np.kron(np.eye(dA), np.eye(dB)-K.conj().T@K)
        xiA_w = np.einsum('abcb->ac', (psqrt(comp)@xi@psqrt(comp)).reshape(dA, dB, dA, dB))
        out = KK@xi@KK.conj().T+np.kron(xiA_w, omega)
        assert abs(np.trace(out).real-1) < 1e-9
        assert tn(np.einsum('abac->bc', out.reshape(dA, dB, dA, dB))-s) < 1e-8
        ratio = tn(out-xi)/(2*np.sqrt(z)+z)
        worst = max(worst, ratio)
    print(f"   400 random cases: D(r)=s exact, max ||out-xi||_1/(2 sqrt z + z) = {worst:.4f} (must be <= 1)")

def ring_gs(L, kind):
    X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex); I = np.eye(2)
    if kind == 'heis': h = (np.kron(X, X)+np.kron(Y, Y)+np.kron(Z, Z)).real/4
    if kind == 'xx': h = (np.kron(X, X)+np.kron(Y, Y)).real/4
    if kind == 'tfim': h = (-np.kron(Z, Z)-0.5*(np.kron(X, I)+np.kron(I, X))).real
    if kind == 'ghz':
        v = np.zeros(2**L); v[0] = v[-1] = 2**-0.5; return v
    import scipy.sparse as sp, scipy.sparse.linalg as spl
    H = sp.csr_matrix((2**L, 2**L))
    hs = sp.csr_matrix(h)
    for j in range(L-1):
        H = H+sp.kron(sp.kron(sp.identity(2**j), hs), sp.identity(2**(L-j-2)), format='csr')
    # closing bond via cyclic shift of the last term
    perm = np.array([((i << 1) & (2**L-1)) | (i >> (L-1)) for i in range(2**L)])
    Hl = sp.kron(sp.identity(2**(L-2)), hs, format='csr').tocoo()
    H = H+sp.csr_matrix((Hl.data, (perm[Hl.row], perm[Hl.col])), shape=(2**L, 2**L))
    w, v = spl.eigsh(H, k=1, which='SA'); return v[:, 0]

def rdm(psi, L, sites):
    T = psi.reshape([2]*L); rest = [i for i in range(L) if i not in sites]
    T = np.transpose(T, list(sites)+rest).reshape(2**len(sites), -1); return T@T.conj().T

def hs_basis(q):
    B = []
    for i in range(q):
        E = np.zeros((q, q), complex); E[i, i] = 1; B.append(E)
    for i in range(q):
        for j in range(i+1, q):
            E = np.zeros((q, q), complex); E[i, j] = E[j, i] = 2**-0.5; B.append(E)
            E = np.zeros((q, q), complex); E[i, j] = -1j*2**-0.5; E[j, i] = 1j*2**-0.5; B.append(E)
    return B

if __name__ == "__main__":
    print("(b) covariance step on ring ground states")
    for kind, L in (('heis', 12), ('xx', 12), ('tfim', 12), ('ghz', 12)):
        psi = ring_gs(L, kind)
        for t, R in ((2, 1), (1, 1), (1, 2)):
            q = 2**t
            if 3*R*t > L: continue
            G = hs_basis(q)
            win = lambda g, i: list(range((g*R+i)*t, (g*R+i+1)*t))       # group g in {0:U,1:V,2:W}
            SW = sum(np.kron(g_, g_) for g_ in G)
            def a_of(Rr, g0):
                wins = [list(range((g0*R+i)*t, (g0*R+i+1)*t)) for i in range(Rr)]
                tot = q**4/Rr
                for i in range(Rr):
                    for j in range(Rr):
                        if i != j: tot += np.trace(rdm(psi, L, wins[i]+wins[j])@SW).real/Rr**2
                return tot
            aR, a2R = a_of(R, 0), a_of(2*R, 0)
            tau = sum(rdm(psi, L, win(0, i)+win(2, j)) for i in range(R) for j in range(R))/R**2
            # measurement record on V: outcomes (alpha,sign) per window; effects F=(I +- G)/(2 q^2)
            effs = [((np.eye(q)+sg*g_)/(2*q*q), sg*q*q*g_) for g_ in G for sg in (+1, -1)]
            Vs = [s for i in range(R) for s in win(1, i)]
            Us = [win(0, i) for i in range(R)]; Ws = [win(2, j) for j in range(R)]
            eta = np.zeros((q*q, q*q), complex); mse = 0.0; ptot = 0.0
            # joint state of (U-window i, V) and (W-window j, V) for all i,j
            rhoUV = [rdm(psi, L, Us[i]+Vs) for i in range(R)]; rhoWV = [rdm(psi, L, Ws[j]+Vs) for j in range(R)]
            rhoV = rdm(psi, L, Vs)
            for outcome in itertools.product(range(len(effs)), repeat=R):
                F = effs[outcome[0]][0]; Xv = effs[outcome[0]][1]/R
                for o in outcome[1:]:
                    F = np.kron(F, effs[o][0])
                XbarV = sum(effs[o][1] for o in outcome)/R
                p = np.trace(rhoV@F).real
                if p < 1e-15: continue
                gU = sum(np.einsum('aibj,ji->ab', rhoUV[i].reshape(q, q**R, q, q**R), F) for i in range(R))/(R*p)
                gW = sum(np.einsum('aibj,ji->ab', rhoWV[j].reshape(q, q**R, q, q**R), F) for j in range(R))/(R*p)
                eta += p*np.kron(gU, gW); ptot += p
                # E||Xbar_U - gamma_U(v)||^2 = E||Xbar_U||^2 - E||gamma_U||^2 ; accumulate E||gamma_U||^2
                mse -= p*np.sum(np.abs(gU)**2)
            mse += aR
            lhs = tn(tau-eta)
            print(f"   {kind:5s} L={L} t={t} R={R}: sum p={ptot:.6f}  a_R-a_2R={aR-a2R:.4f}  MSE(U|v)={mse:.4f} <= 4(a_R-a_2R)={4*(aR-a2R):.4f}: {mse <= 4*(aR-a2R)+1e-9}"
                  f"   ||tau-eta||_1={lhs:.5f} <= q*MSE-type bound {q*np.sqrt(mse*mse):.4f} and <= 4q(a_R-a_2R)={4*q*(aR-a2R):.3f}")
