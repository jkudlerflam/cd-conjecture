"""C186: few-body problem of the general second-order cusp model (PROOF 2), number not conserved.

Units |t2| = 1, spin language (hard-core bosons, no Jordan-Wigner signs), blockade (no two adjacent particles).
    H = sum_j P[ w b+_{j+1} b_j + h.c.  - (b+_{j+2} b_j + h.c.)  + V2 n_j n_{j+2} + V3 n_j n_{j+3} - mu n_j
                 + ( r b+_j (n_{j+2} - n_{j-2}) + h.c. ) ]P,
    V2 = -2 + v,  V3 = -g  (g = 2/3 for every cusp bond),  r = sqrt(v) for every cusp bond (R2), r real by gauge.
The vacuum is an exact eigenstate (the creation term needs a particle at distance 2).

Infinite line, total momentum K. Basis |K;c> = sum_x e^{iKx} T_x |c>, c = sorted positions with the leftmost at 0,
extent max(c) <= R (hard wall in relative coordinates), particle number 1..Nmax. Both truncations are subspace
restrictions, so every eigenvalue is a variational upper bound for the untruncated problem at that K.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, itertools, sys, json, time
from scipy.optimize import minimize_scalar, brentq


def configs(N, R):
    if N == 1:
        return [(0,)]
    out = []
    # gaps >= 2, sum <= R
    def rec(prefix, last, left):
        if left == 0:
            out.append(tuple(prefix))
            return
        for x in range(last + 2, R + 1):
            # remaining particles need 2 sites each
            if x + 2 * (left - 1) > R:
                break
            rec(prefix + [x], x, left - 1)
    rec([0], 0, N - 1)
    return out


class Model:
    def __init__(self, u, phi, v, g=2/3, r=None, mu=0.0):
        self.w = u * np.exp(1j * phi)
        self.V2 = -2.0 + v
        self.V3 = -g
        self.r = np.sqrt(max(v, 0.0)) if r is None else r
        self.mu = mu

    def eps(self, k):
        return 2 * abs(self.w) * np.cos(k - np.angle(self.w)) - 2 * np.cos(2 * k)

    def actions(self, c, Nmax):
        """yield (new_positions_tuple_unsorted_ok, amplitude); diagonal returned separately"""
        S = set(c)
        N = len(c)
        dia = -self.mu * N
        for a, b_ in itertools.combinations(c, 2):
            dd = abs(a - b_)
            if dd == 2:
                dia += self.V2
            elif dd == 3:
                dia += self.V3
        out = []
        w = self.w
        for x in c:
            rest = S - {x}
            # NN hops (x+1, x-1 empty by blockade)
            if (x + 2) not in S:
                out.append((rest | {x + 1}, w))
            if (x - 2) not in S:
                out.append((rest | {x - 1}, np.conj(w)))
            # NNN hops
            if (x + 2) not in S and (x + 3) not in S:
                out.append((rest | {x + 2}, -1.0))
            if (x - 2) not in S and (x - 3) not in S:
                out.append((rest | {x - 2}, -1.0))
            # annihilation (dN = -1)
            if N > 1 and self.r != 0:
                amp = np.conj(self.r) * (((x + 2) in S) - ((x - 2) in S))
                if amp != 0:
                    out.append((rest, amp))
        # creation (dN = +1)
        if N < Nmax and self.r != 0:
            cand = set()
            for x in c:
                cand.add(x + 2); cand.add(x - 2)
            for j in cand:
                if j in S or (j - 1) in S or (j + 1) in S:
                    continue
                amp = self.r * (((j + 2) in S) - ((j - 2) in S))
                if amp != 0:
                    out.append((S | {j}, amp))
        return dia, out

    def decomp(self, Nmax, R, Nmin=1):
        """H(K) = sum_s e^{-iKs} M_s  + diag(dia0) - mu*diag(N); cached per (Nmax, R, Nmin, w, r, V2, V3)."""
        key = (Nmax, R, Nmin, self.w, self.r, self.V2, self.V3)
        if getattr(self, "_cache", None) is not None and self._cache[0] == key:
            return self._cache[1]
        basis = []
        for N in range(Nmin, Nmax + 1):
            basis += configs(N, R)
        index = {c: i for i, c in enumerate(basis)}
        mu_save, self.mu = self.mu, 0.0
        dia0 = np.zeros(len(basis)); Nd = np.zeros(len(basis))
        trip = {}
        for i, c in enumerate(basis):
            dia, out = self.actions(c, Nmax)
            dia0[i] = dia; Nd[i] = len(c)
            for newS, amp in out:
                if len(newS) < Nmin:
                    continue
                nc = sorted(newS)
                s = nc[0]
                key2 = tuple(x - s for x in nc)
                if key2[-1] > R:
                    continue
                j = index.get(key2)
                if j is None:
                    continue
                t = trip.setdefault(s, ([], [], []))
                t[0].append(j); t[1].append(i); t[2].append(amp)
        self.mu = mu_save
        n = len(basis)
        mats = {s: sp.csr_matrix((np.array(v, complex), (r_, c_)), shape=(n, n)) for s, (r_, c_, v) in trip.items()}
        out = (mats, dia0, Nd, basis)
        self._cache = (key, out)
        return out

    def build(self, K, Nmax, R, Nmin=1):
        mats, dia0, Nd, basis = self.decomp(Nmax, R, Nmin)
        H = sp.diags(dia0 - self.mu * Nd).astype(complex).tocsr()
        for s, M in mats.items():
            H = H + np.exp(-1j * K * s) * M
        return H, basis

    def lowest(self, K, Nmax, R, k=1, Nmin=1, check=False):
        H, basis = self.build(K, Nmax, R, Nmin)
        if check:
            assert abs(H - H.getH()).max() < 1e-12
        n = H.shape[0]
        if n < 400:
            e = np.linalg.eigvalsh(H.toarray())[:k]
        else:
            e = np.sort(spla.eigsh(H, k=k, which='SA', tol=1e-11, ncv=max(20, 2 * k + 10))[0])
        return e


def bare_valleys(m):
    ks = np.linspace(-np.pi, np.pi, 20001, endpoint=False)
    e = m.eps(ks)
    idx = np.nonzero((e < np.roll(e, 1)) & (e < np.roll(e, -1)))[0]
    # label by position: valley a has cos k > 0 (near 0), valley b has cos k < 0 (near pi)
    return sorted([ks[i] for i in idx], key=lambda k: -np.cos(k))


def dressed_min(m, Nmax, R1, k0, half=0.4):
    f = lambda k: m.lowest(k, Nmax, R1)[0]
    res = minimize_scalar(f, bracket=None, bounds=(k0 - half, k0 + half), method='bounded',
                          options={'xatol': 1e-6})
    return res.x, res.fun


def threshold(m, Nmax, R1):
    ka0, kb0 = bare_valleys(m)[:2]
    # the window must exclude total momenta near 2k_a, where two undressed particles are as cheap
    half = 0.3 * min(abs(np.angle(np.exp(1j * ka0))), abs(np.angle(np.exp(1j * kb0))),
                     abs(np.angle(np.exp(1j * (ka0 - kb0)))) / 2)
    ka, ea = dressed_min(m, Nmax, R1, ka0, half)
    kb, eb = dressed_min(m, Nmax, R1, kb0, half)
    return (ka, ea), (kb, eb)


def tune_mu(m, Nmax, R1, target=0.0):
    """mu such that the lower dressed valley sits at energy `target` (dilute onset)."""
    def f(mu):
        m.mu = mu
        (ka, ea), (kb, eb) = threshold(m, Nmax, R1)
        return min(ea, eb) - target
    mu0 = float(min(m.eps(np.linspace(-np.pi, np.pi, 20001))))
    a, b = mu0 - 1.5, mu0 + 0.5
    mu = brentq(f, a, b, xtol=1e-9)
    m.mu = mu
    return mu


def tune_phi(m, Nmax, R1, phi0=np.pi/2, span=0.08):
    """retune the hop phase so that the two dressed valleys are degenerate (mu retuned each time)."""
    u = abs(m.w)
    def f(phi):
        m.w = u * np.exp(1j * phi)
        tune_mu(m, Nmax, R1)
        (ka, ea), (kb, eb) = threshold(m, Nmax, R1)
        return ea - eb
    fa, fb = f(phi0 - span), f(phi0 + span)
    if fa * fb > 0:
        raise RuntimeError(f"no valley degeneracy in phi window: {fa} {fb}")
    phi = brentq(f, phi0 - span, phi0 + span, xtol=1e-7)
    m.w = u * np.exp(1j * phi)
    tune_mu(m, Nmax, R1)
    return phi


def binding(u, phi, v, Nmax, R, R1=60, r=None, g=2/3, verbose=True, degenerate=False):
    m = Model(u, phi, v, g=g, r=r)
    if degenerate:
        phi = tune_phi(m, Nmax, R1, phi0=phi)
    mu = tune_mu(m, Nmax, R1)
    (ka, ea), (kb, eb) = threshold(m, Nmax, R1)
    K = ka + kb
    E2 = m.lowest(K, Nmax, R)[0]
    B = ea + eb - E2
    out = dict(u=u, phi=phi, v=v, r=m.r, g=g, Nmax=Nmax, R=R, R1=R1, mu=mu, ka=ka, kb=kb, ea=ea, eb=eb,
               K=K, E2=float(E2), B=float(B))
    if verbose:
        print(json.dumps(out), flush=True)
    return out


if __name__ == "__main__":
    pass
