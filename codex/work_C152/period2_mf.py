"""Fermi-point count for period-2 nearest-neighbour quadratic (BdG) chains, the form a CDW / bond-ordered
Hartree-Fock-Bogoliubov state takes in this class. c_free = (number of sign changes of all BdG bands over the
reduced zone)/4. Sanity: the uniform case (period-2 description of a period-1 chain) must give at most 1."""
import numpy as np
rng = np.random.default_rng(7)
K = np.linspace(-np.pi, np.pi, 1201)[:-1]
def bdg(hR, DR):
    # hR, DR: dict R -> 2x2 blocks <a,0| . |b,R>; h(k) = sum_R hR e^{ikR}
    h = sum(m[None] * np.exp(1j * K * R)[:, None, None] for R, m in hR.items())
    hm = sum(m[None] * np.exp(-1j * K * R)[:, None, None] for R, m in hR.items())
    D = sum(m[None] * np.exp(1j * K * R)[:, None, None] for R, m in DR.items())
    top = np.concatenate([h, D], axis=2)
    bot = np.concatenate([np.conj(np.transpose(D, (0, 2, 1))), -np.transpose(hm, (0, 2, 1))], axis=2)
    return np.linalg.eigvalsh(np.concatenate([top, bot], axis=1))
def model(p):
    muA, muB, t1r, t1i, t2r, t2i, d1r, d1i, d2r, d2i = p
    t1 = t1r + 1j * t1i; t2 = t2r + 1j * t2i; d1 = d1r + 1j * d1i; d2 = d2r + 1j * d2i
    # sites A_n, B_n; bonds (A_n, B_n): t1, d1 ; (B_n, A_{n+1}): t2, d2
    h0 = np.array([[-muA, t1], [np.conj(t1), -muB]]); h1 = np.array([[0, 0], [t2, 0]]); hm1 = h1.conj().T
    # pairing Delta_ij antisymmetric: Delta_{A0,B0} = d1, Delta_{B0,A1} = d2
    D0 = np.array([[0, d1], [-d1, 0]]); D1 = np.array([[0, 0], [d2, 0]]); Dm1 = np.array([[0, -d2], [0, 0]])
    return {0: h0, 1: h1, -1: hm1}, {0: D0, 1: D1, -1: Dm1}
def count(p):
    E = bdg(*model(p)); s = np.sign(E)
    return np.sum(s != np.roll(s, 1, axis=0)) / 4.0
if __name__ == "__main__":
    # sanity: uniform chains written with period 2 (muA = muB, t1 = t2, d1 = d2)
    mx = 0
    for i in range(400):
        q = rng.normal(size=5); p = [q[0], q[0], q[1], q[2], q[1], q[2], q[3], q[4], q[3], q[4]]
        mx = max(mx, count(p))
    print('uniform sanity, max c_free =', mx)
    hist = {}; best = 0
    for i in range(20000):
        p = rng.normal(size=10); c = count(p); hist[c] = hist.get(c, 0) + 1
        if c > best: best, bp = c, p
    print('period-2 random:', dict(sorted(hist.items())), 'max', best, np.round(bp, 3))
