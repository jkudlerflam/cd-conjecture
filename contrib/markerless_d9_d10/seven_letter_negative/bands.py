"""Band-structure diagnostics for a logical Bloch matrix h(k) = mu + T e^{ik} + T^dag e^{-ik}."""
import numpy as np, scipy.linalg as sl, itertools
def fermi_roots(mu, T, E=0.0):
    """all 2n roots z of det(mu - E + T z + T^dag/z) = 0 (scale-invariant linearisation)"""
    n = mu.shape[0]; sc = np.abs(T).max()
    Ts = T/sc; M = (mu - E*np.eye(n))/sc
    A = np.block([[np.zeros((n,n)), np.eye(n)],[-Ts.conj().T, -M]])
    B = np.block([[np.eye(n), np.zeros((n,n))],[np.zeros((n,n)), Ts]])
    return sl.eigvals(A, B)
def margins(ks):
    """min separation of Fermi points (mod 2pi) and of distinct unordered pair sums (mod 2pi)"""
    ks = np.sort(np.mod(ks, 2*np.pi)); m = len(ks)
    d2 = min((ks[(i+1)%m]-ks[i])%(2*np.pi) for i in range(m))
    sums = sorted((ks[i]+ks[j])%(2*np.pi) for i in range(m) for j in range(i, m))
    d4 = min((sums[(i+1)%len(sums)]-sums[i])%(2*np.pi) for i in range(len(sums)))
    return d2, d4
def report(mu, T, E=0.0, tol=1e-7):
    n = mu.shape[0]; z = fermi_roots(mu, T, E)
    on = np.abs(np.abs(z)-1) < tol; ks = np.angle(z[on])
    out = dict(n=n, on_circle=int(on.sum()), max_off=float(np.abs(np.abs(z)-1).max()))
    if on.sum() >= 2:
        out['d2'], out['d4'] = margins(ks)
        out['reciprocal'] = bool(np.allclose(np.sort(np.mod(ks,2*np.pi)), np.sort(np.mod(-ks,2*np.pi)), atol=1e-6))
        out['pi_shift'] = bool(np.allclose(np.sort(np.mod(ks,2*np.pi)), np.sort(np.mod(ks+np.pi,2*np.pi)), atol=1e-6))
    out['ks'] = np.sort(ks)
    return out
