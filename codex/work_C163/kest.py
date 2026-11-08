"""C163 (Claude Opus agent, 2026-11-08): Luttinger-parameter estimator for half-floating liquids on rings.
K_S(L) = (L/4) S_mm(2 pi/L), S_mm(q) = sum_r e^{iqr} <s^m_0 s^m_r>_conn, m = spiral axis of the polar decomposition
(the Theta-odd density that carries the parent U(1) charge).  For a Luttinger liquid with the commensurate vertex
e^{2 i phi} of dimension K, s^m ~ (2/pi) d_x phi gives S_mm(q) -> 2 K |q| / pi, so K_S -> K (XX: K_S = 1 exactly
at every L, since S_zz(q) = 2|q|/pi for |q| < pi in the free band).
usage: python3 kest.py  (self test)"""
import numpy as np
import est as E, tk


def k_est(v, L=16):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    J, f = E.J_f(v)
    m = tk.spiral_parent(J, f)['m']
    h = E.bond_from_JB(J, f)
    d = E.ring_data(h, L, nev=1)
    psi = d['psi']
    mag = tk.magnetization(psi, L)
    Smm = 0.0
    q = 2 * np.pi / L
    Cr = []
    for r in range(L):
        if r == 0:
            c = 1.0
        else:
            C = tk.corr_matrix(psi, L, r)
            c = float(np.real(m @ C @ m))
        Cr.append(c)
    Cr = np.array(Cr) - (m @ mag / L) ** 2
    Smm = float(np.real(np.sum(np.exp(1j * q * np.arange(L)) * Cr)))
    return L * Smm / 4, dict(m=m.tolist(), mag=(mag / L).tolist(), Cr=Cr.tolist(), m0=d['m0'])


if __name__ == '__main__':
    import time
    t = time.time()
    for De in (0.0, -0.5, -0.7, -0.9):
        v = np.array([1, 1, De, 0, 0, 0, 0, 0, 0.0])
        K12 = k_est(v, 12)[0]; K16 = k_est(v, 16)[0]
        print(f'XXZ Delta={De}: K_S(12)={K12:.3f} K_S(16)={K16:.3f}  exact K={tk.luttinger_K(De):.3f}', flush=True)
    print('time', time.time() - t)
