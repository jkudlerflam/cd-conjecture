"""C159: structure-factor test of the 'half-floating' statement (PROOF section 2.3).
For each bond: ring ED ground state, connected correlations of the spin component along the spiral axis m,
structure factor S_mm(k) and the transverse one; prints the peak momenta."""
import sys, time, json
import numpy as np
import est as E
import tk

pts = {
    'T_C68_eta0.1': [1.1, 0.9, 0.4, 0.1, 0.2, 0.5, 0, 0, 0],
    'K_C68_eta0.1': [1.1, 0.9, 0.4, 0, 0, 0.5, 0.1, 0.2, 0],
    'N_C68_eta0.1': [1.1, 0.9, 0.4, 0.1, 0.2, 0.5, 0.3, 0.4, 0.5],
    'U1_XXZ_Dz_hz0.5': [1, 1, 0.4, 0, 0, 0.5, 0, 0, 0.5],
    'U1_XXZ_Dz_h0': [1, 1, 0.4, 0, 0, 0.5, 0, 0, 0],
}
L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
names = sys.argv[2:] if len(sys.argv) > 2 else list(pts)
for name in names:
    t0 = time.time()
    v = np.array(pts[name], float)
    J, f = E.J_f(v)
    par = tk.spiral_parent(J, f)
    m = par['m']
    h = E.bond_from_JB(J, f)
    d = E.ring_data(h, L, nev=2)
    psi = d['psi']
    M = tk.magnetization(psi, L)
    Cs = [tk.corr_matrix(psi, L, r) for r in range(L)]
    ks = 2 * np.pi * np.arange(L) / L
    Smm = []; Sperp = []
    P = np.eye(3) - np.outer(m, m)
    for k in ks:
        Smm.append(sum(np.exp(-1j * k * r) * (m @ Cs[r] @ m - (m @ M) ** 2) for r in range(L)).real)
        Sperp.append(sum(np.exp(-1j * k * r) * (np.trace(P @ Cs[r] @ P) - M @ P @ M) for r in range(L)).real)
    Smm = np.array(Smm); Sperp = np.array(Sperp)
    lev = np.sort(d['lev'].ravel()); gap = lev[1] - lev[0]
    print(json.dumps(dict(name=name, L=L, theta_over_pi=par['theta'] / np.pi, m=np.round(m, 3).tolist(), Delta=par['Delta'],
                          harm=par['harm'], hm=par['hm'], M=np.round(M, 4).tolist(), m0=d['m0'], Lgap=L * gap,
                          k_peak_mm_over_pi=ks[np.argmax(Smm)] / np.pi, Smm=np.round(Smm, 3).tolist(),
                          k_peak_perp_over_pi=ks[np.argmax(Sperp)] / np.pi, time=time.time() - t0)), flush=True)
