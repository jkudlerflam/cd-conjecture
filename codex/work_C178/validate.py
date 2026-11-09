"""Validation of kring.py against C173's full-basis ring ED (ring_ed.py, copied logic via import path)."""
import sys, numpy as np, scipy.linalg as sla
sys.path.insert(0, '../work_C173')
import ring_ed as R
import kring as Kr

def full_spec(L, N, w, theta=0.0, **kw):
    out = []
    for m in range(L):
        K = 2*np.pi*m/L
        H, reps, P, keep = Kr.hamiltonian(L, N, K, w, theta=theta, **kw)
        if H.shape[0]:
            out += list(np.linalg.eigvalsh(H.toarray()))
    return np.sort(out)

w = np.exp(1j*np.pi/2)
for (L, N) in ((14, 3), (16, 4), (15, 4)):
    for theta, bc in ((0.0, (-1)**(N-1)), (np.pi, -(-1)**(N-1))):
        a = full_spec(L, N, w, theta=theta)
        b, H = R.build_ring(L, N, w, bc=bc)
        e = np.linalg.eigvalsh(H.toarray())
        print(f"L={L} N={N} theta={theta:.3f}: dim {len(a)} vs {len(e)}, max|diff| = {np.abs(a-e).max():.2e}")
# phi generic
w2 = 1.0*np.exp(1j*1.3)
a = full_spec(16, 4, w2); b, H = R.build_ring(16, 4, w2, bc=-1); e = np.linalg.eigvalsh(H.toarray())
print("phi=1.3 L=16 N=4 maxdiff", np.abs(a-e).max())
# controls: no blockade, NNN hop w=0 -> boson legs; compare to C173 ring_ed with blockade False (fermion signs) only in spectrum count
a = full_spec(12, 3, 0.0, blockade=False, V2=0, V3=0)
print("two-leg control L=12 N=3 dim", len(a), "E0", a[0])
# n(k) check against C173 obdm at L=16, N=4, twist 0 (bosons periodic = C173 bc=-1)
L, N = 16, 4
Kmin, emin, best = None, 1e9, None
for m in range(L):
    K = 2*np.pi*m/L
    H, reps, P, keep = Kr.hamiltonian(L, N, K, w)
    e, v = np.linalg.eigh(H.toarray())
    if e[0] < emin: emin, Kmin, best = e[0], K, (v[:, 0], reps, P, keep)
psi, reps, P, keep = best
ks, nk = Kr.momentum_distribution(L, N, Kmin, psi, reps, P, keep)
b, H = R.build_ring(L, N, w, bc=-1)
ee, vv = np.linalg.eigh(H.toarray())
print("GS energies", emin, ee[0], ee[1])
if abs(ee[1]-ee[0]) > 1e-8:
    G = R.obdm(b, vv[:, 0], L, N, bc=-1); ks2, nk2 = R.nk(G, L, -1, N)
    print("nk maxdiff", np.abs(np.sort(nk)-np.sort(nk2)).max(), " sum", nk.sum())
print("n(k):", np.round(nk, 3))
# S check: at self-conjugate K the S matrix commutes with H and is unitary, S^2 = 1
for (L, N) in ((16, 4), (18, 3), (20, 4)):
    for m in range(L):
        K = 2*np.pi*m/L
        if abs(np.exp(1j*(2*K - N*np.pi)) - 1) > 1e-9: continue
        H, reps, P, keep = Kr.hamiltonian(L, N, K, w)
        S = Kr.S_matrix(L, N, K, reps, P, keep).toarray(); Hd = H.toarray()
        print(f"L={L} N={N} K/pi={K/np.pi:.3f}: |[H,S]|={np.abs(Hd@S-S@Hd).max():.1e} |S^2-1|={np.abs(S@S-np.eye(len(S))).max():.1e} |S S^+ -1|={np.abs(S@S.conj().T-np.eye(len(S))).max():.1e}")
# n(k) check with a nondegenerate ground state
for (L, N, phi) in ((18, 3, 1.3), (16, 4, 1.1), (20, 5, 0.9)):
    w = np.exp(1j*phi)
    b, H = R.build_ring(L, N, w, bc=(-1)**(N-1))
    ee, vv = sla.eigh(H.toarray())
    G = R.obdm(b, vv[:, 0], L, N, bc=(-1)**(N-1)); ks2, nk2 = R.nk(G, L, 1, N)  # C173's nk picks the half-integer grid for bc*(-1)^(N-1) = -1; pass bc=1 to get e^{ikL}=(-1)^(N-1)
    best = None
    for m in range(L):
        K = 2*np.pi*m/L
        Hk, reps, P, keep = Kr.hamiltonian(L, N, K, w)
        if Hk.shape[0] == 0: continue
        e, v = np.linalg.eigh(Hk.toarray())
        if best is None or e[0] < best[0]: best = (e[0], K, v[:, 0], reps, P, keep)
    ks, nk = Kr.momentum_distribution(L, N, best[1], best[2], *best[3:])
    print(f"nk range [{nk.min():.3f},{nk.max():.3f}]; nk check L={L} N={N} phi={phi}: gap {ee[1]-ee[0]:.3e}, E0 {best[0]:.10f} vs {ee[0]:.10f}, max|nk diff| = {np.abs(nk-nk2).max():.2e}, kgrid diff {np.abs(ks-ks2).max():.1e}")
