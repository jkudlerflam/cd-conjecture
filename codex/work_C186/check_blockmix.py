"""C186: blockmix.py block Hamiltonian (Lb = 4, 12 sites, penalty U) against a real-space open-chain P-projected build
of the mixed model, and dmrg_d against dense ED of the block Hamiltonian."""
import numpy as np, itertools
from blockmix import block_model, LOC
from dmrg_d import dmrg
u, phi, v, mu, r = 1.0, 1.55, 0.3, -2.3, 0.55
w = u*np.exp(1j*phi)
Lb = 4; n = 3*Lb
def full_block(U):
    hb, hs = block_model(w, V2=-2+v, mu=mu, U=U, r=r)
    d = 5; D = d**Lb
    H = np.zeros((D, D), complex)
    I = np.eye(d)
    for j in range(Lb):
        ops = [I]*Lb; ops[j] = hs
        M = ops[0]
        for o in ops[1:]: M = np.kron(M, o)
        H += M
    for j in range(Lb-1):
        M = np.eye(d**j) if j > 0 else np.eye(1)
        M = np.kron(np.kron(M, hb), np.eye(d**(Lb-j-2)))
        H += M
    return H
# real space, P exact, open chain
S = [s for s in itertools.product((0,1), repeat=n) if all(not (s[i] and s[i+1]) for i in range(n-1))]
pos = {s: i for i, s in enumerate(S)}
H = np.zeros((len(S), len(S)), complex)
for a, s in enumerate(S):
    o = lambda j: s[j] if 0 <= j < n else 0
    for j in range(n):
        if s[j]:
            H[a, a] += -mu + (-2+v)*o(j+2) - (2/3)*o(j+3)
            for dj, am in ((1, w), (-1, np.conj(w))):
                t = j+dj
                if 0 <= t < n and not o(t) and not o(t+dj):
                    x = list(s); x[j] = 0; x[t] = 1; H[pos[tuple(x)], a] += am
            for dj in (2, -2):
                t = j+dj; sg = 1 if dj > 0 else -1
                if 0 <= t < n and not o(t) and not o(t+sg):
                    x = list(s); x[j] = 0; x[t] = 1; H[pos[tuple(x)], a] += -1
            amp = r*(o(j+2) - o(j-2))
            if amp: x = list(s); x[j] = 0; H[pos[tuple(x)], a] += amp
        elif not o(j-1) and not o(j+1):
            amp = r*(o(j+2) - o(j-2))
            if amp: x = list(s); x[j] = 1; H[pos[tuple(x)], a] += amp
assert np.abs(H - H.conj().T).max() < 1e-12
e_rs = np.linalg.eigvalsh(H)[:8]
for U in (40.0, 400.0, 4000.0):
    e_b = np.linalg.eigvalsh(full_block(U))[:8]
    print(f"U={U:7.0f}: max |E_block - E_realspace| over 8 lowest = {np.abs(e_b - e_rs).max():.2e}")
Hb = full_block(40.0)
e0 = np.linalg.eigvalsh(Hb)[0]
hb, hs = block_model(w, V2=-2+v, mu=mu, U=40.0, r=r)
out = dmrg(hb, Lb, chis=(16, 32, 64), cutoff=1e-12, onsite=[hs]*Lb)
print(f"DMRG E = {out['E']:.12f}, dense ED E0 = {e0:.12f}, diff = {abs(out['E']-e0):.1e}")
