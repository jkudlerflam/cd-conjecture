"""Commuting nearest-neighbour chains in Bravyi-Vyalyi form: R_L <= log m.

Site space C^3 = (C^1 x C^1) + (C^1 x C^2): sector alpha=0 has left/right
factor dims (1,1), sector alpha=1 has (1,2). The bond acts on the right factor
of site i and the left factor of site i+1, blockwise in (alpha_i, alpha_{i+1}),
with random Hermitian blocks. Checks: bonds sharing a site commute, the
transfer-matrix formula Z_L = Tr T^L, and R_L <= log m (m = 2 sectors).
Output: commuting_results.json."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from zeno_run import _embed

ROOT = Path(__file__).resolve().parent
rng = np.random.default_rng(7)
# basis of a site: index 0 -> sector 0 (r-dim 1); indices 1,2 -> sector 1 with r = 0,1
sec = [0, 1, 1]; rpos = [0, 0, 1]; rdim = [1, 2]
blocks = {}
for a in range(2):
    for b in range(2):
        D = rdim[a]                       # left dims are all 1
        X = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
        blocks[(a, b)] = ((X + X.conj().T) / 2).real if D == 1 else (X + X.conj().T) / 2
h = np.zeros((9, 9), complex)
for i in range(3):
    for j in range(3):
        for i2 in range(3):
            for j2 in range(3):
                if sec[i] != sec[i2] or sec[j] != sec[j2] or j != j2: continue
                h[3 * i2 + j2, 3 * i + j] = np.atleast_2d(blocks[(sec[i], sec[j])])[rpos[i2], rpos[i]]
assert np.allclose(h, h.conj().T)

def embed_c(op, i, j, n):
    return _embed(op.real.reshape(3, 3, 3, 3), i, j, n, 3) + 1j * _embed(op.imag.reshape(3, 3, 3, 3), i, j, n, 3)

h12, h23 = embed_c(h, 0, 1, 3), embed_c(h, 1, 2, 3)
comm = float(np.abs(h12 @ h23 - h23 @ h12).max())
assert comm < 1e-12

_spec = {}
def logZ(n, beta):
    if n not in _spec:
        H = sum(embed_c(h, i, (i + 1) % n, n) for i in range(n))
        _spec[n] = eigh(H, eigvals_only=True)
    ev = _spec[n]
    return float(np.log(np.sum(np.exp(-beta * (ev - ev[0])))) - beta * ev[0])

def logZ_transfer(n, beta):
    T = np.array([[np.sum(np.exp(-beta * np.atleast_1d(np.linalg.eigvalsh(np.atleast_2d(blocks[(a, b)])))))
                   for b in range(2)] for a in range(2)])
    return float(np.log(np.trace(np.linalg.matrix_power(T, n)).real))

rows = []
for L in [2, 3, 4]:
    for beta in [0.5, 2., 8., 32.]:
        lzL, lzD = logZ(L, beta), logZ(2 * L, beta)
        R = 2 * lzL - lzD
        err = abs(lzL - logZ_transfer(L, beta))
        assert err < 1e-9 and R <= np.log(2) + 1e-12
        rows.append(dict(L=L, beta=beta, R=R, transfer_error=err))
        print(L, beta, 'R', round(R, 6), 'log 2', round(np.log(2), 6), 'transfer err', err)
(ROOT / 'commuting_results.json').write_text(json.dumps(dict(commutator=comm, rows=rows), indent=1) + '\n')
