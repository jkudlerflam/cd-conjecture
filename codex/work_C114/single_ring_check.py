"""Evaluate the proved coherent cost using only small single-ring ED."""
from pathlib import Path
import numpy as np
from checks import models, ground, ring_matrix, purity

out = ['# C114 additional single-ring checks', '',
       '[Numerical] The coherent cost is evaluated using the proved formula '
       'q_L/S_L. This does not diagonalize the doubled ring. All single-ring '
       'Hilbert spaces have dimension at most 729.', '',
       '| model | L | q_L | S_L | q_L/S_L | L q_L/(M S_L) |',
       '|:--|--:|--:|--:|--:|--:|']
max_residual = 0.
for name, h, d, M, A in models():
    for L in ([4, 6, 8] if d == 2 else [3, 6]):
        e, psi, residual = ground(ring_matrix(h, d, L))
        max_residual = max(max_residual, residual)
        mat = psi.reshape(-1, d)
        sigma = mat.T@mat.conj()
        b = float(np.trace(h@np.kron(sigma, sigma)).real)
        q = 2*(b-e/L)
        ps = [purity(psi, d, L, r) for r in range(L//2+1)]
        S = 1+2*sum(ps[1:])-(ps[-1] if L % 2 == 0 else 0)
        assert q > 0
        assert q <= 4*M+1e-10
        assert S >= sum(d**(-min(r, L-r)) for r in range(L))-1e-10
        out.append(f'| {name} | {L} | {q:.9f} | {S:.9f} | {q/S:.9f} | {L*q/(M*S):.9f} |')
out += ['', f'[Numerical] Maximum eigenpair residual {max_residual:.3g}. '
        'All dimension and norm inequalities passed. These sizes do not '
        'establish asymptotic powers for the interacting models.']
Path('SINGLE_RING_CHECKS.md').write_text('\n'.join(out)+'\n')
print('\n'.join(out))
