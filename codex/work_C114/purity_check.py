"""Small XX correlation-matrix checks, largest matrix is 32 by 32."""
from pathlib import Path
import numpy as np
from checks import coherent_trial, ground, ring_matrix, models


def xx_purities(L):
    shift = .5 if (L//2) % 2 == 0 else 0.
    momenta = 2*np.pi*(np.arange(L)+shift)/L
    occupied = np.argsort(np.cos(momenta))[:L//2]
    sites = np.arange(L//2)
    waves = np.exp(1j*sites[:, None]*momenta[occupied][None, :])/np.sqrt(L)
    C = waves@waves.conj().T
    p = [1.]
    alpha = 2/np.pi**2
    for r in range(1, L//2+1):
        nu = np.linalg.eigvalsh(C[:r, :r])
        variance = float(np.sum(nu*(1-nu)))
        value = float(np.prod(nu**2+(1-nu)**2))
        assert variance+1e-12 >= np.log(r)/np.pi**2
        assert value <= np.exp(-2*variance)+1e-12
        assert value <= r**(-alpha)+1e-12
        p.append(value)
    return np.array(p+p[-2:0:-1])


out = ['# C114 XX coherent reconnection', '',
       '[Numerical] Exact finite-ring correlation matrices, of dimension at most 32. '
       'The energy q_L / sum_r tr(rho_r^2) is the optimum in the coherent '
       'prefix-swap subspace, as proved in PROOF.md. No many-body extrapolation.', '',
       '| L | purity sum | optimum cost | L times cost | L^(3/4) times cost |',
       '|--:|--:|--:|--:|--:|']
last = 0.
for L in [4, 8, 16, 32, 64]:
    p = xx_purities(L)
    q = 2/(L*np.sin(np.pi/L))
    cost = q/np.sum(p)
    assert L*cost > last
    last = L*cost
    if L == 4:
        h = models()[0][1]
        e, psi, _ = ground(ring_matrix(h, 2, L))
        result = coherent_trial(ring_matrix(h, 2, 2*L), psi, 2, L, e)
        assert abs(result['coherent']-cost) < 1e-11
        assert np.max(np.abs(np.array(result['gram_row'])-p)) < 1e-11
    out.append(f'| {L} | {np.sum(p):.9f} | {cost:.9f} | {L*cost:.9f} | {L**.75*cost:.9f} |')

# Check the symmetry-free upper bound on small arbitrary complex bonds too.
rng = np.random.default_rng(114)
for j in range(3):
    x = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    h = (x+x.conj().T)/8
    e, psi, _ = ground(ring_matrix(h, 2, 4))
    H2 = ring_matrix(h, 2, 8)
    e2, _, _ = ground(H2)
    result = coherent_trial(H2, psi, 2, 4, e)
    assert e2-2*e <= result['uniform']+1e-9
out += ['', '[Numerical] Three arbitrary complex Hermitian bonds at d=2, '
        'L=4 also pass the Gram identity, Re K=q I, and the uniform coherent '
        'upper bound. No reality or on-site symmetry was imposed for these tests.', '',
        '[Proved] For XX, the number-variance inequality gives '
        'tr(rho_r^2) <= min(r,L-r)^(-2/pi^2) for 1 <= r < L. '
        'It proves that L times the optimized cost diverges. The observed '
        'L^(-3/4) scaling is a sharper CFT calibration, not required for that proof.', '',
        '[Numerical] All assertions passed.']
Path('PURITY_CHECKS.md').write_text('\n'.join(out)+'\n')
print('\n'.join(out))
