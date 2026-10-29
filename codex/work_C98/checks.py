"""Small local checks for C98. NumPy/SciPy only; no network or external data."""
import math
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import digamma, gammaln
from scipy.integrate import quad


def ring_matrix(h, d, length, periodic=True):
    """Explicit sparse tensor-product Hamiltonian, used only at small sizes."""
    size = d**length
    indices = np.arange(size, dtype=np.int64)
    powers = d ** np.arange(length, dtype=np.int64)
    digits = (indices[:, None] // powers) % d
    rows, cols, vals = [], [], []
    for i in range(length if periodic else length-1):
        j = (i+1) % length
        local = digits[:, i] * d + digits[:, j]
        for a, b in zip(*np.nonzero(np.abs(h) > 1e-14)):
            chosen = indices[local == b]
            target = chosen + (a//d-b//d)*powers[i] + (a%d-b%d)*powers[j]
            rows.append(target)
            cols.append(chosen)
            vals.append(np.full(chosen.size, h[a, b], dtype=complex))
    return sparse.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(size, size)).tocsr()


def ground(h, d, length, periodic=True):
    matrix = ring_matrix(h, d, length, periodic)
    # Deterministic generic vector avoids special symmetry subspaces.
    v0 = np.random.default_rng(98+length).normal(size=matrix.shape[0])
    e, vec = eigsh(matrix, k=1, which='SA', tol=2e-12, v0=v0)
    residual = np.linalg.norm(matrix @ vec[:, 0] - e[0]*vec[:, 0])
    assert residual < 2e-9, residual
    return float(e[0]), residual


def xxz_bulk(delta):
    if abs(delta-1) < 1e-14:
        return 0.25-math.log(2)
    gamma = math.acos(delta)
    # Stable form of sinh((pi-gamma)x)/(sinh(pi x) cosh(gamma x)).
    def integrand(x):
        if x == 0:
            return (math.pi-gamma)/math.pi
        return (2*math.exp(-2*gamma*x) * (-math.expm1(-2*(math.pi-gamma)*x))
                / ((-math.expm1(-2*math.pi*x))*(1+math.exp(-2*gamma*x))))
    integral, err = quad(integrand, 0, np.inf, epsabs=1e-12, epsrel=1e-12)
    assert err < 1e-9
    return delta/4-math.sin(gamma)*integral


def embed_pair(op, i, j, n=4, d=2):
    size = d**n
    mat = np.zeros((size, size), dtype=complex)
    for col in range(size):
        digits = [(col // d**k) % d for k in range(n)]
        b = digits[i]*d+digits[j]
        for a in range(d*d):
            row = col+(a//d-digits[i])*d**i+(a%d-digits[j])*d**j
            mat[row, col] += op[a, b]
    return mat


def main():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)/2
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)/2
    sz = np.diag([1., -1.])/2
    xx = np.kron(sx, sx)+np.kron(sy, sy)
    zz = np.kron(sz, sz)
    out = ['# C98 checks', '',
           '[Numerical] Small-ring ED, exact-formula evaluations and finite-matrix identities. '
           'No large-size interacting extrapolation or new Bethe-ansatz computation.', '',
           '## XXZ ring calibration', '',
           'h=SxSx+SySy+Delta SzSz, M=1/2, d=2. Q_L=6L(e0 L-E_L)/(pi^2 M).', '',
           '| Delta | L | E_L | Q_L | residual |', '|:--|--:|--:|--:|--:|']
    energies = {}
    max_residual = 0.
    for delta in [-.5, 0., .5, 1.]:
        bulk = xxz_bulk(delta)
        for length in [4, 6, 8, 10, 12]:
            e, residual = ground(xx+delta*zz, 2, length)
            energies[delta, length] = e
            max_residual = max(max_residual, residual)
            q = 6*length*(bulk*length-e)/(math.pi**2*.5)
            out.append(f'| {delta:g} | {length} | {e:.12f} | {q:.9f} | {residual:.2g} |')
            if delta == 0:
                assert abs(e+1/math.sin(math.pi/length)) < 1e-10
        if delta == 0:
            assert abs(bulk+1/math.pi) < 1e-12
    out += ['', '## Finite sharp-bound failures', '']
    e4 = energies[1., 4]
    deficit = 4*xxz_bulk(1.)-e4
    sharp = math.pi**2/(12*4)
    doubling = energies[1., 8]-2*e4
    double_bound = math.pi**2/(8*4)
    assert deficit > sharp and doubling > double_bound
    out += [f'Heisenberg L=4 deficit {deficit:.12f}, proposed no-remainder ceiling {sharp:.12f}.',
            f'Heisenberg E(8)-2E(4)={doubling:.12f}, proposed no-remainder ceiling {double_bound:.12f}.',
            '', '## XX open-chain and rank checks', '',
            '| L | E_open exact | ED error | half-chain Schmidt rank |', '|--:|--:|--:|--:|']
    for length in [4, 6, 8, 10, 12]:
        op = .5-.5/math.sin(math.pi/(2*(length+1)))
        e, _ = ground(xx, 2, length, periodic=False)
        assert abs(e-op) < 1e-10
        # Half-filled spin-ring JW parity selects periodic/AP momenta.
        shift = .5 if (length//2) % 2 == 0 else 0.
        mom = 2*math.pi*(np.arange(length)+shift)/length
        occ = np.argsort(np.cos(mom))[:length//2]
        fourier = np.exp(1j*np.arange(length)[:, None]*mom[None, :])/math.sqrt(length)
        block = fourier[:length//2, occ]
        eigen = np.linalg.eigvalsh(block @ block.conj().T)
        assert eigen.min() > 1e-12 and eigen.max() < 1-1e-12
        out.append(f'| {length} | {op:.12f} | {abs(e-op):.2g} | {2**(length//2)} |')
    out += ['', f'Open-minus-ring limiting cost = 1/2-1/pi = {.5-1/math.pi:.12f}.',
            f'Product-of-two-ring seam cost tends to 2/pi = {2/math.pi:.12f}.',
            '', '## Fundamental SU(N) ring calibration', '',
            'h=P, M=1. Vacuum-compatible sizes only. Q_L uses the known Bethe bulk energy.', '',
            '| N | L | E_L | Q_L | Q_infinity |', '|--:|--:|--:|--:|--:|']
    for n, lengths in [(2, [4, 6, 8]), (3, [3, 6, 9]), (4, [4, 8])]:
        p = np.zeros((n*n, n*n))
        for a in range(n):
            for b in range(n):
                p[a*n+b, b*n+a] = 1
        bulk = 1+2*(digamma(1/n)+np.euler_gamma)/n
        qinf = 2*(n-1)/(n*math.log2(n))
        for length in lengths:
            e, residual = ground(p, n, length)
            max_residual = max(max_residual, residual)
            q = 6*length*(bulk*length-e)/(math.pi**2*math.log2(n))
            out.append(f'| {n} | {length} | {e:.12f} | {q:.9f} | {qinf:.9f} |')
    out += ['', '## Symmetric SU(N) inherited continuum formulas', '',
            'Formula checks only. The all-N,k proof and physical inputs are inherited from C14/C35.', '',
            '| N | k | log_2 d | A/(2pi M) | A/(pi M log_2 d) |', '|--:|--:|--:|--:|--:|']
    for n, k in [(2,1), (2,2), (3,1), (3,2), (7,2), (8,2), (10,10), (100,100)]:
        logd = (gammaln(n+k)-gammaln(n)-gammaln(k+1))/math.log(2)
        harmonic = sum(1/j for j in range(1, k+1))
        ratio = k*(n*n-1)/(n*(n+k)*harmonic)
        q = 2*ratio/logd
        assert q < 1+1e-12
        out.append(f'| {n} | {k} | {logd:.9f} | {ratio:.9f} | {q:.9f} |')
    out += ['', '## Uniform partial-copy swap', '',
            'Four-site density-matrix checks use the thermodynamic nearest-neighbor marginal.', '',
            '| model | theta | measured cost per doubled bond | formula | error |',
            '|:--|--:|--:|--:|--:|']
    swap = np.array([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]], dtype=complex)
    for name, h, corr in [
        ('Heisenberg', xx+zz, [(.25-math.log(2))/3]*3),
        ('XX', xx, [-1/(2*math.pi), -1/(2*math.pi), -1/math.pi**2])]:
        rho = np.eye(4, dtype=complex)/4
        for s, r in zip([sx, sy, sz], corr):
            rho += 4*r*np.kron(s, s)
        assert np.linalg.eigvalsh(rho).min() >= -1e-12
        rho4 = embed_pair(rho, 0, 1) @ embed_pair(rho, 2, 3)
        assert abs(np.trace(rho4)-1) < 1e-12
        h4 = embed_pair(h, 0, 1)+embed_pair(h, 2, 3)
        e = np.trace(h @ rho).real
        prefactor = -e+4*e*e/3 if name == 'Heisenberg' else 1/math.pi+4/math.pi**3
        for theta in [0., math.pi/8, math.pi/4, math.pi/2]:
            gate = math.cos(theta)*np.eye(4)+1j*math.sin(theta)*swap
            unitary = embed_pair(gate, 0, 2) @ embed_pair(gate, 1, 3)
            cost = np.trace(h4 @ unitary @ rho4 @ unitary.conj().T).real-2*e
            predicted = math.sin(2*theta)**2*prefactor
            assert abs(cost-predicted) < 1e-12
            out.append(f'| {name} | {theta:.9f} | {cost:.12f} | {predicted:.12f} | {abs(cost-predicted):.2g} |')
    out += ['', '## Dyadic identity and concavity obstruction', '']
    # Model energy sequence, not a Hamiltonian construction.
    for xi in [1., 10., 1000.]:
        delta = lambda length: xi/(length+xi)
        length = 8
        terms = [(2*delta((2.**j)*length)-delta((2.**(j+1))*length))/2.**(j+1)
                 for j in range(45)]
        assert abs(sum(terms)-delta(length)) < 1e-12
        lgrid = np.arange(1., 1001.)
        seq = -delta(lgrid)
        assert np.max(np.diff(seq, n=2)) <= 1e-12
        out.append(f'xi={xi:g}: dyadic identity error {abs(sum(terms)-delta(length)):.2g}; '
                   f'deficit <=1, Casimir-like coefficient K={xi:g}.')
    out += ['', f'Maximum reported ring eigenpair residual {max_residual:.3g}.',
            '', '[Numerical outcome] All assertions passed. The two no-remainder Heisenberg '
            'inequalities are violated as expected. The asymptotic logarithmic conjecture '
            'was not falsified.']
    Path('CHECKS.md').write_text('\n'.join(out)+'\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
