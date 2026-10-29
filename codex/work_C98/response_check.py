"""Check the exact two-seam response identity on two three-site XX rings."""
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from checks import embed_pair

sx = np.array([[0, 1], [1, 0]], dtype=complex)/2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)/2
h = np.kron(sx, sx)+np.kron(sy, sy)
length = 3
n = 2*length
base = sum((embed_pair(h, j+offset, (j+1)%length+offset, n=n)
            for offset in [0, length] for j in range(length)),
           np.zeros((2**n, 2**n), dtype=complex))


def switch(j):
    k = (j+1)%length
    return (embed_pair(h, j, k+length, n=n)+embed_pair(h, j+length, k, n=n)
            -embed_pair(h, j, k, n=n)-embed_pair(h, j+length, k+length, n=n))


vx, vy = switch(length-1), switch(0)
corners = [np.linalg.eigvalsh(base+s*vx+t*vy) for s,t in [(0,0),(1,0),(0,1),(1,1)]]
assert np.max(np.abs(corners[0]-corners[3])) < 1e-12
assert np.max(np.abs(corners[1]-corners[2])) < 1e-12


def free(energy, beta):
    return energy[0]-np.log(np.exp(-beta*(energy-energy[0])).sum())/beta


def susceptibility(s, t, beta):
    energy, u = np.linalg.eigh(base+s*vx+t*vy)
    weights = np.exp(-beta*(energy-energy[0]))
    weights /= weights.sum()
    x, y = u.conj().T @ vx @ u, u.conj().T @ vy @ u
    de = energy[None, :]-energy[:, None]
    dp = weights[:, None]-weights[None, :]
    coeff = beta*(weights[:, None]+weights[None, :])/2
    np.divide(dp, de, out=coeff, where=np.abs(de)>1e-10)
    meanx = np.dot(weights, x.diagonal()).real
    meany = np.dot(weights, y.diagonal()).real
    return np.sum(coeff*x*y.T).real-beta*meanx*meany


nodes, weights = leggauss(24)
nodes, weights = (nodes+1)/2, weights/2
out = ['# Two-seam response check', '',
       '[Numerical] Two three-site XX rings, six qubits total. These odd rings test the '
       'finite-temperature matrix identity, not a vacuum Casimir coefficient. '
       'A 24 by 24 Gauss-Legendre quadrature integrates the Kubo susceptibility.', '',
       '| beta | F(2L)-2F(L) | half integrated chi | error | min sampled chi | max sampled chi |',
       '|--:|--:|--:|--:|--:|--:|']
for beta in [2., 8.]:
    data = np.array([[susceptibility(s,t,beta) for t in nodes] for s in nodes])
    integral = .5*np.einsum('i,j,ij', weights, weights, data)
    defect = free(corners[1], beta)-free(corners[0], beta)
    err = abs(integral-defect)
    assert err < 1e-7, err
    out.append(f'| {beta:g} | {defect:.12f} | {integral:.12f} | {err:.2g} | '
               f'{data.min():.9f} | {data.max():.9f} |')
out += ['', '[Numerical outcome] All assertions passed. Susceptibility need not be positive '
        'on non-vacuum-compatible finite rings.']
Path('RESPONSE_CHECK.md').write_text('\n'.join(out)+'\n')
print('\n'.join(out))
