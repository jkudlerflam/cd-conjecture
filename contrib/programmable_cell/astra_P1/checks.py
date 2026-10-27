"""Small independent consistency checks of the analytic constructions."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.linalg import eigh, null_space

rng = np.random.default_rng(9182)
omega = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
Pomega = np.outer(omega, omega.conj())
Qomega = np.eye(4) - Pomega
norm = lambda a: float(np.linalg.norm(a, 2))
lines = []

# Signed-penalty counterexample, one syndrome, T=1, target gap 1.
for J in (1, 8, 10):
    hs = J/2 * np.array([[1., -1.], [-1., 1.]]) + np.diag([4., -2.])
    lines.append(f'signed example J={J}: eigenvalues={np.linalg.eigvalsh(hs)}')

# General positive ancillary dilation of a random target, including Omega/syndrome mixing.
D, dc, N = 4, 5, 2
Gw = np.eye(D, dtype=complex)[:, :N]
G = np.kron(omega[:, None], Gw)
Qc = np.eye(4*D) - G@G.conj().T
x = rng.normal(size=(4*D, 4*D)) + 1j*rng.normal(size=(4*D, 4*D))
A = Qc @ (x@x.conj().T / (4*D) + np.eye(4*D)) @ Qc
values, vectors = eigh(A)
E = np.zeros((dc*D, D), complex)
E[:D, :] = np.eye(D)
Eps = np.kron(np.eye(4), E)
result = np.zeros_like(A)
max_kernel_residual = 0.
max_unitarity_error = 0.
for val, vec in zip(values, vectors.T):
    if val < 1e-8:
        continue
    # For coefficients vec[ps,w], np.linalg.svd gives workspace kets Vh rows.
    a, singular, vh = np.linalg.svd(vec.reshape(4, D), full_matrices=False)
    rank = sum(singular > 1e-12)
    ws = vh[:rank, :].T
    ws = np.column_stack([ws, null_space(ws.conj().T)])
    source_basis = E@ws
    destinations = np.zeros_like(source_basis)
    for l in range(rank):
        destinations[(l+1)*D, l] = 1
    for k in range(D-rank):
        destinations[k+1, rank+k] = 1
    su = np.column_stack([source_basis, null_space(source_basis.conj().T)])
    du = np.column_stack([destinations, null_space(destinations.conj().T)])
    U = du@su.conj().T
    vlocal = np.zeros((4, dc), complex)
    for l in range(rank):
        vlocal[:, l+1] = singular[l] * a[:, l]
    vlocal = vlocal.reshape(-1)
    h = val*np.outer(vlocal, vlocal.conj())
    Uf = np.kron(np.eye(4), U)
    hf = np.kron(h, np.eye(D))
    result += Eps.conj().T@Uf.conj().T@hf@Uf@Eps
    max_kernel_residual = max(max_kernel_residual, norm(hf@Uf@Eps@G))
    max_unitarity_error = max(max_unitarity_error, norm(U.conj().T@U-np.eye(dc*D)))
lines += [f'positive dilation reconstruction error={norm(result-A):.6e}',
          f'positive dilation maximum term/code residual={max_kernel_residual:.6e}',
          f'positive dilation maximum unitary residual={max_unitarity_error:.6e}']
assert norm(result-A) < 1e-10
assert max_kernel_residual < 1e-10

# The boundary Green function has b distinct visible poles at fixed gap and d=4.
mu, hop = 3., .5
for b in (2, 5, 12):
    one = mu*np.eye(b)-hop*(np.eye(b, k=1)+np.eye(b, k=-1))
    e, u = eigh(one)
    k = np.arange(1, b+1)
    exact_e = mu-2*hop*np.cos(k*np.pi/(b+1))
    exact_w = 2/(b+1)*np.sin(k*np.pi/(b+1))**2
    lines.append(f'boundary-pole example b={b}: energy error={norm(e-exact_e):.3e}, residue error={norm(abs(u[0,:])**2-exact_w):.3e}')

with open('checks.txt', 'w') as f:
    f.write('\n'.join(lines)+'\n')
print('\n'.join(lines))
