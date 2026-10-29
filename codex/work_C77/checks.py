"""Small algebra checks for C77. No thermodynamic central-charge estimate."""
import itertools
import json
from pathlib import Path
import numpy as np


def projected_chain(q, L, N, t=1., V=.23):
    basis = [x for x in itertools.product(range(q + 1), repeat=L)
             if sum(a != 0 for a in x) == N]
    index = {x: i for i, x in enumerate(basis)}
    h = np.zeros((len(basis), len(basis)))
    word_errors = 0
    for i, x in enumerate(basis):
        h[i, i] = V * sum(x[j] != 0 and x[j+1] != 0 for j in range(L-1))
        for j in range(L-1):
            if (x[j] == 0) != (x[j+1] == 0):
                y = list(x)
                y[j], y[j+1] = y[j+1], y[j]
                y = tuple(y)
                h[index[y], i] -= t
                word_errors += tuple(a for a in x if a) != tuple(a for a in y if a)
    return h, word_errors

q, L, N = 2, 6, 3
h_color, word_errors = projected_chain(q, L, N)
h_charge, _ = projected_chain(1, L, N)
spectrum_error = np.max(np.abs(np.linalg.eigvalsh(h_color) -
                              np.repeat(np.linalg.eigvalsh(h_charge), q**N)))
assert word_errors == 0
assert spectrum_error < 1.e-12

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
pauli = [X, Y, Z]
gamma, u, D, b = .17, .4, .5, .23
h = ((1+gamma)*np.kron(X, X) + (1-gamma)*np.kron(Y, Y)
     + u*np.kron(Z, Z) + D*(np.kron(X, Y)-np.kron(Y, X))
     + b/2*(np.kron(Z, I)+np.kron(I, Z)))
comm = lambda a, c: a @ c - c @ a
hermitian_error = np.linalg.norm(h-h.conj().T)
parity_error = np.linalg.norm(comm(h, np.kron(Z, Z)))
charge_map = np.column_stack([comm(h, np.kron(p,I)+np.kron(I,p)).reshape(-1)
                             for p in pauli])
charge_singular = np.linalg.svd(charge_map, compute_uv=False)
assert hermitian_error < 1.e-12 and parity_error < 1.e-12
assert min(charge_singular) > 1.e-8

# Purely thermodynamic functions compatible with entropy and velocity counts.
# These are not constructed lattice models.
logd = np.log(2.)
thermo = []
for c in [1., 1.5, 10., 100.]:
    v = 1.
    A = np.pi*c/(3*v)
    tau = logd/A
    T = tau*np.logspace(-6, 6, 1000)
    s = logd*(T/tau)/np.sqrt(1+(T/tau)**2)
    assert np.all(s > 0) and np.all(s <= logd)
    assert np.all(np.diff(s) > 0)
    assert abs(s[0]/T[0]/A-1) < 1.e-10
    thermo.append({'formal_c':c, 'linear_entropy_slope':A, 'crossover':tau})

out = dict(hardcore_check=dict(q=q, L=L, N=N, colored_dimension=len(h_color),
           charge_dimension=len(h_charge), word_violations=int(word_errors),
           spectrum_max_error=float(spectrum_error)),
           candidate_algebra=dict(gamma=gamma,u=u,D=D,b=b,
           hermiticity_error=float(hermitian_error),parity_error=float(parity_error),
           additive_uniform_charge_singular_values=charge_singular.tolist()),
           formal_thermodynamic_checks=thermo)
Path('checks.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
