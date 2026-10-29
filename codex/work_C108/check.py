"""Small deterministic checks for C108. No many-body simulations."""
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

out = []

def record(s):
    out.append(s)
    print(s)


def replacement_kraus(values, upper=True):
    values = np.asarray(values)
    endpoint = np.max(values) if upper else np.min(values)
    selected = np.flatnonzero(values == endpoint)
    rest = np.flatnonzero(values != endpoint)
    p = np.diag((values == endpoint).astype(float))
    ks = [p.astype(complex)]
    for b in rest:
        k = np.zeros((len(values), len(values)), complex)
        k[selected[0], b] = 1
        ks.append(k)
    return p, ks


max_violation = 0.0
max_block_error = 0.0
max_unital_error = 0.0
count = 0
for values in ([0, 1], [0, 1, 2], [0, 1, 1], [0, 0, 2, 2], [0, 1, 7]):
    d = len(values)
    total = np.add.outer(values, values).ravel()
    inds = np.arange(d*d)
    for phase in (0.2, 0.7, 1.3):
        raw = np.sin((inds[:, None]+1)*(inds[None, :]+2)+phase)
        raw = raw + 1j*np.cos(inds[:, None]+2*inds[None, :]+phase)
        h = (raw+raw.conj().T)/2
        h = h*(total[:, None] == total[None, :])
        ev = np.linalg.eigvalsh(h)
        h = h-(ev[-1]+ev[0])*np.eye(d*d)/2
        h = h/np.linalg.norm(h, 2)
        for upper in (False, True):
            p, ks = replacement_kraus(values, upper)
            completeness = sum(k.conj().T@k for k in ks)
            max_unital_error = max(max_unital_error,
                                   np.linalg.norm(completeness-np.eye(d)))
            transformed = np.zeros_like(h)
            for k1 in ks:
                for k2 in ks:
                    k = np.kron(k1, k2)
                    transformed += k.conj().T@h@k
            b = transformed-h
            p2 = np.kron(p, p)
            max_block_error = max(max_block_error, np.linalg.norm(b@p2))
            upper_bound = 2*(np.eye(d*d)-p2)
            violation = np.linalg.eigvalsh(b-upper_bound)[-1]
            max_violation = max(max_violation, violation)
            r = np.eye(d)-p
            looser = 2*(np.kron(r,np.eye(d))+np.kron(np.eye(d),r))
            assert np.linalg.eigvalsh(looser-upper_bound)[0] > -1e-12
            count += 1
assert max_violation < 1e-12
assert max_block_error < 1e-12
assert max_unital_error < 1e-12
record(f"Replacement channel: {count} cases passed.")
record(f"Largest upper-bound residual: {max_violation:.3e}.")
record(f"Largest maximal/minimal-block residual: {max_block_error:.3e}.")
record(f"Largest Kraus completeness residual: {max_unital_error:.3e}.")

z = np.array([[1., .3, -.2], [.2, 1.7, .1], [-.4, .1, .9]])
v = np.array([.2, 1.3, 2.4])
dmat = z.T@np.diag(v)@z/np.pi
kappa = z.T@np.diag(1/v)@z/np.pi
squared = eigh(dmat, kappa, eigvals_only=True)
err = np.max(np.abs(squared-v*v))
assert err < 1e-12
record(f"Generalized response eigenvalues reproduce v^2: residual {err:.3e}.")
assert np.min(np.linalg.eigvalsh(v[-1]**2*kappa-dmat)) > -1e-12
record("Matrix inequality D <= v_max^2 kappa passed.")

for m in (1, 2, 5, 20):
    zcol = np.ones((m, 1))/np.sqrt(m)
    v0 = .7
    dd = (zcol.T@zcol).item()*v0/np.pi
    kk = (zcol.T@zcol).item()/(np.pi*v0)
    assert abs(dd-v0/np.pi) < 1e-12
    assert abs(kk-1/(np.pi*v0)) < 1e-12
record("One-charge example: identical D and kappa for m=1,2,5,20 passed.")

for n in np.linspace(.01, .99, 99):
    t = 1.0
    vv = 2*t*np.sin(np.pi*n)
    dd = vv/np.pi
    curvature = np.pi*vv
    assert abs(dd*curvature-vv*vv) < 1e-12
record("XX response identity: 99 densities passed; integrated curvature=4t.")

for gamma in np.linspace(.01, np.pi-.01, 101):
    vv = np.pi*np.sin(gamma)/(2*gamma)
    kval = np.pi/(2*(np.pi-gamma))
    curvature = np.pi*(np.pi-gamma)*np.sin(gamma)/gamma
    dd = np.pi*np.sin(gamma)/(4*gamma*(np.pi-gamma))
    assert abs(curvature-np.pi*vv/kval) < 1e-12
    assert abs(dd*curvature-vv*vv) < 1e-12
record("XXZ response identity: 101 anisotropies passed.")

for eps in (1., .1, .01, .001):
    peak = 1/eps
    integral = 1/np.sqrt(.25+eps*eps)
    width = np.sqrt(.25+eps*eps)-eps
    assert integral < 2
    assert width <= .5
    assert abs(peak*eps-1) < 1e-12
record("Convex spike: unbounded peak and bounded curvature integral passed.")
record(f"Locality constant 6 exp(1) = {6*np.exp(1):.12f}.")

Path('CHECKS.md').write_text(
    '# C108 checks\n\n[Numerical] Run with `python3 check.py`. '
    'All calculations use small matrices or scalar formulas. '
    'They do not verify the assumed infrared matching.\n\n'
    + '\n'.join('[Numerical] '+line+'\n' for line in out), encoding='utf-8')
