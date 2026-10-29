"""Small numerical checks of C99's exact identities, not a CFT simulation."""
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.])

def entropy(p):
    p = np.real(np.asarray(p))
    p = p[p > 1e-14]
    return -np.dot(p, np.log(p))

def reduced(psi, keep, n):
    rest = [j for j in range(n) if j not in keep]
    mat = psi.reshape([2]*n).transpose(list(keep)+rest).reshape(2**len(keep), -1)
    return mat @ mat.conj().T

def S(psi, keep, n):
    return entropy(np.linalg.eigvalsh(reduced(psi, keep, n)))

def tensor(ops):
    result = np.array([1.])
    for op in ops:
        result = np.kron(result, op)
    return result

# Full four-spin active chain, with internal evolution on both sides of the cut.
J = 1.0
p = 0.25
pair = np.array([0, np.sqrt(p), 1j*np.sqrt(1-p), 0], dtype=complex)
psi = np.kron(np.kron([1., 0.], pair), [1., 0.])
H = np.zeros((16, 16), complex)
for j in range(3):
    for pauli in (X, Y):
        ops = [I]*4
        ops[j] = ops[j+1] = pauli
        H += J/2*tensor(ops)
dt = 1e-5
rate = (S(expm(-1j*dt*H)@psi, [0, 1], 4)-S(expm(1j*dt*H)@psi, [0, 1], 4))/(2*dt)
expected = np.sqrt(3)/2*J*np.log(3)
active_error = abs(rate-expected)
assert active_error < 1e-7, active_error
H_periodic = H + J/2*(tensor([X, I, I, X])+tensor([Y, I, I, Y]))
casimir_energy_error = abs(np.linalg.eigvalsh(H_periodic)[0]+2*J/np.sin(np.pi/4))
assert casimir_energy_error < 1e-12

# Pure coherent Ising spectator, open boundaries. Correlations equal tanh(k)^r.
n = 8
bits = ((np.arange(2**n)[:, None] >> np.arange(n-1, -1, -1)) & 1)
spins = 1-2*bits
correlation_error = 0.
spectator_entropy_error = 0.
for k in (0.2, 0.8, 1.7, 3.):
    chi = np.exp(k/2*np.sum(spins[:, :-1]*spins[:, 1:], axis=1))
    chi = chi/np.linalg.norm(chi)
    probs = abs(chi)**2
    for r in range(1, n):
        actual = np.dot(probs, spins[:, 0]*spins[:, r])
        correlation_error = max(correlation_error, abs(actual-np.tanh(k)**r))
    energy = np.sum(bits, axis=1)  # Delta=1 onsite spectator Hamiltonian.
    initial_entropy = S(chi, list(range(4)), n)
    for t in (0.13, 1.2, 3.):
        value = S(chi*np.exp(-1j*t*energy), list(range(4)), n)
        spectator_entropy_error = max(spectator_entropy_error, abs(value-initial_entropy))
assert correlation_error < 1e-12
assert spectator_entropy_error < 1e-12
xi = -1/np.log(np.tanh(np.array([1., 2., 4., 6.])))
assert np.all(np.diff(xi*rate) > 0)

# Euclidean local filtering, computed directly from normalized Schmidt weights.
M = 0.7
filter_error = 0.
filter_rows = []
for rank in (2, 8, 64, 512):
    def filter_entropy(tau):
        weights = np.r_[0.5*np.exp(-2*M*tau),
                        np.full(rank, 0.5/rank*np.exp(2*M*tau))]
        weights /= weights.sum()
        return entropy(weights)
    derivative = (filter_entropy(dt)-filter_entropy(-dt))/(2*dt)
    expected_filter = M*np.log(rank)
    filter_error = max(filter_error, abs(derivative-expected_filter))
    filter_rows.append((rank, derivative, expected_filter))
assert filter_error < 1e-8

# Distant MI under an actual uniform nearest-neighbor XX-Pauli chain.
n = 7
A, E, B = 0, 1, 5
bits = ((np.arange(2**n)[:, None] >> np.arange(n-1, -1, -1)) & 1)
state = np.zeros(2**n, complex)
state[(bits[:, A] == 0) & (bits[:, E] == bits[:, B])] = 1.
state /= np.linalg.norm(state)
H = np.zeros((2**n, 2**n), complex)
for j in range(n-1):
    ops = [I]*n
    ops[j] = ops[j+1] = X
    H += tensor(ops)
evals, evecs = np.linalg.eigh(H)
coef = evecs.conj().T @ state

def mi(t):
    evolved = evecs @ (np.exp(-1j*t*evals)*coef)
    return S(evolved, [A], n)+S(evolved, [B], n)-S(evolved, [A, B], n)
mi_error = max(abs(mi(t)-entropy([np.sin(t)**2, np.cos(t)**2])) for t in (0., 0.1, np.pi/8, 0.7))
t = np.pi/8
mi_rate = (mi(t+dt)-mi(t-dt))/(2*dt)
mi_expected = 2*np.sin(t)*np.cos(t)*np.log(np.cos(t)**2/np.sin(t)**2)
assert mi_error < 1e-12
assert abs(mi_rate-mi_expected) < 1e-8

# Multisector local-join integral, numerical quadrature over an infinite interval.
epsilon = 7.
sectors = [(0.5, 0.03), (2., 0.5), (1.2, 4.)]
Acoef = sum(c*v for c, v in sectors)
integral = sum(quad(lambda t: c*v*v/(3*(epsilon**2+v*v*t*t)), 0, np.inf,
                    epsabs=1e-12, epsrel=1e-12)[0] for c, v in sectors)
local_join_error = abs(6*epsilon/np.pi*integral-Acoef)
assert local_join_error < 1e-10
lo, hi = 0.2, 12.
window_integral = sum(quad(lambda t: c*v*v/(3*(epsilon**2+v*v*t*t)), lo, hi)[0]
                      for c, v in sectors)
window_formula = sum(c*v/(3*epsilon)*(np.arctan(v*hi/epsilon)-np.arctan(v*lo/epsilon))
                     for c, v in sectors)
assert abs(window_integral-window_formula) < 1e-12

lines = [
    '# C99 checks', '',
    '[Numerical] All assertions passed in python3 checks.py. These are small matrix and quadrature checks of the analytic proofs, not evidence for or against the general Casimir bound.', '',
    f'Four-spin active-chain entropy derivative: {rate:.12g}. Analytic value: {expected:.12g}. Absolute error: {active_error:.3g}.',
    f'Periodic four-spin XX ground energy agrees with -2J/sin(pi/L) to {casimir_energy_error:.3g}. The asymptotic Casimir expansion is analytic, not fitted.',
    f'Coherent Ising spectator: maximum correlation error {correlation_error:.3g}; maximum entropy change under onsite real-time evolution {spectator_entropy_error:.3g}.',
    'Spectator correlation lengths at k=1,2,4,6: '+', '.join(f'{v:.9g}' for v in xi)+'.',
    f'Euclidean one-qubit filtering: maximum derivative error {filter_error:.3g}.',
    f'Seven-spin uniform chain: maximum distant-MI identity error {mi_error:.3g}; rate at t=pi/8 is {mi_rate:.12g}, predicted {mi_expected:.12g}.',
    f'Multisector local-join positive integral: recovered A={6*epsilon/np.pi*integral:.12g}, expected {Acoef:.12g}, error {local_join_error:.3g}. Finite-window arctangent identity also passed.', '',
    '[Proved analytically, not established by finite sampling] Divergence of xi_Z r_0 and of the Euclidean rate M log N follows from their displayed formulas. Finite calculations check signs, factors, tensor assignments, and the cancellation of within-half terms.', '',
    '[Not tested] No uniform continuum matching, general capacity lower bound, or microscopic violation of A<=C M log d has been established.',
]
Path(__file__).with_name('CHECKS.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
