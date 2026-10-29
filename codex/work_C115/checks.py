"""Small deterministic checks for C115. Run from this directory only."""
from pathlib import Path
import math
import numpy as np

lines = ["# C115 checks", "",
         "[Numerical] Deterministic algebra and scalar sums only. Largest matrix is 16 by 16. No many-body infrared extrapolation.", ""]


def record(message):
    lines.extend(["[Numerical] PASS. " + message, ""])


ident = np.eye(2, dtype=complex)
sx = np.array([[0., 1.], [1., 0.]], dtype=complex)
sy = np.array([[0., -1j], [1j, 0.]], dtype=complex)
sz = np.diag([1., -1.]).astype(complex)


def embed_two(a, i, b, j):
    value = np.ones((1, 1), dtype=complex)
    for k in range(4):
        value = np.kron(value, a if k == i else b if k == j else ident)
    return value


swap = np.zeros((16, 16))
for left in range(4):
    for right in range(4):
        swap[4*right+left, 4*left+right] = 1.

worst = 0.
for epsilon in [0.1, 0.01, 0.001]:
    tf, ts = 1/(1+epsilon), epsilon/(1+epsilon)
    bond = tf*(embed_two(sx, 0, sx, 2)+embed_two(sy, 0, sy, 2))/2
    bond += ts*(embed_two(sx, 1, sx, 3)+embed_two(sy, 1, sy, 3))/2
    values, vectors = np.linalg.eigh(bond)
    error = abs((values[-1]-values[0])/2-1)
    worst = max(worst, error)
    assert error < 1e-12
    assert np.linalg.norm(swap @ bond @ swap-bond) < 1e-12
    for n in [0, -1]:
        mat = vectors[:, n].reshape(4, 4)
        rho_left = mat @ mat.conj().T
        rho_right = mat.T @ mat.conj()
        assert np.linalg.norm(rho_left-rho_right) < 1e-12
    for onsite in [np.kron(sx, sy), .3*np.kron(sz, ident)+.7*np.kron(sy, sx)]:
        coboundary = np.kron(onsite, np.eye(4))-np.kron(np.eye(4), onsite)
        assert np.linalg.norm(swap @ coboundary @ swap+coboundary) < 1e-12
        gauged = np.linalg.eigvalsh(bond+coboundary)
        assert (gauged[-1]-gauged[0])/2 >= 1-1e-12
    vf, vs = 2*tf, 2*ts
    assert vs < 1 < vf and abs(vf+vs-2) < 1e-12
record(f'Fast-plus-slow XX bond has half-spread one, equal extremal marginals and the swap gauge witness. Maximum spectral error {worst:.3e}. N(1)=1 and A=2 throughout.')

u, x, v = 1., 4., 2.
assert v*x > u*x and u*x/v == 2.
record('At v=2u the energy window includes dimensions only through x/2, not through x.')

for onsite in [sx, sy, sz, .3*sx+.2*sy-.7*sz]:
    coboundary = np.kron(onsite, ident)-np.kron(ident, onsite)
    assert abs(coboundary[0, 3]) < 1e-14
record('A two-site transition changing both site letters has zero matrix element against every tested coboundary, consistent with the analytic M>=t witness for the spatial-cell vacuum hop.')

tau = .25
trace_values = []
for epsilon in [1e-1, 1e-2, 1e-3, 1e-4]:
    vs = 2*epsilon/(1+epsilon)
    bound = 1/(-np.expm1(-2*np.pi*tau*vs))
    count = math.ceil(1/vs)  # m>=0 with m<1/vs, x=u=1
    trace_values.append(bound)
    lines.append(f'[Numerical] epsilon={epsilon:g}, v_s={vs:.9g}, oscillator trace lower bound={bound:.9g}, strict-window state lower count={count}.')
assert np.all(np.diff(trace_values) > 0)
lines.append('')
record('The lower trace and window count grow while d=4, M/u=1 and N(u)=1 stay fixed. These are evaluations of proved lower bounds, not finite-L trace fits.')


def log_boson_q(aspect, radius):
    levels = np.arange(1., 241.)
    charges = np.arange(-160., 161.)
    oscillator = -2*np.sum(np.log1p(-np.exp(-2*np.pi*aspect*levels)))
    theta_n = np.sum(np.exp(-np.pi*aspect*charges**2/radius**2))
    theta_w = np.sum(np.exp(-np.pi*aspect*charges**2*radius**2))
    return oscillator+np.log(theta_n)+np.log(theta_w)


worst_modular = 0.
for radius in [0.5, 1., np.sqrt(2.), 3., 10.]:
    for aspect in [.25, .5, .75, 1., 2., 4.]:
        logz = np.pi*aspect/6+log_boson_q(aspect, radius)
        logz_dual = np.pi/(6*aspect)+log_boson_q(1/aspect, radius)
        worst_modular = max(worst_modular, abs(logz-logz_dual))
        assert abs(logz-logz_dual) < 1e-11
    assert log_boson_q(.5, radius) >= np.pi/4-1e-12
bxx = log_boson_q(.5, np.sqrt(2.))
record(f'Compact-boson modular identity and log Q(1/2)>=pi/4 at five radii. Maximum residual {worst_modular:.3e}; b_XX={bxx:.12g}.')

gamma = .5
indices = np.arange(200.)
geometric = np.sum(2**(-(1-gamma)*indices))
assert abs(geometric-(2+np.sqrt(2))) < 1e-13
linear = np.sum((indices+1)*2**(-indices))
assert abs(linear-4) < 1e-13
assert abs(4*np.e*geometric-4*np.e*(2+np.sqrt(2))) < 1e-12
record(f'Dyadic constants: C_W=4e(2+sqrt(2))={4*np.e*geometric:.12g}; C_S=16e(2+sqrt(2))/pi={16*np.e*geometric/np.pi:.12g}; linear-window sum=4.')

stacks = [[1.], [.5, .5], [.7, .2, .1], [1/8]*8, [2.**(-k) for k in range(1, 12)]]
for hoppings in stacks:
    speeds = 2*np.array(hoppings)
    mass = sum(hoppings)
    edges = np.unique(np.r_[0., speeds])
    integral = sum((b-a)*np.count_nonzero(speeds >= (a+b)/2) for a, b in zip(edges[:-1], edges[1:]))
    assert abs(integral-2*mass) < 1e-12
    for cutoff in np.r_[speeds, .5*speeds, 2*speeds]:
        number = np.count_nonzero(speeds >= cutoff)
        assert number <= len(hoppings)
        assert cutoff*number <= 2*mass+1e-12
record('Exact step-function integration on five XX stacks gives integral N du=A=2M and N(u)<=min(m,2M/u), including velocity endpoints.')

# This is a stipulated envelope model, not a fit to campaign spectra.
for n in range(1, 101):
    kappa, scale = .8, 2.
    payload_speed = scale*np.exp(-kappa*n)
    for cutoff in np.geomspace(1e-36, 4., 120):
        number = 1+(n if payload_speed >= cutoff else 0)
        envelope = 1+max(0., np.log(scale/cutoff)/kappa)
        assert number <= envelope+1e-12
record('A stipulated n-payload exponential envelope plus one length mode obeys the logarithmic N(u) bound for n<=100. This checks the conditional algebra only.')

documented_upper_logs = {1: -13.8, 2: -20.8, 3: -29.5, 4: -38.1, 23: -205.2}
for n, logcell in documented_upper_logs.items():
    logphysical = logcell+np.log10(n+1)
    assert logphysical < -13
    lines.append(f'[Numerical, inherited inputs] n={n}, log10 upper physical velocity coefficient={logphysical:.6f}, after multiplying the recorded cell velocity by b=n+1. With the proved M>=t witness, v_max/M is at most this coefficient times lambda^2 g^2 alpha_n^4 in the recorded effective description.')
lines.extend(['', '[Numerical] All assertions passed. No universal interacting counting estimate is tested or proved by this script.', ''])
Path('CHECKS.md').write_text('\n'.join(lines))
print('\n'.join(lines))
