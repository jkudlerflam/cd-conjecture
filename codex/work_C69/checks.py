"""Elementary checks for C69. Run only from this work folder."""
from pathlib import Path
import math


def log_cosh(x):
    y = abs(x)
    return y + math.log1p(math.exp(-2.0 * y)) - math.log(2.0)


sectors = [(0.5, 0.2), (1.0, 1.7), (3.0, 4.0)]
ell = 7.0
c = sum(ca for ca, _ in sectors)
A = sum(ca * va for ca, va in sectors)
B = sum(ca / va for ca, va in sectors)


def entropy_increment(t):
    return sum(ca / 6.0 * log_cosh(math.pi * va * t / (2.0 * ell))
               for ca, va in sectors)


def entropy_rate(t):
    return sum(math.pi * ca * va / (12.0 * ell)
               * math.tanh(math.pi * va * t / (2.0 * ell))
               for ca, va in sectors)


errors = []
for t in [0.0, 0.1, 1.0, 10.0, 100.0]:
    step = 1e-4 * max(1.0, t)
    derivative = (entropy_increment(t + step) - entropy_increment(t - step)) / (2.0 * step)
    errors.append(abs(derivative - entropy_rate(t)))
assert max(errors) < 1e-8
expected_rate = math.pi * A / (12.0 * ell)
assert abs(entropy_rate(10000.0) - expected_rate) < 1e-14
assert entropy_rate(0.0) == 0.0
q = 4.0
observation_time = q * ell / min(v for _, v in sectors)
finite_time_fraction = 1.0 - 2.0 * math.log(2.0) / (math.pi * q)
assert entropy_increment(observation_time) / observation_time >= finite_time_fraction * expected_rate

# The log anomaly has canceled from entropy increments.
scale = 11.0
for t in [0.1, 1.0, 100.0]:
    scaled = sum(ca / 6.0 * log_cosh(math.pi * va * (scale * t) / (2.0 * scale * ell))
                 for ca, va in sectors)
    assert abs(scaled - entropy_increment(t)) < 1e-13
    scaled_rate = sum(math.pi * ca * va / (12.0 * scale * ell)
                      * math.tanh(math.pi * va * (scale * t) / (2.0 * scale * ell))
                      for ca, va in sectors)
    assert abs(scale * scaled_rate - entropy_rate(t)) < 1e-13

v_H = c / B
assert min(v for _, v in sectors) <= v_H <= max(v for _, v in sectors)
assert abs(A * B - c * c) > 1.0
assert A * B >= c * c

# Check the exponential-series majorant used in the LR path count.
max_lr_ratio = 0.0
for t in [0.0, 0.01, 0.2, 1.0, 3.0]:
    x = 6.0 * t
    terms = [1.0]
    for n in range(1, 300):
        terms.append(terms[-1] * x / n)
    for r in [1, 2, 5, 20]:
        tail = math.fsum(terms[r:])
        bound = math.exp(-r + x * math.e)
        assert tail <= bound * (1.0 + 1e-14)
        max_lr_ratio = max(max_lr_ratio, tail / bound)

for n in [1, 2, 10, 1000]:
    assert math.log1p(1.0 / n) >= 1.0 / (n + 1.0)
    assert 3.0 / math.log1p(1.0 / n) <= 3.0 * (n + 1.0) * (1.0 + 1e-14)

# Literal pi/6 convention with u=v/2 must reproduce pi/3.
u_B = sum(ca / (va / 2.0) for ca, va in sectors)
assert abs(math.pi * B / 3.0 - math.pi * u_B / 6.0) < 1e-14

lines = [
    '# C69 elementary checks',
    '',
    '[Numerical] All assertions in checks.py passed. These are formula checks, not a microscopic CFT simulation.',
    '',
    f'Maximum centered-difference error for the strip entropy derivative: {max(errors):.6g}.',
    f'Multisector late-time entropy rate: {expected_rate:.12g}.',
    f'Finite-time lower fraction at q=4: {finite_time_fraction:.12g}; the integrated strip rate satisfies the bound.',
    f'Weighted harmonic speed: {v_H:.12g}, between 0.2 and 4.',
    f'Maximum sampled exponential-tail / LR-majorant ratio: {max_lr_ratio:.6g}.',
    '',
    '| Constant | Value |',
    '|:--|--:|',
    f'| SIE conditional coefficient 216/pi | {216.0 / math.pi:.12g} |',
    f'| Thermal coefficient 3/pi | {3.0 / math.pi:.12g} |',
    f'| Thermal coefficient at alpha=0.9 | {3.0 / (0.9 * math.pi):.12g} |',
    f'| Literal pi/6 coefficient 6/pi | {6.0 / math.pi:.12g} |',
    f'| Spatial doubling coefficient 3/log 2 | {3.0 / math.log(2.0):.12g} |',
    f'| Heisenberg-required coefficient pi/log 2 | {math.pi / math.log(2.0):.12g} |',
    f'| Conservative LR coefficient 6e | {6.0 * math.e:.12g} |',
    '',
    '[Proved in RESULT.md, tested here] Equal spatial scaling leaves entropy increments invariant and divides rates by the scale. The pi/6 convention with u=v/2 agrees with the ordinary nonchiral pi/3 formula. The one-site spatial bound is at most 3(n+1) log d/alpha.',
    '',
    '[Assumed/open] These checks do not verify uniform lattice preparation, dynamical CFT convergence, or the conjectured universal bound on A.',
]
Path(__file__).with_name('CHECKS.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines[:9]))
