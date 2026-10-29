"""C118 full-spectrum utilities, adapted from C116 after reading its code."""
from collections import defaultdict
import numpy as np
from scipy.linalg import eigvalsh

def spectrum(model, length):
    d = 3 if model == 'SU3' else 2
    sectors = defaultdict(list)
    powers = d**np.arange(length, dtype=np.int64)
    for state in range(d**length):
        digits = (state // powers) % d
        counts = tuple(int(np.count_nonzero(digits == a)) for a in range(d-1))
        sectors[counts].append(state)
    energies = []
    maxdim = 0
    trace = 0.0
    trace2 = 0.0
    for basis in sectors.values():
        maxdim = max(maxdim, len(basis))
        index = {state: j for j, state in enumerate(basis)}
        h = np.zeros((len(basis), len(basis)))
        for col, state in enumerate(basis):
            digits = (state // powers) % d
            for i in range(length):
                j = (i+1) % length
                a, b = int(digits[i]), int(digits[j])
                if model == 'SU3':
                    swapped = state + (b-a)*int(powers[i]) + (a-b)*int(powers[j])
                    h[index[swapped], col] += 1.0
                else:
                    if model == 'Heisenberg':
                        h[col, col] += (0.25 if a == b else -0.25)
                    if a != b:
                        swapped = state + (b-a)*int(powers[i]) + (a-b)*int(powers[j])
                        h[index[swapped], col] += 0.5
        assert np.max(np.abs(h-h.T)) < 1e-13
        ev = eigvalsh(h, overwrite_a=False, check_finite=False)
        trace += float(np.trace(h))
        trace2 += float(np.sum(h*h))
        energies.extend(ev.tolist())
    energies = np.sort(energies)
    assert len(energies) == d**length
    assert abs(np.sum(energies)-trace) < 1e-8*max(1,abs(trace))
    assert abs(np.dot(energies, energies)-trace2) < 1e-8*max(1,trace2)
    return energies, maxdim

