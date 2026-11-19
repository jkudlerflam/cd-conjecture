"""Exact checks of the residue-word expansion and the gate example in PROOF.md.

Run with python3 -B check_resolvent.py. All arithmetic is rational.
"""

from fractions import Fraction as F
from itertools import product


FLAGS = 3
SIZE = 1 << FLAGS
ENERGY = [state.bit_count() for state in range(SIZE)]
PORTS = (0b011, 0b110)


def apply_ports(vector, ports):
    return [sum((vector[s ^ mask] for mask in ports), F(0)) for s in range(SIZE)]


def perturbation_energies(ports, max_order):
    # Ground component of the eigenvector is fixed to one at every coupling.
    waves = [[F(i == 0) for i in range(SIZE)]]
    energies = [F(0)]
    for order in range(1, max_order + 1):
        source = apply_ports(waves[-1], ports)
        energies.append(source[0])
        wave = [F(0)]
        for i in range(1, SIZE):
            feedback = sum(energies[j] * waves[order-j][i] for j in range(1, order+1))
            wave.append((feedback-source[i])/ENERGY[i])
        waves.append(wave)
    return energies


def contour_energies(ports, max_order):
    # Trace of zeta times the resolvent around the unique low level.
    result = [F(0)]
    for order in range(1, max_order + 1):
        total = F(0)
        for choices in product(ports, repeat=order):
            for start in range(SIZE):
                state = start
                path = [ENERGY[start]]
                for mask in choices:
                    state ^= mask
                    path.append(ENERGY[state])
                degree = path.count(0) - 2
                if state != start or degree < 0:
                    continue
                coefficients = [F(1)] + [F(0)] * degree
                for energy in filter(None, path):
                    series = [-F(1, energy**(j+1)) for j in range(degree+1)]
                    coefficients = [
                        sum(coefficients[k] * series[j-k] for k in range(j+1))
                        for j in range(degree+1)
                    ]
                total += coefficients[degree]
        result.append(total)
    return result


def eliminate(matrix, rhs=None):
    # All matrices in this check are positive definite.
    n = len(matrix)
    pivots = []
    for j in range(n):
        pivot = matrix[j][j]
        assert pivot > 0
        pivots.append(pivot)
        for i in range(j+1, n):
            if not matrix[i][j]:
                continue
            factor = matrix[i][j]/pivot
            for k in range(j, n):
                matrix[i][k] -= factor*matrix[j][k]
            if rhs is not None:
                rhs[i] -= factor*rhs[j]
    if rhs is not None:
        for i in reversed(range(n)):
            rhs[i] = (rhs[i] - sum(matrix[i][j]*rhs[j] for j in range(i+1, n)))/pivots[i]
    return pivots


def boundary_green(left, right):
    # Pair mass six. A single path has unit hopping and mass three.
    dimension = left*right
    matrix = [[F(0) for _ in range(dimension)] for _ in range(dimension)]
    for i in range(left):
        for j in range(right):
            p = i*right + j
            matrix[p][p] = F(6)
            if i:
                matrix[p][p-right] = matrix[p-right][p] = F(-1)
            if j:
                matrix[p][p-1] = matrix[p-1][p] = F(-1)
    source = [F(i == 0) for i in range(dimension)]
    eliminate(matrix, source)
    return source[0]


def main():
    through = 6
    left = contour_energies(PORTS[:1], through)
    right = contour_energies(PORTS[1:], through)
    both = contour_energies(PORTS, through)
    for ports, energies in [(PORTS[:1], left), (PORTS[1:], right), (PORTS, both)]:
        assert energies == perturbation_energies(ports, through)
    linked = [ab-a-b for ab, a, b in zip(both, left, right)]
    assert linked[1:] == [0, 0, 0, -F(1, 4), 0, F(3, 8)]
    print("connected energy coefficients at orders 1 to 6:", ", ".join(map(str, linked[1:])))
    lengths = range(1, 6)
    gram = [[boundary_green(a, b) for b in lengths] for a in lengths]
    assert gram == list(map(list, zip(*gram)))
    pivots = eliminate(gram)
    print("exact positive pivots of the gate Gram matrix:", len(pivots))


if __name__ == "__main__":
    main()
