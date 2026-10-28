"""Independent finite Fourier checks for ../sources/BRANCH_SOURCES.md."""
from pathlib import Path
import json
import numpy as np

length = 11
momenta = 2*np.pi*np.arange(length)/length
branches = momenta[[1, 4, 8]]
radius = len(branches)-1
fourier = np.exp(1j*np.outer(np.arange(length), momenta))/np.sqrt(length)
records = []
for a, ka in enumerate(branches):
    coefficients = {0: 1+0j}
    denominator = 1.0
    for b, kb in enumerate(branches):
        if b == a:
            continue
        factor = {0: .5, 1: -.25*np.exp(-1j*kb), -1: -.25*np.exp(1j*kb)}
        new = {}
        for r, coefficient in coefficients.items():
            for s, value in factor.items():
                new[r+s] = new.get(r+s, 0)+coefficient*value
        coefficients = new
        denominator *= np.sin((ka-kb)/2)**2
    coefficients = {r: c/denominator for r, c in coefficients.items()}
    def p(k):
        return sum(c*np.exp(1j*r*k) for r, c in coefficients.items())
    interpolation_error = max(abs(p(kb)-(1 if a == b else 0)) for b, kb in enumerate(branches))
    l1 = sum(abs(c) for c in coefficients.values())
    densities = []
    for x in range(length):
        matrix = np.zeros((length, length), dtype=complex)
        for r, c in coefficients.items():
            matrix[x, (x+r) % length] += c/2
            matrix[(x-r) % length, x] += c/2
        assert np.max(np.abs(matrix-matrix.conj().T)) < 1e-12
        densities.append(matrix)
    max_vertex_error = 0.0
    for qindex, q in enumerate(momenta):
        oq = sum(np.exp(-1j*q*x)*matrix for x, matrix in enumerate(densities))
        vertex = fourier.conj().T @ oq @ fourier
        expected = np.zeros_like(vertex)
        for kindex, k in enumerate(momenta):
            expected[kindex, (kindex+qindex) % length] = (p(k)+p(k+q))/2
        max_vertex_error = max(max_vertex_error, float(np.max(np.abs(vertex-expected))))
    assert interpolation_error < 1e-12
    assert l1 <= 1/denominator+1e-12
    assert max_vertex_error < 1e-11
    records.append(dict(branch=a, interpolation_error=float(interpolation_error),
                        coefficient_l1=float(l1), l1_bound=float(1/denominator),
                        maximum_vertex_error=max_vertex_error))
result = dict(length=length, branches=len(branches), transferred_momenta=length,
              records=records, passed=True)
Path(__file__).with_name('source_results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
