"""Checks the noncommuting-matrix identity in IR_SOURCE_REVIEW.md."""
import json
from pathlib import Path
import numpy as np

rng = np.random.default_rng(20261028)
results = []
for n in (2, 4, 6):
    velocity = np.linspace(.7, 1.4, n) * np.where(np.arange(n) % 2, -1, 1)
    z = np.diag(1+.01*rng.random(n))
    coupling = rng.normal(size=(n, n))
    coupling = .01*(coupling+coupling.T)/2
    lambda_z = np.linalg.solve(z, coupling @ z)
    u = np.eye(n)+.01*rng.normal(size=(n, n))
    xi = .35
    def bubble(x):
        return np.diag((-1j+velocity*x)/(1j+velocity*x)/(4*np.pi*np.abs(velocity)))
    a, b = bubble(xi), bubble(0)
    identity = np.eye(n)
    zminus2 = np.linalg.inv(z @ z)
    ca = np.linalg.solve(identity+a@lambda_z, a) @ zminus2
    cb = np.linalg.solve(identity+b@lambda_z, b) @ zminus2
    factor = np.linalg.solve(identity+a@lambda_z, a-b) @ np.linalg.inv(identity+lambda_z@b) @ zminus2
    error = float(np.max(np.abs(ca-cb-factor)))
    assert error < 1e-12
    predicted_middle = velocity*xi/(2*np.pi*np.abs(velocity)*(1j+velocity*xi))
    middle_error = float(np.max(np.abs(np.diag(a-b)-predicted_middle)))
    assert middle_error < 1e-12
    physical = u @ factor @ u.T
    minimum = float(np.linalg.svd(physical, compute_uv=False).min())
    assert minimum > 0
    results.append(dict(branches=n, factorization_error=error,
                        middle_formula_error=middle_error,
                        physical_smallest_singular_value=minimum))
report = dict(matrix_checks=results,
              bounded_ratio_counterexample_at_scale_minus_1000=float(2**(-.05*1000)),
              passed=True)
Path(__file__).with_name('ir_source_results.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
