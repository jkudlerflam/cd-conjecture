"""Numerical checks of the new elementary bounds and normalization."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from check_c129 import energy, bond_M
from bethe_scan import solve, DELTAS

sx = np.array([[0, 1], [1, 0]])/2
sy = np.array([[0, -1j], [1j, 0]])/2
sz = np.diag([.5, -.5])
eye = np.eye(2)
paulis = [sx, sy, sz]
max_M_error = 0.
for jp, dm, jz, B in [(1, 0, 0, 0), (1, .3, 0, 0), (1, 0, 1.5, 0),
                       (1, .3, -.99, .8), (-.7, .9, -2., 3.), (0, 0, 2, .4)]:
    h = jp*(np.kron(sx,sx)+np.kron(sy,sy))+dm*(np.kron(sx,sy)-np.kron(sy,sx))
    h += jz*np.kron(sz,sz)-B/2*(np.kron(sz,eye)+np.kron(eye,sz))
    def spread(p):
        o = sum(p[k]*paulis[k] for k in range(3))
        return np.ptp(np.linalg.eigvalsh(h+np.kron(o,eye)-np.kron(eye,o)))/2
    for initial in [[0,0,0], [.4,-.8,.9]]:
        opt = minimize(spread, initial, method='Nelder-Mead',
                       options=dict(maxiter=4000, xatol=1e-10, fatol=1e-12))
        max_M_error = max(max_M_error, abs(opt.fun-bond_M(np.hypot(jp,dm),jz,B)))
assert max_M_error < 2e-8

count, bound_violation, endpoint_error, iterative_error = 0, 0., 0., 0.
for L in range(2, 7):
    for n in range(1, L//2+1):
        for delta in [-1, -.999999, -.999, -.9, 0, 1, 3, 10]:
            small, large = energy(L, n, delta), energy(2*L, 2*n, delta)
            M = bond_M(1, delta)
            actual = large-2*small
            k, m = 2*n, 2*L
            Ak = np.sin(k*np.pi/m)/np.sin(np.pi/m)
            free_bound = (k-Ak)+max(delta,0)*(k*k-Ak*Ak)/m
            cross_bound = (delta+1)*k*(k-1)/(m-1)
            coarse_dilute = 2*np.pi**2*n**3/L*M/L
            bound_violation = max(bound_violation, actual-free_bound, actual-cross_bound,
                                  actual-coarse_dilute)
            assert actual <= min(free_bound, cross_bound, coarse_dilute)+2e-9
            if delta == -1:
                endpoint_error = max(endpoint_error, abs(small+L/4), abs(large+L/2))
            count += 1
for delta in DELTAS:
    a, _ = solve(8, 4, delta)
    b, _ = solve(8, 4, delta, iterative_threshold=0)
    iterative_error = max(iterative_error, float(abs(a-b)))
assert iterative_error < 2e-10

# Compare the explicitly requested 36 supplied points, to their printed accuracy.
reference = Path('../../contrib/rigorous_toolbox_map/d2_sanity_output.txt').read_text().splitlines()[:36]
ours = json.loads(Path('sanity_check.json').read_text())
reference_error = 0.
for line in reference:
    cols = line.split()
    jz, dm, B = map(float, cols[:3])
    refM = float(cols[3].split('=')[1])
    for L, refQ in zip([2,4,6,8], map(float,cols[4:8])):
        row = next(r for r in ours if r['Jz']==jz and r['DM']==dm and r['B']==B and r['L']==L)
        reference_error = max(reference_error, abs(row['Q']-refQ), abs(row['M']-refM))
assert reference_error <= .000051
result = dict(variational_cases=count, largest_bound_violation=bound_violation,
              ferromagnetic_energy_error=endpoint_error, optimized_M_error=max_M_error,
              iterative_solver_energy_error=iterative_error, rounded_survey_max_error=reference_error)
Path('variational_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
