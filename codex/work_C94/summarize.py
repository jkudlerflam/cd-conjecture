"""Summarize only C94's generated checks."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
c=json.loads((HERE/'canonical_results.json').read_text())
t=json.loads((HERE/'thermal_results.json').read_text())
s=json.loads((HERE/'coefficient_scaling.json').read_text())
lines=['# C94 numerical checks','','[Numerical] These are finite logical vacuum-clock calculations, with quantum contact sources and distributed marker energies. Full C88 auxiliary traces are treated analytically in PROOF.md, not diagonalized here.','',
 f"Canonical comparisons: {len(c['rows'])}. Thermal comparisons: {len(t['rows'])}. Maximum joint matrix dimensions: {max(r['dimension'] for r in c['rows'])} canonical and {t['metadata']['max_dimension']} thermal.",'',
 f"Largest canonical eigenpair residual: {max(r['residual'] for r in c['rows']):.5g}. Largest bath linear-solve residual: {max(r['linear_residual'] for r in c['kernels']):.5g}. All canonical defect costs are positive.",'',
 'The table reports max(0,-Delta R)/(J W). The canonical case is (30) against (15,15), N=30, v=1. The thermal case is (2,2) against (4), N=4, s/t=0.65, T=J.','',
 '| J/t | Canonical observed epsilon | Thermal observed epsilon |','| :-- | --: | --: |']
for j in [.04,.02,.01,.005,.0025]:
 a=next(r for r in c['rows'] if r['J']==j and r['N']==30 and r['lengths']==[30])
 b=next(r for r in t['rows'] if r['J']==j and r['T_over_J']==1)
 lines.append(f"| {j} | {a['epsilon_observed']:.9g} | {b['epsilon_observed']:.9g} |")
lines+=['',f"Largest thermal residual-difference cutoff error bound: {max(r['delta_R_cutoff_error_bound'] for r in t['rows']):.5g}.",'',
 '[Numerical] At the separate free vacancy crossing, the logical response grows as 1/T.','',
 '| T/t | D_wrong-D_good, t=1 | T(D_wrong-D_good) |','| :-- | --: | --: |']
for r in t['crossing']:lines.append(f"| {r['T']} | {r['D_difference']:.10g} | {r['T_D_difference']:.10g} |")
lines+=['',f"Analytic limit for the last column: {t['metadata']['crossing_slope']:.12g}.",'',
 '[Numerical] Half-density coefficient checks use the exact bath resolvent. They do not form the larger joint interacting matrices at N=9.','',
 '| N=v | Defect | W | (B_wrong-B_good)/W |','| :-- | :-- | --: | --: |']
for r in s:lines.append(f"| {r['N']} | {r['label']} | {r['W']} | {r['delta_B_over_W']:.12g} |")
lines+=['','[Scope] The data support decreasing normalized errors on these finite systems. They do not bound a thermodynamic supremum or an all-order remainder uniform in defect length. The b=15 joint calculations have high occupancy. No half-density b=15 interacting trace is claimed.','',
 '[Reproduction] Run python3 -B checks.py, python3 -B coefficient_scaling.py, python3 -B crossing_certificate.py, then python3 -B summarize.py from this folder. The inherited modules are local copies of the read C89 implementations, with the local import adjusted.','']
(HERE/'CHECKS.md').write_text('\n'.join(lines))
