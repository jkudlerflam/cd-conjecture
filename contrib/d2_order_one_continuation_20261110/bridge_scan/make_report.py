#!/usr/bin/env python3
"""Render a compact, reproducible report of complete results only."""
import json
import math
from pathlib import Path
from summarize import HERE, summarize

out=summarize();name='b_r5_2_d1_10_q0_hm7_8'
chi_checks=[]
for snap in sorted((HERE/'validation_snapshots').glob('*.json')):
    a=json.loads(snap.read_text())
    current=HERE/'dmrg'/f'{a["name"]}_L{a["L"]}_P{a["parity"]}_s{a["seed"]}.json'
    if not current.exists():continue
    b=json.loads(current.read_text())
    if b['status']!='done' or b['chi']<=a['chi']:continue
    chi_checks.append(dict(name=a['name'],L=a['L'],parity=a['parity'],seed=a['seed'],chis=[a['chi'],b['chi']],
        delta_E=b['E_physical']-a['E_physical'],delta_S_fractional=b['S_fractional']-a['S_fractional'],
        delta_S_central3=b['S_central3']-a['S_central3'],both_sweep_gates=a['convergence']['sweep_converged'] and b['convergence']['sweep_converged']))
chi_checks.sort(key=lambda x:(x['L'],x['parity'],x['chis'][0]))
(HERE/'chi_convergence.json').write_text(json.dumps(chi_checks,indent=2)+'\n')
rows=sorted([r for r in out['rows'] if r['name']==name],key=lambda r:r['L'])
pairs=[p for p in out['pairs'] if p['name']==name]
lines=['# Comparable-coupling bridge continuation','',
'The promising finite-size candidate has exact couplings `(Jx,Jy,Jz,Dz,hz)=(1/2,-1/2,5/8,1/10,-7/8)`.',
'Its literal coefficient spread is R=8.75, in the registered exploratory range 3<R<=10. Multiplying the Hamiltonian by 40 gives integer coefficients `(20,-20,25,4,-35)`. It has six expanded Pauli strings, four grouped terms, and the exact relation Jx+Jy=0. This cancellation is disclosed; no basis change is used to reduce the advertised hierarchy. Exchange singular values are 5/8, 3/5, 2/5 (spread 25/16), and the antisymmetric-to-symmetric Frobenius ratio is 0.149854.','',
'The preregistered closest-density rule selected this point from the original L24 grid. The exact same field and couplings are used at every larger size; there is no size-dependent retuning. A surrounding floating phase and robustness to q=Jx+Jy perturbations have not been established.','',
'This is a finite-size numerical lead, not a demonstration of an asymptotic c>1 or conformality. The parent\'s independent small periodic chains currently look near c=1, so the geometry dependence and remaining size corrections must stay explicit.','',
'| L | chi | central density | density range | fractional S | central-three S | profile c (1/4,3/8) | far XX RMS | both sector sweep gates |',
'|---|---|---|---|---|---|---|---|---|']
for r in rows:
    lines.append(f"| {r['L']} | {r['chi']} | {r['density']:.8f} | {r['density_bulk_range']:.8f} | {r['S_fractional']:.9f} | {r['S_central3']:.9f} | {r['profiles']['0.25']['c']:.5f}, {r['profiles']['0.375']['c']:.5f} | {r['long_correlation_rms']['X']:.8f} | {r['both_sectors_converged']} |")
lines+=['','The chi column refers to the retained lower-energy sector. Both-sector bond dimensions and source records are explicit in `summary.json`.','',
'| Size pair | fractional-window c slope | central-three c slope |','|---|---|---|']
for p in pairs:lines.append(f"| {p['sizes'][0]}→{p['sizes'][1]} | {p['c_fractional']:.6f} | {p['c_central3']:.6f} |")
if chi_checks:
    lines+=['','| L, parity | chi step | delta E | delta fractional S | delta central-three S |','|---|---|---|---|---|']
    for c in chi_checks:lines.append(f"| {c['L']}, {c['parity']} | {c['chis'][0]}→{c['chis'][1]} | {c['delta_E']:.3g} | {c['delta_S_fractional']:.3g} | {c['delta_S_central3']:.3g} |")
pending=[r for r in rows if r['running_upgrades']]
if (HERE/'largest_size_preservation.json').exists():
    largest=next((r for r in rows if r['L']==192),None)
    if largest is not None and largest['chi']<256 and largest not in pending:pending.append(largest)
if pending:lines+=['','A higher-chi validation is pending completion at '+', '.join('L='+str(r['L']) for r in pending)+'. The tables retain the latest completed state, with snapshot provenance in `summary.json`.']
lines+=['','The correlation origin is at site L/4 and endpoints extend through 3L/4; the far-correlation RMS uses the latter half of that separation array. It is an order diagnostic, not a fitted critical exponent. Sector splittings remain variational diagnostics.']
external64=list((HERE.parent/'independent_sizes').glob(name+'_L64_P*_s17.json'))
if len(external64)==2:
    rr=[json.loads(p.read_text()) for p in external64]
    if all(r['status']=='done' and r['convergence']['sweep_converged'] for r in rr):
        b=min(rr,key=lambda r:r['E_physical']);assert b['couplings']==rows[0]['couplings']
        anchors={r['L']:r for r in rows};a,c=anchors[48],anchors[96]
        low={k:6*(b[k]-a[k])/math.log(64/48) for k in ['S_fractional','S_central3']}
        high={k:6*(c[k]-b[k])/math.log(96/64) for k in ['S_fractional','S_central3']}
        lines+=['','The parent independently checked L=64 with both parities at chi192. The fixed-field adjacent slopes are '
            f"48→64: {low['S_fractional']:.6f}/{low['S_central3']:.6f} and 64→96: {high['S_fractional']:.6f}/{high['S_central3']:.6f} "
            '(fractional/central-three). See [intermediate-size records](../independent_sizes/).']
ring=HERE.parent/'bridge_ring_summary.json'
if ring.exists():
    r=json.loads(ring.read_text())['rows'][-1]
    lines+=['',f"The independent periodic {r['sizes'][0]}→{r['sizes'][1]} entropy slope is {r['c_half_size']:.6f}; "
        f"the L={r['sizes'][1]} profile estimate is {r['profiles'][-1]:.6f}. These small-ring estimates remain near one, while the larger open-chain entropy estimates are higher; the difference in geometry and accessible sizes remains unresolved. "
        'See [periodic summary](../bridge_ring_summary.json).']
lines+=['',
'## Controls and provenance','',
'The unchanged prior C190 full microscopic DMRG wrapper is imported with new output and checkpoint directories. Both exact parities are solved. Repeated XX and symmetric r2 plumbing controls reproduce the prior energy, entropy and density exactly; see `plumbing_controls.json`. The previous contribution is unchanged.','',
'`PROTOCOL.md` predates all output. `promoted.json` freezes the original six-point promotion. `ADAPTIVE_AMENDMENT.md` separately declares two post-screen safeguards and six finer fields. `R3_EXTENSION.md` declares the ten-field r3 scan and bounded density refinement. `LEAD_VALIDATION.md`, `LARGE_SIZE_EXTENSION.md` and `L192_CHI_VALIDATION.md` declare follow-up checks before their output.','',
'The original grid has 36 Hamiltonians plus plumbing controls; two density bisections were added. Five weak single-size promotions selected for L48, six targeted L24 refinements and ten r3 L24 fields remain unrun and explicitly deferred after the search focused on this lead. These five promotions are untested at matched size, not already ruled out. A DEFERRED log entry is not a computed result.','',
'All completed numerical records, including negative and contradictory results, remain in `dmrg/`; `summary.json` includes field-energy concavity checks. `validation_snapshots/` preserves earlier chi JSON records. The old L96 even chi192 NPZ was overwritten and is unavailable; its raw JSON and pre-upgrade diagnostics remain. All eight final lead states have hash-verified JSON/NPZ snapshots in `final_state_preservation.json`. Later preservation is documented in `STATE_PRESERVATION.md`; the parent also created durable Git-ignored workspace copies, verified by [local checkpoint manifest](../local_checkpoint_manifest.json). No cluster job, commit or push was performed by this scan.','',
'## Other completed matched-size points','',
'| Point | fractional c24→48 | central-three c24→48 |','|---|---|---|']
for p in out['pairs']:
    if p['name']!=name:lines.append(f"| `{p['name']}` | {p['c_fractional']:.6f} | {p['c_central3']:.6f} |")
lines+=['','The near-preferred R3.25 point has smooth density but a size slope near one. The r2,D3/16,h=-5/8 entropy peak decreases at L48 on the lower-energy branch; its high-entropy odd sector is higher in energy. Neither is a positive c>1 lead.','',
'![Full initial field scan](field_scan.png)','',
'![Matched-size slopes](size_slopes.png)','']
(HERE/'REPORT.md').write_text('\n'.join(lines))
print(HERE/'REPORT.md')
