#!/usr/bin/env python3
"""Summarize completed periodic-ring records without interpreting gaps as CFT tests."""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
source = HERE / 'preferred_rings.jsonl'
records = [json.loads(line) for line in source.read_text().splitlines()]
rows = []
for name in sorted({r['name'] for r in records}):
    rs = sorted((r for r in records if r['name'] == name), key=lambda r: r['L'])
    assert len({r['L'] for r in rs}) == len(rs)
    for a, b in zip(rs, rs[1:]):
        rows.append(dict(name=name, sizes=[a['L'], b['L']],
            parities=[a['selected_parity'], b['selected_parity']],
            c_half_size=3*(b['entropy_half']-a['entropy_half'])/math.log(b['L']/a['L']),
            profiles=[a['c_profile'], b['c_profile']],
            max_ritz_residual=max(x for r in (a,b) for s in r['sectors'] for x in s['residuals']),
            momentum_resolution_pass=all(s['momentum_resolution_pass'] for r in (a,b) for s in r['sectors'])))
out = dict(source=str(source.relative_to(HERE)), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
           note='PBC c_eff=3 Delta S_half/log(L2/L1). Finite sizes do not certify a CFT or exclude dilute sectors.', rows=rows)
(HERE / 'preferred_ring_summary.json').write_text(json.dumps(out, indent=2)+'\n')
lines = ['# Independent periodic-ring check', '',
    'Both attractive models were solved in the full microscopic qubit Hilbert space,',
    'separately in both parity sectors. Translation orbits of returned ground',
    'vectors were completed before selecting entropies. The lowest momentum',
    'index in the ground manifold is the fixed selection rule. All selected',
    'ground states here have even parity and momentum zero.', '',
    '| Point | Sizes | Half-chain size slope | Profile slopes |',
    '|---|---:|---:|---:|']
for r in rows:
    a,b=r['sizes'];p,q=r['profiles']
    lines.append(f"| {r['name']} | {a}→{b} | {r['c_half_size']:.6f} | {p:.6f}→{q:.6f} |")
lines += ['',
    'The size estimator is `3 ΔS_half / log(L2/L1)` for periodic chains.',
    'The earlier same-code L16→20 controls give approximately 0.993 (XX),',
    '0.499 (Ising), and 1.940 (two independent XX chains, a range-two control).',
    'These preferred attractive rings do not independently support an',
    'above-one entropy coefficient at these lengths. They also do not rule',
    'out longer crossover scales or a dilute additional sector. The open-chain',
    'entropy excess must survive larger lengths and order diagnostics.', '',
    'Raw Lanczos spacings remain explicitly distinguished from gaps above a',
    'certified complete ground manifold; completeness of the low spectrum is',
    'not established. No gap or Casimir-energy veto is applied.', '',
    'The exact field is zero. The rational point has coefficients',
    '`(Jx,Jy,D,Jz)=(14,-5,11,-12)/24`, with raw spread 2.8.',
    'The reference is `phi=1.2, Delta=.8, V=-2, t=1` in the stated',
    'fermion convention. It has raw spread approximately 2.656.', '',
    'Reproduce with `prepare_rings.py`, the registered `ed_followup.py` command',
    'in the continuation protocol, and `analyze_rings.py`. Per-sector residuals,',
    'raw/resolved entropies, coupling metrics and source hashes are retained',
    'in `preferred_rings.jsonl` and `preferred_ring_summary.json`.', '']
bridge_path=HERE/'bridge_rings.jsonl'
if bridge_path.exists():
    rs=sorted([json.loads(line) for line in bridge_path.read_text().splitlines()],key=lambda r:r['L'])
    bridge_rows=[]
    lines += ['## R=8.75 bridge', '',
        'The same registered calculation for `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`',
        'also selects even parity and momentum zero at every completed length.', '',
        '| Sizes | Half-chain size slope | Profile slopes |', '|---:|---:|---:|']
    for a,b in zip(rs,rs[1:]):
        c=3*(b['entropy_half']-a['entropy_half'])/math.log(b['L']/a['L'])
        row=dict(sizes=[a['L'],b['L']],c_half_size=c,profiles=[a['c_profile'],b['c_profile']],
            max_ritz_residual=max(x for r in (a,b) for s in r['sectors'] for x in s['residuals']),
            momentum_resolution_pass=all(s['momentum_resolution_pass'] for r in (a,b) for s in r['sectors']))
        bridge_rows.append(row)
        lines.append(f"| {a['L']}→{b['L']} | {c:.6f} | {a['c_profile']:.6f}→{b['c_profile']:.6f} |")
    lines += ['', 'These small rings are close to the one-mode control. The larger open-chain',
        'entropy excess therefore has no independent periodic-ring confirmation',
        'yet. A longer crossover length is possible; a c>1 CFT is not established.', '']
    (HERE/'bridge_ring_summary.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(bridge_path.read_bytes()).hexdigest(),rows=bridge_rows),indent=2)+'\n')
(HERE / 'RINGS.md').write_text('\n'.join(lines))
print(json.dumps(out, indent=2))
