#!/usr/bin/env python3
"""Compact live status; only completed records are treated as measurements."""
from collections import Counter
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
for folder in ('validation', 'attractive_scan/results', 'bridge_scan/dmrg'):
    records = []
    for path in sorted((HERE / folder).glob('*.json')):
        r = json.loads(path.read_text())
        if 'history' in r:
            records.append((path, r))
    print(folder, dict(Counter(r.get('status') for _, r in records)))
    for path, r in records:
        if r.get('status') == 'done':
            continue
        h = r.get('history', [])
        s = h[-1] if h else {}
        keys = ['sweep_global', 'chi', 'E_physical', 'S_fractional', 'trunc', 't']
        print(' ', path.stem, {k: s.get(k) for k in keys})
