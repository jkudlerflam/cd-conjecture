#!/usr/bin/env python3
"""Keep optimizer logs under 1 MB: retain every final entry and, per start, the last 5 logged
improvements (the trajectory of 'best so far' points is otherwise thousands of lines)."""
import json, glob, os, sys
for f in sys.argv[1:]:
    rows = [json.loads(l) for l in open(f)]
    keep = []
    by_start = {}
    for r in rows:
        if r.get('final'):
            keep.append(r)
        else:
            by_start.setdefault(r['start'], []).append(r)
    for s, lst in by_start.items():
        keep.extend(lst[-5:])
    keep.sort(key=lambda r: (r['start'], r.get('final', False), r.get('eval', 10 ** 9)))
    with open(f, 'w') as g:
        for r in keep:
            g.write(json.dumps(r) + '\n')
    print(f, len(rows), '->', len(keep), os.path.getsize(f) // 1024, 'KB')
