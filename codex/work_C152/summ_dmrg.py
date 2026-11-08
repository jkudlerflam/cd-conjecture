"""Summarize dmrg_pt.py outputs: label, L, c fits (parity 0, 1, averaged) on central 25-75 and 12-88 windows, S_mid, truncation."""
import json, sys
recs = []
for fn in sys.argv[1:]:
    for line in open(fn):
        try: recs.append(json.loads(line))
        except Exception: pass
recs.sort(key=lambda r: (r['label'], r['L']))
for r in recs:
    a = r['c_25_75']; b = r['c_12_88']
    print(f"{r['label']:40s} L={r['L']:3d} c25={a[0]:6.2f} {a[1]:6.2f} {a[2]:6.2f}  c12={b[0]:6.2f} {b[1]:6.2f} {b[2]:6.2f}  Smid={r['S_mid']:.3f} tr={r['trunc']:.1e} dE={r['dE_last']:.1e}")
