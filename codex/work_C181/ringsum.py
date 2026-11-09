import json, sys
for f in sys.argv[1:]:
    for l in open(f):
        r = json.loads(l)
        g = r['gs']
        print(f"{f[8:-4]:10s} phi={r['phi']:.4f} gsK={g[0]:3d} E={g[1]:.5f} M={g[2]:+.2f} |", " ".join(f"(K{x[0]},{x[1]-g[1]:.4f},M{x[2]:+.1f})" for x in r['low'][1:5]), f"{r['secs']:.0f}s")
