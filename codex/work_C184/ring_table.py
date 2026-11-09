import json, sys, collections, math
f = sys.argv[1]
d = collections.defaultdict(dict)
for l in open(f):
    r = json.loads(l); d[(r['Delta'], r['V'])][r['L']] = r
Vs = sorted({k[1] for k in d}); Ds = sorted({k[0] for k in d})
def diag(r):
    lev = r['levels']; e0 = lev[0]
    dpar = abs(r['E0_par'][0] - r['E0_par'][1])
    # first gap above the two-fold (if ordered) or one-fold
    g1 = lev[1] - e0; g2 = lev[2] - e0
    return dpar, g1, g2, r['S']
print("Delta  V   | L*dpar(12) L*dpar(16) | L*g1(12) L*g1(16) | L*g2(16) | S12 S16")
for D in Ds:
    for V in Vs:
        if 16 not in d[(D, V)] or 12 not in d[(D, V)]: continue
        a = diag(d[(D, V)][12]); b = diag(d[(D, V)][16])
        print(f"{D:4.1f} {V:5.1f} | {12*a[0]:8.3f} {16*b[0]:8.3f} | {12*a[1]:8.3f} {16*b[1]:8.3f} | {16*b[2]:8.3f} | {a[3]:.3f} {b[3]:.3f}")
