"""Item 1: entry-by-entry comparison of the referee rebuild (ref_bonds.py) with the
author's bond9.bond / bond10.bond at random generic parameter points."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.dirname(HERE))
import ref_bonds as rb
import bond9, bond10

assert bond9.LET == ["L0", "L1", "R0", "R1", "W", "A", "B", "XA", "XB"], bond9.LET
assert bond10.LET[-1] == "V"

rng = np.random.default_rng(20261112)
out = {"h9": [], "h10": [], "h10_with_K_over_J": []}
for trial in range(12):
    par = rb.random_par(rng, mobile=True)
    # author's parameter names
    pa = dict(lam=par["lam"], g=par["g"], Delta=par["Delta"], D=par["D"], nu=par["nu"], f=par["f"],
              phiA=par["phiA"], psiA=par["psiA"], psiB=par["psiB"], td=par["td"],
              J=par["J"], t=par["t"], mu=par["mu"], K=par["K"])
    w, o, bad = bond9.bond(pa)
    H9_author = w + 0.5 * (np.kron(np.diag(o), np.eye(9)) + np.kron(np.eye(9), np.diag(o))) + par["K"] * bad
    H9_ref = rb.h9(par)
    out["h9"].append(float(np.abs(H9_author - H9_ref).max()))
    H10_author = bond10.bond(pa)
    H10_ref = rb.h10(par)
    diff = np.abs(H10_author - H10_ref)
    out["h10"].append(float(diff.max()))
    if trial == 0:
        idx = np.argwhere(diff > 1e-12)
        names = rb.NAMES
        print("h10 mismatching entries at trial 0:")
        for (r, c) in idx:
            print("   <%s %s| h |%s %s>  author %.6f  referee %.6f   (J=%.4f K=%.4f)" % (
                names[r // 10], names[r % 10], names[c // 10], names[c % 10],
                H10_author[r, c].real, H10_ref[r, c].real, par["J"], par["K"]))
    # hypothesis: the code uses K where the README formula says J*K
    par2 = dict(par); par2["K"] = par["K"] / par["J"]
    out["h10_with_K_over_J"].append(float(np.abs(H10_author - rb.h10(par2)).max()))

print("max |h9_author - h9_referee| over 12 points:", max(out["h9"]))
print("max |h10_author - h10_referee| over 12 points:", max(out["h10"]))
print("same with referee K -> K/J (README's J*K read as the code's K):", max(out["h10_with_K_over_J"]))
# count off-diagonals at a generic point
off = H9_ref - np.diag(np.diag(H9_ref))
print("nonzero off-diagonal entries of h9:", int((np.abs(off) > 1e-14).sum()))
off10 = H10_ref - np.diag(np.diag(H10_ref))
print("nonzero off-diagonal entries of h10:", int((np.abs(off10) > 1e-14).sum()))
json.dump(out, open(os.path.join(HERE, "compare_bonds.json"), "w"), indent=1)
