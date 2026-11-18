"""Negative control for ring_invariants.py: a bond with one extra, deliberately wrong
term must produce violations, otherwise the checker is vacuous."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_bonds as rb, ring_invariants as ri
rng = np.random.default_rng(7)
par = rb.random_par(rng, mobile=True)
out = {}
# (a) a second head rule (R_a, W) <-> (W, L_a): should break F and the boundary pattern
H = rb.h9(par).copy()
H[rb.R0 * 9 + rb.W, rb.W * 9 + rb.L0] += 1; H[rb.W * 9 + rb.L0, rb.R0 * 9 + rb.W] += 1
out["second_head_rule"] = ri.run(9, 4, H, ["F", "Q", "heads", "occ", "bmask"])["viol"]
# (b) a context-free bit flip with a token: breaks nothing in F but breaks Q
H = rb.h9(par).copy()
H[rb.A * 9 + rb.R0, rb.A * 9 + rb.R1] += 1; H[rb.A * 9 + rb.R1, rb.A * 9 + rb.R0] += 1
out["charge_violating_flip"] = ri.run(9, 4, H, ["F", "Q", "heads", "occ", "bmask"])["viol"]
# (c) C84-style contextual walker creation (L_a, R_0) <- (R_0, R_a): creates an (L,R) pair
H = rb.h9(par).copy()
H[rb.L0 * 9 + rb.R0, rb.R0 * 9 + rb.R0] += 1; H[rb.R0 * 9 + rb.R0, rb.L0 * 9 + rb.R0] += 1
out["contextual_walker"] = ri.run(9, 4, H, ["F", "Q", "heads", "occ", "bmask"])["viol"]
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, "ring_invariants_control.json"), "w"), indent=1)
