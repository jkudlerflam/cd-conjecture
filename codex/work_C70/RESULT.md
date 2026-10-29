# C70 result

[assumed, unresolved target] C60 (17) remains neither proved nor refuted for the full single-channel routing. I did not construct the requested simultaneous bounded-load expansion or exhibit a legal nonvanishing counterexample tree.

[proved, local extension] PROOF.md Sections 2-4 remove C64's independence restriction for one quartic substitution. The child operator acts on the entire exterior function, including affine ancestor Taylor fields, weighted covariance factors, and determinant rows or columns. Every derivative term has a source/sink record covering the child cut. Consumed weights concatenate paths as in C61. A derivative reaching an exterior field ends at that field's later contraction scale and gives an additional favorable gain.

[proved, representative formula] With a child source at t, a line at s carrying b ancestor moments, and their sinks at h, the missing exterior-field derivative has scale factor

    2^[(h-t)+b(h-s)] = 2^(h-s) 2^[(s-t)+b(h-s)].

Its cut counts are one on the child interval and b+1 on the parent interval. It does not create a deficient cut. At b=2 it can produce a third exterior derivative, so this local result alone does not establish the uniform load bound.

[proved, surviving-port decomposition] Section 5 splits a quartic remainder into components carrying surviving external momentum factors and a component using only the other independent child momenta. A surviving factor makes a later quartic contact subtraction trivial. At a quadratic ancestor it lowers the residual Taylor order. This is an exact decomposition before ancestor coordinate moments are allocated.

[assumed, remaining gap] Those incoming moments differentiate both the surviving factors and their Taylor-integral kernels. One must reorganize these terms before absolute values and retain their spatial paths throughout the reorganization. The local formulas do not prove that this can be done simultaneously with bounded endpoint moments, field derivatives, coefficient mass and permitted anchor incidence. No global majorant, source estimate or thermal conclusion follows from this run.

[numerical, finite exact-arithmetic checks] `python3 -B checks.py` passes 72 affine product-rule checks, 36 transposition checks, 14 surviving-port decompositions, 30 homogeneous-remainder identities, 837 cut-record checks, and 39 Taylor cancellations. These are local algebra checks, not an exhaustive legal-tree test. Output is in checks.json.

[record] PROOF.md contains the derivations and explicit limitations. PROGRESS.md and STATE.md provide the handoff. All writes are confined to work_C70. No independent review or subagents were used.
