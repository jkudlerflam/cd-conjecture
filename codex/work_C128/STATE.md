# C128 local handoff

2026-11-06.

[Proved] New XX estimate: log Ztilde_L<=log 3+beta v[csc(pi/L)-L/pi], from C125's positive crossed product and convex midpoint quadrature. The integral is I_L=pi/[2 sin(pi/L)]-L/2. R is at most twice this expression. With M=v/2 the all-even-L>=4 ratio coefficient is 16sqrt(2)-64/pi=2.2555842822, and including L=2 gives 8-16/pi=2.9070418211. A literal but weaker scalar Fock product has S=2, V'=4pi M, K=log 3-2 log(1-exp(-1)). It is not a frequency-preserving spectral injection.

[Proved] The explicit nonnegative column construction uses alternating K_0,K_1 and a two-column matrix B=K_0 K_1. The positive operator T=sqrt(K_0^(1/2) K_1 K_0^(1/2)) reproduces the even-ring traces, but is not asserted entrywise nonnegative. Continuous histories have shared event times across neighboring columns.

[Proved] A column-event count generating function is bounded by exp[beta(z-1)] for all requested XXZ models. This does not justify independent temporal clocks. The latter comparison is singular because a paired jump has weight order t, versus order t^2 for independent column events. Conditional bridge rates are unbounded. All-frequency linear crossed decay fails already for XX, inherited from C125.

[Proved] The unchanged XX normalized-trace bound fails for actual Heisenberg: at L=4 its beta slope is 3-4 log 2=0.22741128, versus sqrt(2)-4/pi=0.14097402 in the XX expression at the same M=1/2.

[Open] A general scalar Fock majorant, and finite C(d), are not settled. PROOF equation (18) states the missing interacting partition-function comparison with the valid linear-low-frequency/logarithmic-high-frequency rate and a controlled zero-frequency prefactor. Equation (19) proves its sufficient scalar consequence. No proof of the many-state comparison and no lattice counterexample to every possible set of majorant constants were obtained.

[Numerical] checks.py reproduces 192 small-Trotter cases and 35 exact XX cases. Maximum relative contraction error 5.58e-12. Explicit tails are retained for the slowly convergent L=2 product. Generous interacting candidate passes samples, unchanged XX estimate fails in 24 samples. CHECKS.md records scope, numerical failures and the corrected overflow in the first tail computation.

[Work record] All writes stayed in this folder. Project STATE.md is outside the authorized write area and was not updated. run_typhon.slurm is ready and unsubmitted. No network shell, cluster connection, sub-agent, independent review, or concurrent-run result was used. C127's specification was read, but its work was not read or awaited.
