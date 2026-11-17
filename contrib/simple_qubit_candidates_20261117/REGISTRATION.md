# Bounded simplicity-first diagnostic campaign

Written 2026-11-17 before new DMRG outputs in this contribution.

We retain the C190/C195 symmetric-line production point at V=40,
kappa=1, hz=-19.976047250613814 as the evidence benchmark. Its field is
not the rounded STATUS value or the untuned models.py reference field.
The full onsite field acts on both endpoints of an open chain.

The two initial low-hierarchy choices are V=8 and V=12, kappa=1. C177's
old bands and open/infinite results motivate these choices, while its
convergence failures and unresolved V=8 FES reading prevent treating
them as verified phases. We first test small rational a near 1, and
permit a bounded reduction to the existing dilute-gas window near
0.95--0.97. This is sequential density bracketing, not a Cartesian scan.
Start with L=32, both parity sectors; do not tune a separately at each L.

Use the literal C190 parity solver in its real staggered computational
gauge, retaining the original common physical Pauli basis for all counts.
Validate the literal matrix and full endpoint fields against ED at L=8
in both sectors; validate entropy and observables, not energy alone.
Benchmark a length-64, chi-128 allocation and fixed local-solve work
before extending runs. Run at most two DMRG jobs simultaneously, each
with one BLAS thread. Keep checkpoints outside Git.

Pilot stages use chi=16,32 with a mixer followed by chi=48,64 without it.
For selected points use chi=64,96,128 and lengths 32,64; extend the most
informative case only if convergence/runtime justify it. At least three
sweeps per unmixed stage and at most six; stopping requires delta E<2e-9,
max delta Sbar over W=1,6,12<5e-4, delta n_mid<2e-4 and local residual<1e-8.
Every limit-hit run remains flagged. A last-stage entropy change above
0.01 or last discarded weight above 1e-7 flags chi convergence. A low
residual alone is insufficient.

Controls are one XX chain and two decoupled XX legs at matched sizes,
near density 0.19, with the chemical potential midway between finite-chain
levels. Compare with exact free-fermion entropies, using each leg's own
Jordan--Wigner map. Use the same solver, windows and comparable chi
stages. Report both parity sectors for candidates, lower-energy branch,
all parity pairings and the density drift. Record central density,
adjacent occupation, valley weights and connected separation-two pair
correlations. Do not fit a Luttinger exponent from a few short distances.

Selection metric: independent real parameters (also after scale fixing),
nonzero physical Pauli monomials, raw coefficient spread with Delta=1,
observed tolerance in a/kappa, and strength of actual convergence/size
evidence. No entropy reading at L<=64 establishes a thermodynamic CFT.
No new point is called c>1 established without stronger, converged
infrared evidence. A failure boundary or resource-limited unresolved
point is an acceptable reported result; never prefer a larger apparent
fit over calibrated evidence. The epsilon=10^-3 sixth-order construction
remains a separate mechanism study and is not a practical candidate here.
