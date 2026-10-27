# C12 progress

The task is to find a simple one-site translation-invariant nearest-neighbor Hermitian chain with total central charge c > d - 1, preferably d <= 64. This run is closed book. No pre-existing files are read or listed, and no external sources are consulted. Work and generated files stay in the current directory. The supplied facts are taken as input.

## Initial strategy

Screen frustration-free cells, squeezed-space models, emergent constraints, and spontaneously formed ladders before committing to numerical work. A successful construction must have an explicit two-site interaction, a controlled relativistic infrared theory, and two calibrated central-charge estimators. Degeneracy and a quadratic ferromagnetic mode are not central charge. Negative results will be reported as such.

The request to remain closed book supersedes the usual project-memory reading and updating workflow. PROGRESS.md and RESULT.md will carry this run's record.

## First screening

1. A frustration-free ferromagnetic SU(q) cell has symmetric ground space Sym^b(C^q), dimension binomial(b+q-1,q-1). Projection of an end generator gives the collective generator divided by b. This produces a large representation, not b independent flavors. An SU(q)_b WZW point would be attractive because c=b(q^2-1)/(b+q), but the projected first-order bond is only bilinear and does not supply the integrable higher-degree interaction. This is a possible target, not a constructed counterexample.
2. A fixed injective MPS parent has a cell ground space bounded by D^2 once the length exceeds its injectivity scale. Its two edge spaces do not grow with cell length. Colored-path constraints can have growing degeneracy, but boundary access is through stack operations, not independent gapless channels.
3. A rigid rung carrying k independent q-state chains costs d=q^k and gives c=k(q-1), which never exceeds d-1. Unrolling the rung into a long cell moves the problem to transmitting k independent channels through its single boundary bond.
4. Large entropy from an SU(q) ferromagnet, extensively degenerate constraints, or a quadratic band edge cannot be counted as a relativistic central charge. These need explicit z tests.

No candidate has yet passed the screen. The main remaining analytic possibility is an inexpensive large-representation cell tuned to a higher-level WZW point. Its missing ingredients are period selection and independently tunable projected interactions.

## A direct dilute-gas route at d=3

Consider local occupation states 0,1,2, with a conserved total occupation N. A general small family has atomic hopping 10 <-> 01, pair hopping 20 <-> 02, and pair conversion 11 <-> 20,02. Three stable incommensurate minima, one atomic and two molecular, would potentially give three Luttinger modes. Ordinary two-atom conversion is oscillatory when molecular minima occur at nonzero momenta; the momentum-neutral four-atom conversion must also be tested by its scaling dimension. Stability against larger bound clusters and phase separation is essential.

Created pair_screen.py, an exact two-particle relative-coordinate diagonalization. In 3,000 random parameter samples for equal-sign, reflection-even conversion, no interior pair-band global minimum below the two-atom threshold was found. The relative coordinate is truncated at 64 sites, and K is sampled at 65 points. This is a negative numerical screen, not a general exclusion theorem. The first launch failed because this environment provides python3 rather than python; the rerun completed.

An odd-parity conversion (opposite signs for 11 -> 20 and 11 -> 02) changes the pair-channel form factor from cos(K/2) to sin(K/2). It is Hermitian and one-site translation invariant. A second two-body scan is running for this case.
