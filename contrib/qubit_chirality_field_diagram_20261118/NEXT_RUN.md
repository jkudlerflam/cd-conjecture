# One discriminating next calculation (not launched)

The most useful next calculation is an **operator-resolved chirality
measurement on verified large-xi states near h=1.7465**, including a
momentum-resolved excitation estimate. It must distinguish the actual
same-parity, physical momentum-pi chirality mode from low opposite-parity
levels and from oscillatory boson composites. Another entropy fit alone
does not do this.

## First use the already-produced states

Obtain the collaborator's existing HDF5 states, full NPZ transfer spectra
and solver block logs, with hashes and parameter records. First request
h=1.7465, chi=128,160,200,256 on the pin branch and the available
unpinned/complex comparison. Then use h=1.7405 and 1.750 as neighboring
diagnostic controls. The states were not in this checkout; nothing in
this contribution claims to have contracted them.

Use a separate independent contraction implementation for the first
small-state checks. Verify energy, Z, bond expectations, physical parity,
the real-frame sign maps and transfer-cell length two. Retain the
physical translation/gauge action so momentum pi is not confused with
momentum zero after folding into the two-site cell.

1. Resolve several leading transfer eigenvalues in each parity sector,
   including phases and left/right eigenvectors. Project local Z, literal
   bond, B2 and kappa onto them. Mark the near-unit cat eigenvalues
   explicitly and check that the identified modes are stable when the
   number of computed eigenvalues increases. Do not assign sector xi0
   to the Ising field solely because kappa is parity even.
2. Measure both anchors of the connected kappa correlator, subtracting
   the complex disconnected part before averaging. Retain r=8--1024,
   the exact signs and phases, and separate the unmodulated staggered
   component from cos(Q r) components. Fit several windows, for example
   32--128, 48--192 and 64--256 when those windows lie below the relevant
   operator correlation length. Fit Q from data/transfer phases and allow
   exponents to vary; do **not** impose eta=1/4. Compare a pure Ising
   component with Ising-times-boson composite alternatives.
3. In the same state obtain same-parity variational excitation energies
   and kappa spectral weights near physical momentum pi, and comparison
   boson-sensitive excitations near their measured soft momentum. Use
   both signs of small momentum offsets, e.g. 0.005,0.01,0.02,0.04 pi.
   Fit velocity and residual mass jointly, allowing unequal velocities.
   Repeat the smallest set at two chi values; distinguish a positive
   extrapolated speed from a finite-chi mass or a velocity-collapse fit.

An unmodulated eta near 1/4 with stable chi/window behavior and a distinct
positive-speed chirality branch would be independent support for the
boson x Ising interpretation. A second Gaussian tower or a different
operator exponent would change that interpretation. Failure to converge
is an inconclusive result, not evidence against Ising. Even successful
finite-scale tests would leave sufficiently weak first-order behavior open.

## Resources and stopping rules

This is a proposed protocol, not authorization to restart production.
Start with one worker, one BLAS thread and an 8 GiB memory ceiling.
Benchmark transfer/operator contractions at chi<=128 with a **15-minute
wall limit**, then allow at most **60 minutes total** for measurements of
existing states, including the excitation pilot. Those are ceilings, not
claimed runtime predictions. Save partial measurements and terminate at
the cap; do not self-resubmit. If the excitation pilot cannot complete
within the remaining budget, stop and report its measured cost.

For a state to enter an exponent or velocity comparison, require its
supplied logs (or an explicitly authorized future reoptimization) to show
successive mixer-free block changes below 1e-11 in energy density,
5e-6 in entropy, and 5e-6 in **both** chirality and B order, plus stable
operator-resolved transfer lengths to 1%. Repeat at least one independent
initial branch; a persistent order-orientation difference must be retained
even when energies agree to 1e-8. Check alternative unit cells in any
future state production before declaring a global minimum.

If missing states/logs or failed convergence force new MPS production,
first submit a concrete small proposal: the fixed h=1.7465 point,
chi=128,200,256, two starts (chirality seed and real/bond seed), with
1 core, 8 GiB and a **6 CPU-hour total cap**. The original submission's
1--6 CPU-hour ladder estimate is not an independent benchmark; stop after
the first stage if measured cost cannot fit the remaining cap. Do not
extend to chi=512 or scan fields without a new resource decision.

The lower boundary near 1.668 and the upper boundary above 1.758 deserve
later dedicated scaling. They are secondary to resolving whether this
very simple fixed bond actually supports the proposed extra critical
sector. No existing negative search on another Hamiltonian excludes it.
