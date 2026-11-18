# Bounded three-flavour continuation

Keep physical d=9 and the frozen marker sector of the preceding audit.
Use only the literal controls N, N_L, independent LL/RR swaps, and
P_M tensor |L1><L1|. First verify full-code projections; a free symmetric
logical onsite matrix is not an admissible substitute.

Benchmark the n=3 finite pair coefficient before any search. The first
benchmark took .038 seconds, one BLAS thread, including the generated
quartic on the 64-dimensional pair code. No DMRG or many-body size scan.

Use bounded algebraic searches in the five-dimensional real tridiagonal
space. Retain a reproducible finite search budget and simplify any seed
to rational controls where the root and flow margins permit. Compute
all quartic momentum channels in the original orthonormal band basis.
Include the coupled three-channel one-loop flow, rather than applying
three independent two-channel tests. Report an obstruction or an
undecided step if the available calculation does not establish stability.

Every conclusion retains finite-K, port-order, sector-selection and
infrared assumptions. Even a successful n=3 point is neither an all-n
family nor a c>d-1 result. All changes for this round stay in this folder.

Completed budget: 108 local solves, in twelve-start batches, with at
most 1500 evaluations for the first root design and 500 for each later
design. The generic-channel flow test excludes 48 points with extra or
nearly resonant momentum channels; 60 generic points reach a prescribed
twentyfold-growth threshold. A finite rounding ladder stops at relative
denominator 5000 and supplies the rational primary point. Exact-root
certification and the analytic growth-cone argument for that point do
not depend on the optional search being exhaustive.
