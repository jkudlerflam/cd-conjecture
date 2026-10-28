# C17 progress

Updated 2026-10-27. Work is confined to this folder. No network or large computations.

## Scope and inputs

Read the project STATE.md, ../../specs/_d2_background.md, ../../STATUS.md,
../../specs/C17_d2_mechanisms.md, REPORT.md sections 3.4 and 8, C6's
parametrization, C7's kinetic selector, and C8's qubit FF, boundary-rank,
exact-encoding and symmetry-rank results. There was no existing local progress.
C15 and C16 PROGRESS.md did not exist at the initial check.

The task is a mechanism audit, not a proof of the full qubit conjecture.
Use [proved], [assumed], [heuristic], and [numerical] explicitly. Literature
identifications supplied by the repository are assumptions here, not new checks.

## Findings so far

1. The FF witness K_FF(2)=4 excludes longest-zero-segment selection only.
   A new elementary observation excludes a literal qubit version of C7's
   conserved one-site marker: [H,M_i]=0 for every site and a nontrivial qubit
   projector forces a diagonal binary classical Hamiltonian. Such a chain has
   an energy-minimizing cycle of period at most two. It cannot run a walker.
   Mobile or extended markers, frustrated selection and emergent block codes
   are not excluded.
2. The embedding-index ceiling is conditional. Exact microscopic SU(2)
   reduces the qubit bond to the Heisenberg family; emergent SU(2) has no
   established level ceiling. LSM matching is a congruence, not an upper bound.
3. Exact U(1) reduces the entire bond class to twisted XXZ. The two-body t-V
   equation will distinguish high-momentum bound states from a stable dilute
   pair fluid. Away from U(1), conversion relevance needs the full Luttinger
   matrix and momentum/locality selection rules.
4. Completed Abelian occupation-code partons give c<=1 conditionally. Discrete
   gauge projection and neutral cosets are not covered by this count.
5. Literal even/odd microscopic decoupling is impossible with nonzero NN
   connected interactions. Dual sublattices can occur, e.g. two dual Ising
   chains with total c=1. Emergent decoupling is a separate question.

## Checks completed

checks.py and checks_output.txt now contain reproducible numpy/scipy checks.
All passed. Four-site linear algebra leaves only Jzz,Bz when a one-site marker
is conserved, and only the three fields for literal even/odd decoupling.
Forty random classical binary bonds on rings L=2..10 respect the period-two
minimum-cycle bound. Three two-body bound energies match the half-line formula
to 1.1e-12 or better. The weighted four-parton qubit code has affine rank one.

PBC full-spin ED at L=6,8,10,12,14 calibrates two estimators. XX gives c_S=1.02771
at L=14 and c_Cas=0.999204 from a five-size E/L fit including L^-4. Critical
Ising gives 0.502011 and 0.499976. These are controls, not candidate scans.

C15 and C16 now have progress files, read after the checks. C15 independently
has the same K>=1/2 counterqualification and atom-pair conversion formula;
leave its detailed Z2 phase analysis to that run. C16 studies compactness,
dilute magnons, blockade constraints and operator counts; do not duplicate its
parameter scan or attempt an exhaustion theorem here.

## Completed report

RESULT.md now covers all requested mechanisms, with tagged derivations and a
specific proposed numerical test for each. Its final section ranks the open
routes. No counterexample or full exclusion was obtained. The strongest new
elementary exclusion is the conserved one-site qubit marker, alongside a
physical-tensor-factor proof of the sublattice obstruction. Neither extends
to arbitrary emergent degrees of freedom.

Important qualifications retained in the report:

- A general period-p quantum cell has no established one-mode ceiling.
- Exact microscopic SU(2) is Heisenberg; an emergent SU(2) has no justified
  on-site embedding-index ceiling. Standard spin LSM parity matching allows
  arbitrarily high odd levels under its usual action and cannot upper-bound c.
- Conversion dimensions require locality, momentum and the full Gaussian
  covariance matrix. Even K_a=K_p=1/2 gives Delta_conv=5/2, not relevance.
- Parton subtraction applies only with complete constraints and a physical
  coset identification; discrete gauging does not itself reduce c.
- A nonlocal duality can produce two decoupled Ising sublattices from an
  allowed NN qubit chain, total c=1. This does not contradict physical
  sublattice nondecoupling.
- C15's latest singular pairing-only corner and C16's mixed linear/cubic
  free touching are cited as complementary warnings about effective range
  and nonconformal endpoints. Their new computations were not rerun here.

## Handoff

All requested deliverables are present: RESULT.md, checks.py,
checks_output.txt, PROGRESS.md and local STATE.md. The ranked first target is
an isolated additional neutral critical sector beside a c=1 mode, including
an emergent SU(2)_2 possibility. No global exhaustion theorem was assumed.

No parent STATE.md or other file outside this folder was edited, as explicitly
required by the task. No network, delegation, DMRG, or large computation was
used. There is no unfinished computation. Further work requires selecting a
new microscopic candidate or proving a new dynamical bound, rather than
extending the existing rank/anomaly arguments by assertion.
