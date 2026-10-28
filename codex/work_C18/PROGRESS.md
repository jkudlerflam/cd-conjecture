# C18 progress

## Final handoff

Partial construction completed. RESULT.md contains the full derivation,
explicit rules, all-sector bounds, exact counts for implemented models,
reproduction commands, numerical record, and ranked next steps.

[proved] The C8 full-sweep clock admits two permanently idle qubit ports
at dimension 15w+4. Its own hopping energy selects any interior length
n>=2, with exact midpoint margin gamma_len and period b=n+2. Product
initialization has a separate finite hopping gap. Bare positive Bell
checks with a common exact history kernel give scalar second order.

[proved construction, numerical spectra] The richer NN check seeing a
partner and an auxiliary qubit gives a literal d0=64, b=4 cell with an
exact two-state code and strictly negative second-order ZZ coefficient.
This is not an unbounded-c construction. The full route-A energy
compiler and its work alphabet remain open, so no d0 for unbounded c
is reported. Prefer C19's new uniform Jacobi-filter target for that
compiler; do not redo its target synthesis.

All four Python scripts passed, all compiled, and the corrected clock
suite passed its final rerun (26.8 seconds). The final n=40 clock gap is
3.4710652271014505e-27, evaluated without cancellation. All generated
files are local to this folder. No network or external modifications.

Next work: a finite syndrome-energy menu using the richer boundary
interface, then its C8 sweep implementation with initialized logical
regions and weighted event timing. Preserve the fixed clock graph.

## Chronological record

2026-10-27. Read the requested STATUS, minimal_d0 README, C7 selector,
C10 finite-program construction and obstruction, programmable_cell README,
C9 spectral interface, and the current C8 progress and clock notes.
The project STATE.md was read. It is outside the writable scope and will
not be modified. C19/PROGRESS.md does not exist at the requested path.

## Scope and inherited results

[proved in C8, used as input] C8 now has a nine-letter initialized binary
clock and a 15-letter full-sweep extension. Both have one complete legal
path at each length, conserved local defects, and convex hopping ground
energies. These are not new C18 results and will not be duplicated as
research claims. C8 is working on the work-machine interpreter/readout.

C18 will concentrate on (i) a two-idle-port direct-sum interface and its
exact alphabet count, (ii) hopping initialization energy and the weighted
history version of the route-A spectral interface, (iii) a quantitative
selection margin and finite-coupling compatibility, and (iv) the precise
remaining route-A compiler obstruction. No universal machine or d0 near
45 is assumed. All new code and outputs will remain in this folder.

No computation has yet been run. The intended final verdict may be a
partial construction with a located obstruction, as allowed by the task.

## Step 1. Two-port wrapper and an exact-code obstruction

[proved, write-up in preparation] Replace the scalar C8 marker by two
qubit letters P_L,P_R, with legal boundary pattern P_R P_L and a clock
segment between P_L and P_R. The dimension is 15w+4 for the full-sweep
clock, or 9w+4 for the binary clock. The period is b=n+2, not n+1.
All gates act only on bulk work registers, including boundary moves.
Defects remain conserved and the fragment estimate remains K-4J per
forbidden bond. The chemical potential is split equally over the two
port letters. This is an interface and selector, not a route-A compiler.

[proved, write-up in preparation] Positive Bell-check conjugates plus
positive port-idle energies cannot give nontrivial route-A exchange if
they all annihilate the propagated code at every time. On each port
syndrome F_a G the own Bell check is identity, the other Bell check and
port-idle energies vanish, and propagation acts only on the clock.
The full source resolvent is scalar on the logical code. This obstruction
includes variable positive time weights and finite hopping strength.
It does not exclude non-frustration-free energy programs or richer
port-touching terms such as syndrome-conditioned workspace energies.

Next: quantify selector perturbation tolerance and weighted histories,
then independently check the wrapper rules and the obstruction numerically.

## Step 2. Quantitative selection and completed first checks

[proved] For the two-port period b=n+2 and midpoint chemical potential,
write a_n=2J[1-cos(pi/(C_n+1))]. The least wrong-cell excess is exactly
Gamma_len=(a_(n-1)-2a_n+a_(n+1))/2. Both adjacent lengths attain it.
The least competing energy-density difference is Gamma_len/(b+1).
For C_n=(2n+4)2^n+n-4, Gamma_len is asymptotic to 9a_n/8, rather
than a power of b. A small positive energy program preserves selection
if its target-cell energy shift is less than b Gamma_len/(b+1).

[proved] An endpoint pin of strength J on a C-vertex hopping path has
lowest energy -2J cos(pi/(C+1/2)). Its lift above the unpinned path is
asymptotic to J pi^2/C^3. Thus initialization is a smaller scale than
the C^-2 clock gap and cannot be dropped from the hierarchy.

[numerical] clock_checks.py passed 6,069 local rewrite contexts, all
15^n open words through n=4, and all 17^L port/clock ring words at
L=3,4. Illegal-component spectra were checked through n=3. The sole
legal components have 9,30,79,192 vertices. Selection was checked at
n=2,...,40 against n=1,...,100, and endpoint pins through C=5124.
An initial direct cosine subtraction lost the gap at n=40; the code
was corrected to a stable sine product before the final rerun.

[numerical] interface_checks.py passed. A two-port 288-dimensional
cell with variable positive Bell weights has scalar second-order bond
with nonscalar residual 1.10e-14. A three-time port-idle SWAP example
has positive mean energy 1/4 and limiting port compression norm 1/3,
so arbitrary positive Bell-frame averaging does not preserve route A's
zero first order. Richer syndrome-conditioned energies give nontrivial
port resolvents with the proved weighted-history error bound.

Next: write the complete proof and the restricted uniform target
specification. The missing repeated-bond energy compiler will remain
explicit rather than being counted as an implemented work alphabet.

## Step 3. Literal small route-A cell and live C19 handoff

[proved construction, numerical spectra] routeA_toy.py supplies a literal
15*4+4=64-dimensional site alphabet. At n=2, the outward sweep swaps
work qubits c_1,c_2 while leaving partner qubits s_1,s_2 and both ports
idle. The positive boundary checks are (1-Pi_Bell) tensor (1+Z_c/4),
which are NN because s and c are factors of the same bulk site.
A right-end I-sweep check pins c_2=0 initially. The complete ground
code is two dimensional and the selector is unchanged. Complete legal
cell spectra are obtained from four exact Bell blocks, total dimension
1920. The gap is .000941163372332468. The second-order ZZ coefficient
is -.011377003596467483 (J=eta=1). Exact two-cell eigenvalues in an
invariant 3600-dimensional sector have cubic residuals at lambda
.004,.002,.001. This is a small diagonal interaction, not a high-c model.

[proved specification, numerical JW check] target_spec.py generates a
free-fermion logical bond from b alone, using 3r+1 C9 operator slots.
Checks for r=1,2,3 give residuals below 4e-16. Large b values in these
checks are input integers only, not simulated physical chain lengths.
No runtime or space theorem for a compiled interpreter is inferred.

Live-source recheck: C19/PROGRESS.md now exists and was read. C19 has
already obtained a stronger uniform target using Jacobi filters and
explicit reflection dilation, with parity and exchange-generated
fields. Its RESULT should be the preferred compiler input. C18 will
not duplicate that synthesis. The modest C9 target specification here
is retained as an auxiliary checked formula, not presented as the
remaining target-side bottleneck. C19 also records the live C18
static-Bell obstruction. The remaining issue is the spatial energy
compiler and its actual work alphabet, not an abstract uniform target.

## Step 4. Complete derivation written

RESULT.md now contains the explicit two-port rule table and count,
all-sector proof, weighted history and initialization bounds, exact
midpoint margins, positive-program tolerance, Bell-check no-go proof,
and full d=64 small-cell specification. The toy ZZ term is proved
strictly negative by an entrywise-positive semigroup comparison, not
only by its measured value. The final program discussion uses C19's
new uniform target as the preferred input and lists the remaining
physical compiler obligations. No unbounded-c d0 is claimed.

Final audit in progress: rerun the corrected stable clock-gap evaluation,
compile all four Python files, and check the saved claims and file list.

## Final audit completed

Verified the endpoint-pin boundary condition, both period conventions,
clock vertex/edge counts, exact midpoint excess and density margin,
all illegal-sector comparisons, port-idle source identities, and the
scope of the Bell-check theorem. The latter explicitly requires an
active check for each port. General unitary work gates are described
by path-operator equivalence, not falsely as unbranched basis graphs.
The small cell's nonzero ZZ sign has an analytic proof and its 00
matrix element checks as -3/(5/4+5/4)=-1.2.

Final checks passed. No unresolved test failure remains. The only
corrected numerical issue was direct subtraction of nearly equal
cosines in the first clock-gap evaluation; the saved final data use
the stable sine product. Full two-cell ED was not performed; the
reported finite-lambda tests use exact invariant source sectors.
The unbounded-route compiler and d0 remain unresolved as stated in
RESULT.md, rather than being silently assumed.
