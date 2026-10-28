# Research continuation: a protected first-order counterexample argument

## Decision and status

**The strongest route found here is a species-conserving, first-order
encoded chain using a sparse incidence graph. The continuation develops
it into a candidate refutation at the project's theoretical-physics
standard.** Read [the complete counterexample argument](PROPOSED_COUNTEREXAMPLE.md)
first. Its predicted total Virasoro charge grows as \(p^3\), while its
microscopic local dimension grows at most as \(p^2\log p\).

The new finite constructions include exact simultaneous ports, an exact
charge-preserving preparation circuit, a reversible scheduler with only
\(O(q\log q)\) local states, and an initialized positive history core.
The existing mobile wrapper applies to this core. A graph of girth six
forbids every flavour-changing quartic exactly, and explicit momenta
avoid the remaining dangerous resonances.

The conformal conclusion uses standard local elimination of massive
machinery and Gaussian infrared stability, audited in [IR_AUDIT.md](IR_AUDIT.md).
It is not a new rigorous thermodynamic theorem or a measured numerical
violation. The finite compiler and the complete composition need
independent end-to-end review. The stronger claim of unbounded central
charge at one fixed local dimension is not established by this family;
route A remains the separate approach to that goal.

The [microscopic bond specification](MICROSCOPIC_HAMILTONIAN.md) assembles
the local Hilbert space, selector, wrapper, compensated hopping and marker
motion in one place. It also fixes the first-order chemical potential
explicitly. This remains a parameterized construction: no complete
autonomous controller table, selected violating prime, numerical local
dimension, or verified finite coupling is supplied.

[Executable controller work](../controller_instantiation/README.md) now
provides the explicit encoding circuit, finite dispatcher and length-counter
tables, an initialization-format repair, and a fully saved nonviolating
microscopic bond fixture. Autonomous program generation, verified program
initialization and endpoint readout still need to be composed into the
actual high-channel core. The fixture is not a violating numerical example.

## Provenance

This work continues the repository at commit `5360803`, after pulling the
collaborators' latest C14 result. It was requested by Pace and carried
out by Codex (GPT-6 family), with three parallel research workers, on
2026-10-27 local time / 2026-10-28 UTC. The finite-dimensional arguments
were cross-reviewed by separate AI workers and rerun locally. No human
review and no thermodynamic numerical verification are claimed.

## 1. Why prioritize this route

| Route | Existing strength | Decisive unresolved requirement | Priority |
| --- | --- | --- | --- |
| Sparse species-graph, first order | Exact ports, initialized \(O(q\log q)\) core, protected free reference and mobile composition | Independent audit of the compiler/composition and the stated local infrared inputs | Candidate physics-level refutation |
| Route A, second order | Arbitrarily many logical modes through qubit ports; exact parity now repaired | Fixed-alphabet realization of the boundary spectral data and mobile completion | First for unbounded \(c\) at fixed \(d\) |
| Generic proof of c-d | Restricted classes and scale-dependent bounds | A new principle bounding conformal onset or neutral sectors by \(d\) alone | Lower until such a principle is found |
| More small-d scans | Useful exclusions and calibrations | A credible new mechanism, beyond previously failed finite-size signals | Lower without a specific target |

C14 refutes a proposed velocity-weighted replacement, not c-d itself.
It removes one possible general proof strategy without resolving the
original conjecture. Published finite-size universal-Hamiltonian results
also cannot simply be called a counterexample: preserving a finite
gadget does not by itself establish its thermodynamic CFT.

## 2. New exact ingredients

### Sparse simultaneous port encoding

[proved] For any simple bipartite graph with \(u\) upper vertices,
\(v\) lower vertices and \(N\) edges, put \(q=u+v\). The explicit
isometry in [exact_route](../exact_route/README.md) obeys

\[
W^\dagger F_e^L W=\epsilon f_e^\dagger,\qquad
W^\dagger F_e^R W=\epsilon\left(1-\frac1{2N}\right)f_e^\dagger P,
\qquad \epsilon=\frac1{2\sqrt{uN}}.
\tag{1}
\]

It preserves all species charges and fits as a representation in
\(N+4q\) data sites. Each occupation basis input has only \(u+N\)
branches per encoding stage. This is a concrete reversible recipe,
not an appeal to an arbitrary exponentially large matrix square root.
The physical hopping bond compresses to ordinary logical fermion
hopping with factor

\[
\eta=\frac{1-1/(2N)}{4uN}>0.
\tag{2}
\]

Dense checks construct the entire joint isometry, including signs and
charges; errors are below \(7\times10^{-16}\) in the recorded small
cases. This checks algebra, not a central charge.

### Exact local preparation and initialized history core

[proved by construction] [CIRCUIT](../exact_route/CIRCUIT.md) implements
the sparse isometry using neutral reversible computation and local
species-preserving routing. [EXACT_AMPLITUDES](../exact_route/EXACT_AMPLITUDES.md)
prepares its amplitudes exactly using an indexed algebraic rotation
opcode. The menu has fixed cardinality and its local alphabet remains
\(O(q\log q)\).

[proved using standard reversible computation] The
[oracle scheduler](../exact_route/ORACLE_SCHEDULER.md) keeps its classical
trajectory independent of quantum inputs. The explicit
[moving-head interface](../exact_route/SWEEP_INTERFACE.md) overlays only
one charged record, so the interpreter costs a constant factor rather
than a second factor of q. The already existing fifteen-letter clock
executes this ordered sweep; its full graph was independently checked
on all 50,625 words at length four.

[proved by composition] The [initialized-core theorem](../torus_bottleneck/EXISTENTIAL_CORE.md)
adds trusted scratch initialization, a reversible format verifier,
all-time FAIL penalties and cyclic, flagged two-endpoint readout.
Its complete code has dimension \(2^N\), exact zero energy and a positive
finite-cell gap. The on-site dimension is at most \(Cq(1+\log_2q)\).
No numerical value of the compiler constant C or fully expanded
universal-machine transition table is claimed.

### A graph that removes the dangerous quartics

[proved] In [torus_bottleneck](../torus_bottleneck/README.md), take the
point-line incidence graph of the projective plane over \(\mathbb F_p\),
where \(p\) is prime. With \(s=p^2+p+1\),

\[
u=v=s,\quad q=2s,\quad N=(p+1)s.
\tag{3}
\]

The graph has girth six. Every nonzero integer circulation has l1 norm
at least six. Hence any species-neutral fermion monomial with four or
fewer fields conserves each flavour. Flavour-changing vertices have
free scaling dimension at least three. Their irrelevance persists
when the full scaling-dimension form is bounded below by
\(\kappa\) times the free form, with \(\kappa>2/3\).

The remaining low-dimension backscattering mismatches can be separated by

\[
k_a=\frac\pi4+\frac{a\pi}{12(N+1)},\qquad
\delta=\frac\pi{6(N+1)}.
\tag{4}
\]

The target has one gapless complex fermion per edge. Under the standard
Gaussian infrared identification its Virasoro central charge is \(N\).
The graph and momentum statements are exact; transfer from a microscopic
Hamiltonian to that infrared theory remains a separate requirement.

### Route A compatibility repairs

[proved] [graded_route_a](../graded_route_a/README.md) identifies why
C9's diagonal syndrome prescription does not preserve logical parity
for odd hopping factors, and repairs it with off-diagonal even-to-odd
syndrome weights and a parity-equivariant dilation. For the free target,
the number of operator channels is linear in the number of flavours,
and the required whole-cell memory is \(n+O(\log n)\) qubits. Hopping-only
weights have a closed-form square root and can be aligned by a linear
number of Majorana rotations.

[proved] [separated_port_cells](../separated_port_cells/README.md) then
realizes independent left and right spectral weights as
\(H_L+H_R\), with no term touching both ports. The gadgets share the
logical register and need not commute; their exact product code and
gap follow from commuting boundary-code projectors. Full two-cell
resolvent checks agree below \(4\times10^{-15}\).

These repairs strengthen route A. They do not replace its missing
spatial compiler.

## 3. Symmetry inheritance through the existing mobile wrapper

[proved, conditional on the stated core] The existing C8 mobile wrapper
can preserve a species torus, not just parity. Suppose its input is a
\(d_{\rm core}\)-state uniform positive core

\[
H_{\rm core}=\sum_j h_{0,j,j+1}+a_1+b_k,
\]

with the complete desired zero-energy ground code and a positive gap at
the selected length. In particular, require \(h_0,a,b\ge0\) and
\(H_{\rm core}V=0\), so every term annihilates the code. Require each
term to commute with the appropriate sum of
on-site charge operators, and each port operator \(F_e\) to have charge
\(w_e\). Give every added selector, arrow, vacancy and helix label
zero species charge. On occupied wrapper letters retain the input
core's representation.

The inheritance proof is termwise:

1. Diagonal legality, density and endpoint-tag projectors are neutral.
2. Vacancy hopping exchanges a vacancy with an unchanged occupied
   internal label, so it preserves all total species charges.
3. Lifted \(h_0,a,b\) preserve charge by the input hypothesis. Conditioning
   them on neutral tags or selector projectors does not change that.
4. The marker conversion is \(|H_2\rangle\langle H_0|\) on a vacancy,
   tensored with neutral selector moves. Both vacancy states are neutral.
5. The equalizer's helix qubits and endpoint compensators are neutral.
   Multiplying \(F_e\) by them retains its charge \(w_e\).
6. A link \(F_e^R\otimes(F_e^L)^\dagger+\mathrm{h.c.}\) has zero
   total charge. The marker number is conserved independently.

Thus the full construction in
[selector_wrapper](../../codex/work_C8/selector_wrapper.md) and
[helix_equalizer](../../codex/work_C8/helix_equalizer.md) inherits the
species symmetry exactly. The local dimension is

\[
d_{\rm mobile}=96d_{\rm core}+48.
\tag{5}
\]

This constant overhead does not spoil the graph's superlinear mode
count. Applying it to the new initialized core gives the conditional
infrared comparison

\[
\frac{N+1}{d_{\rm mobile}}
\ge\frac{(p+1)(p^2+p+1)+1}
{192C(p^2+p+1)[1+\log_2(2(p^2+p+1))]+48}
\longrightarrow\infty.
\tag{6}
\]

Here \(N+1\) is the *conditional* payload-plus-length central charge.
The finite core supplies the positive uniform-core hypothesis; a
whole-cell isometry alone would not. Exact algebraic preparation avoids
code-splitting approximation errors. Equation (6) becomes a violation
under the infrared inputs stated in the complete argument.

The common density field also has a local implementation under this
hypothesis. If total upper-species charge restricts on one core to
\(Q_{\rm logical}+C_{\rm core}\), the full wrapped charge restricts to
\(\sum_j Q_j+C_{\rm core}M\), where \(M\) is the number of cells.
Both allowed cell lengths contain the same number of occupied core
letters, and the extra vacancies are neutral. Therefore adding a
chemical potential for that physical charge and subtracting
\(C_{\rm core}\) times the same coefficient from the marker chemical
potential produces the desired common logical density field. For the
explicit sparse encoding with neutral computational work,
\(C_{\rm core}=4u\). Loaded circuit heads transport existing charge;
they do not introduce an extra offset.

## 4. Infrared transfer and the highest-value remaining review

The [infrared audit](IR_AUDIT.md) specializes C13 to exact species
symmetry. Pairing is forbidden. With length density in
\([5/12,7/12]\), all dangerous mixed length/payload backscatterings have
momentum mismatch at least \(\pi/6\), independent of N. The smaller
payload mismatch remains \(\pi/[6(N+1)]\). Every other non-Gaussian
resonance starts at dimension at least three. There is therefore a
strict Gaussian stability neighborhood at each fixed p and cell length.

The [scale hierarchy](SCALE_HIERARCHY.md) includes all port probabilities,
wrapper norms and small core gaps before choosing a fixed nonzero
coupling. Connected virtual corrections start at second order, whereas
both retained energy scales are first order. With standard local
massive-sector integration, the relative corrections tend to zero and
enter that neighborhood regardless of their signs. No global spectral
band isolation or ring-length-dependent coupling is used.

The strongest next independent check is the complete compiler and
positive-core composition, followed by the local excursion argument.
At a strict mathematical standard an all-order uniform low-energy
reduction and the interacting CFT identification remain unproved. At
the project's stated physics standard they are explicit effective-field-
theory inputs to a candidate refutation, not another unidentified
quartic coupling that must be tuned. An expanded alphabet would make
the candidate more readily auditable and give a concrete first
violating dimension.

## 5. Reproduction and review

The adjacent contributions contain their complete NumPy-only scripts
and recorded outputs. Run them from the repository root:

```sh
python3 contrib/exact_route/verify.py
python3 contrib/torus_bottleneck/checks.py
python3 contrib/graded_route_a/checks.py
python3 contrib/separated_port_cells/checks.py
python3 contrib/graded_route_a/sweeping_interface_checks.py
python3 contrib/exact_route/sweep_checks.py
python3 contrib/exact_route/amplitude_checks.py
python3 contrib/exact_route/helper_checks.py
python3 contrib/graded_route_a/clock_schedule_checks.py
python3 contrib/research_continuation/ir_audit_checks.py
```

The graph results and sparse isometry were independently rederived by
two workers and rerun by another. The separated-port proof and numerical
implementation were reviewed independently of their author. The graded
repair was checked both directly and through the separate-port variant.
These are reproducible finite-dimensional results; no new entropy,
Casimir, or finite-size spectrum estimate of a microscopic CFT is being
reported.

The exact finite results and the physics-level conformal inference are
separate claims. The latter should be challenged through the explicit
assumptions and composition in [PROPOSED_COUNTEREXAMPLE](PROPOSED_COUNTEREXAMPLE.md),
not inferred from the small-matrix residuals alone.
