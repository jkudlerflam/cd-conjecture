# A sparse species graph gives a simpler protected high-channel target

## Claim

Replace the complete bipartite species graph in
[the imported torus proposal](../chat_followup/torus_encoding.md) by a
bipartite graph of girth at least six. The following finite-dimensional
and integer statements hold.

1. Every charge-neutral fermion monomial with at most four fields
   conserves each logical flavour separately. This conclusion needs no
   momentum assumption.
2. Point-line incidence graphs of projective planes over a prime field
   give a family with port dimension
   `q = 2(p²+p+1)` and `N = (p+1)(p²+p+1)` logical complex fermions.
   Consequently `N = Θ(q^(3/2))` is still superlinear in the port dimension.
3. In the Gaussian field-theory description, equally spaced Fermi
   momenta in `(π/4, π/3)` suffice to eliminate every potentially relevant
   non-Gaussian vertex near the free point. The smallest needed mismatch
   is `π/[6(N+1)]`. No four-flavour Sidon condition is needed.
4. The joint completely positive compression exists for any such graph.
   If its maximum degree is `D`, the elementary sufficient compression
   coefficient improves to `ε ≤ 1/(4qD)`, replacing the original
   `1/(4qN)` estimate.

These statements specify the code and its intended effective theory.
The companion [finite-core theorem](EXISTENTIAL_CORE.md) now constructs
its initialized, species-preserving history Hamiltonian with local
dimension `d_core ≤ Cq(1+log₂q)` and exact endpoint compressions.
Consequently `N/d_core → ∞`. Turning this ratio into `c/d → ∞` for the
full periodic chain also uses the mobile wrapper, finite-coupling
infrared transfer and treatment of ring lengths; those composition
steps are outside this graph-and-code note.

## Status

[proved] The graph construction, charge selection, compression bound,
and momentum inequalities below. [derived, Gaussian field theory] The
stability neighborhood in the space of scaling dimensions.
[proved in the companion file] An initialized finite core with exact
ports and `O(q log q)` local dimension. The final microscopic c–d
conclusion requires the separate route composition and infrared audit.

The graph gives a simpler protected target for that composition. The complete bipartite proposal gives more modes per port,
but has charge-neutral rectangle conversions that require additional
momentum detuning. The sparse graph removes those quartics exactly while
retaining enough modes to beat any fixed `Cq log q` alphabet overhead.

## How it was established and checked

- Author: Codex, GPT-6 family, 2026-10-27, in the project research session.
- Method: direct graph, charge-lattice, Choi-matrix, and bosonization
  calculations. No external literature result is required for the graph
  statements.
- Checks: `checks.py` enumerates the projective-plane incidence matrices,
  verifies their degrees, unique pairwise incidence and absence of four
  cycles, and enumerates every neutral integer vector of l1 norm at most
  four for the smallest graph. It also checks the Choi norm bound on a
  six-cycle, and the finite list of low-dimension momentum relations.
- Independent human review: none. The numerical checks validate examples;
  the general statements rely on the proofs below.

## 1. Graph-valued charges

Let `G=(L ⊔ R,E)` be a finite simple connected bipartite graph with `q`
vertices and `N` edges. The physical species register has orthonormal
basis `|s>` for `s∈L⊔R`. For an edge `a=(l,r)`, use

\[
F_a=|r\rangle\langle l|,\qquad w_a=e_r-e_l.
\]

The logical creation operator `f_a†` has charge `w_a`. The sign that is
`+1` on `L` and `-1` on `R` makes every `F_a` odd. Thus exact physical
parity remains available as a subgroup of the species torus. Total
`R`-species charge is logical particle number up to the code's constant
reference charge.

Let `W` be the oriented incidence matrix. A nonzero integer vector
`z∈ker W` is an integer circulation. Orient each edge according to the
sign of its nonzero coefficient. Conservation at every vertex implies
that the support contains a directed cycle: follow any outgoing edge;
conservation prevents a dead end, and finiteness forces a repeated
vertex. Removing repeated vertices gives a simple undirected cycle.
Therefore, if `g` is the girth,

\[
0\ne z\in\ker_{\mathbb Z}W
\quad\Longrightarrow\quad
\|z\|_1\ge g.
\tag{1}
\]

Alternating `+1,-1` around a shortest cycle attains equality. The kernel
has rank `N-q+1`, since the connected incidence matrix has rank `q-1`.

For a fermion monomial, let `z_a` be creations minus annihilations of
flavour `a`, summing both chiralities. Its torus charge is `Wz`, and
`||z||₁` is no greater than its number of fields. If `g≥6`, neutrality
and at most four fields imply `z=0`. Every flavour number is then
individually conserved by that monomial. This does not assert exact
microscopic `U(1)^N`: neutral six-field cycle conversions are allowed.

## 2. An explicit superlinear family

Fix a prime `p`. A point is a one-dimensional subspace of `F_p³`, and a
line is the kernel of a nonzero covector, also taken modulo nonzero
scalars. Connect point `[x]` to line `[ℓ]` exactly when `ℓ·x=0`.

There are `M=p²+p+1` points and equally many lines. Each line contains
`p+1` points, because its two-dimensional kernel has `(p²-1)/(p-1)`
one-dimensional subspaces. Any two distinct points span a unique
two-dimensional subspace and hence lie on exactly one common line.
Thus a four-cycle is impossible. Three noncollinear points and their
three joining lines give a six-cycle. In particular,

\[
q=2M,\qquad N=(p+1)M,\qquad D=p+1,\qquad g=6.
\tag{2}
\]

The graph is connected: two points share a line, and every line contains
a point. No appeal to the existence of projective planes at arbitrary
orders is made; the explicit prime-field construction suffices.

The companion finite-core theorem has `d_core≤Cq(1+log₂q)`, hence

\[
\frac N {d_{\rm core}}\ge
\frac{p+1}{2C[1+\log_2(2(p^2+p+1))]}
\longrightarrow\infty.
\tag{3}
\]

The charge register need not have power-of-two dimension: label queries
use `ceil(log₂q)` binary digits. If a particular compiler requires a
power-of-two register, padding to the next power costs less than a
factor two and does not alter (3). The finite-core format verifier removes unwanted padded or malformed
states from its low-energy space.

## 3. A finite and explicit momentum condition

Set

\[
k_a=\frac\pi4+\frac{a\pi}{12(N+1)},\qquad a=1,\ldots,N,
\qquad t_a=(2\cos k_a)^{-1}.
\tag{4}
\]

The free target with hopping `-t_a` and common chemical potential `+1`
has Fermi points `±k_a`. Its velocities are `2t_a sin k_a=tan k_a`, so
`1<v_a<sqrt(3)`. Distinct flavour backscattering carries momentum
`2(k_a-k_b)`, whose distance from `2πZ` is at least

\[
\delta=\frac{\pi}{6(N+1)}.
\tag{5}
\]

Single-flavour backscattering has momentum `2k_a∈(π/2,2π/3)`. The
same-sign two-flavour combination has `2(k_a+k_b)∈(π,4π/3)`, at distance
greater than `2π/3` from `2πZ`. Those and density/current terms exhaust
the neutral quadratic and nonzero quartic vertices by Section 1.
Coincident fermion fields either vanish or require derivatives; they
cannot lower a scaling dimension.

The arithmetic progression in (4) has many exact higher-order momentum
relations. This causes no problem for the neighborhood described next:
the associated non-Gaussian vertices have dimension strictly above two.

## 4. A uniform Gaussian stability neighborhood

Use the conventions of the imported torus note. For integer net chiral
creation vectors `u,v`, write `m=u-v`, `n=u+v`. Thus `m≡n (mod 2)`
coordinatewise, the charge is `Wn`, and the momentum is `m·k`. At the
free point,

\[
\Delta_0(m,n)=\tfrac14(\|m\|_2^2+\|n\|_2^2).
\tag{6}
\]

For a neutral vertex with `n≠0`, (1) gives `||n||₁≥6`. Coordinatewise
integrality and parity give

\[
m_a^2+n_a^2\ge 2|n_a|,
\qquad \Delta_0(m,n)\ge3.
\tag{7}
\]

For even `n_a`, the inequality follows from `n_a²≥2|n_a|`; for odd
`n_a`, `m_a²≥1` and `n_a²+1≥2|n_a|`.

If `n=0`, then `m=2u` and `Δ₀=||u||²`. The only nonzero integer
vectors with `Δ₀≤2` are `u=±e_a` and `u=±e_a±e_b`, `a≠b`; their
momenta were bounded in Section 3. Every other vertex has `Δ₀≥3`.

Consequently, a Gaussian scaling-dimension form `Δ` obeying

\[
\Delta(m,n)\ge\kappa\Delta_0(m,n),\qquad\kappa>2/3,
\tag{8}
\]

has no symmetry-allowed, nonoscillating non-Gaussian vertex of dimension
at most two, provided the finite mismatches in Section 3 persist. In
the usual parity-symmetric form

\[
\Delta=\tfrac14(m^T K m+n^T K^{-1}n),
\]

the sufficient condition is `2/3<λ_min(K)` and `λ_max(K)<3/2`.
More generally (8) states the needed condition on the full dimension
quadratic form, without assuming that mixed terms vanish.

This is a statement inside bosonization. It does not prove that a
particular microscopic system flows into this neighborhood. Currents
and their quadratic products remain; they are the Gaussian couplings
and velocity renormalizations, not missing non-Gaussian vertices. Under
shifts `|δk_a|<δ/8`, the smallest opposite-sign two-flavour mismatch
remains at least `δ/2`.

## 5. The compression is unchanged, with a better bound

Define `Γ_L(F_a)=f_a†` and `Γ_R(F_a)=f_a†P`, assign adjoints their
adjoint images, and all other matrix units zero. Define the unital map
on the two ports exactly as in the imported note,

\[
\Phi(A\otimes B)=\tau(A)\tau(B)I+
\epsilon[\tau(B)\Gamma_L(A)+\tau(A)\Gamma_R(B)],
\qquad\tau(A)=\operatorname{Tr}(A)/q.
\tag{9}
\]

In the Choi matrices of each `Γ`, every nonzero block has norm one and
the block support is the graph adjacency matrix `A_G`, up to a
permutation or transpose. For a block vector with components `x_s`,
the vector of output component norms is bounded componentwise by
`A_G (||x_s||)_s`. Therefore

\[
\|C_{\Gamma_L}\|,\|C_{\Gamma_R}\|
\le\|A_G\|\le D.
\tag{10}
\]

The depolarizing term has Choi matrix `I/q²`, and the perturbation in
(9) has norm at most `2εD/q`. For `0<ε≤1/(4qD)`,

\[
C_\Phi\ge I/(2q^2)>0.
\tag{11}
\]

Stinespring consequently gives the same simultaneous left/right port
compression as before. Its equivariant environment embedding uses one
occupation data site per edge and two dual-species strings. The whole
abstract encoding therefore fits in `N+2q` data sites, with reference
charge

\[
\sum_{l\in L}\deg(l)e_l+2\sum_s e_s.
\]

This is only representation and Hilbert-space accounting; it is not a
nearest-neighbour parent-Hamiltonian construction.

The simultaneous [sparse-isometry construction](../exact_route/README.md)
now gives a more explicit realization and larger compression coefficients
for these same graphs. Its covariance and sequential attenuation were
independently checked during this work. Equation (11) remains a separate
elementary existence bound; generic Choi square-root synthesis is no
longer necessary for this target.

## 6. A sharply limited obstruction to borrowing the Bell-check compiler

The usual maximally entangled vector
`|Ω>=q^(-1/2)Σ_s |s,s>` on two fundamental species registers is not
invariant, even up to phase, under independent simultaneous species
rotations for `q>1`: its summands have distinct weights `2e_s`.
Thus its rank-one projector does not preserve the torus.

More generally, a pure vector of definite species charge on two such
registers lies either in a one-dimensional `|s,s>` sector or in the
two-dimensional span of `|s,t>,|t,s>`. Its Schmidt rank is at most two.
For `q>2`, no full-Schmidt-rank torus-invariant rank-one Bell check exists
on those two registers. The same conclusion holds if neutral local
ancillas are fixed and are not counted as extra charged species.

This does not rule out the torus proposal. A dual species register has
the opposite weights and admits the invariant Bell vector; alternatively,
a dual can be represented by the `(q-1)`-site missing-species string
already used by the proposal. It does show that a standard one-bond
fundamental/fundamental Bell-check adaptation is not the missing proof.
A compiler must either implement those duals locally or use the
flagged final-time port preparation described in the imported note.

## 7. Conditional extension to the mobile length mode

Assume the exact marker conservation and the physical-momentum/gluing
dictionary of [C13, Sections 4–5](../../codex/work_C13/RESULT.md). This
assumption is essential: the physical marker field does not have the
compactification of an extra unconstrained fermion. At a possible
uniform resonance, C13's dictionary instead gives

\[
m\cdot k-2\pi h\nu=0\pmod {2\pi},\qquad
\Delta_{0,\mathrm{res}}=
\tfrac14(\|m\|^2+\|n\|^2)+h^2,
\quad h\in\mathbb Z.
\tag{12}
\]

Here `ν` is the long-cell filling, `Wn=0`, and marker conservation
forbids vertices that change length-particle number. For `n≠0`, (7)
already makes (12) at least three. For `n=0`, put `m=2u`; the only
nonzero integer tuples with `||u||²+h²≤2` are those of Section 4,
the pure length harmonics `u=0,h=±1`, and the mixed tuples
`u=±e_a,h=±1`.

Choose

\[
\nu=\frac{\sqrt2}{12}\in(1/10,1/8).
\tag{13}
\]

Pure length mismatch is `2πν>π/5`. Mixed momenta
`2k_a±2πν` have distance greater than `π/4` from `2πZ`.
Thus the smallest low-dimension cell-coordinate mismatch remains the
payload value (5); adding this length mode does not require a new
`N^-2` spacing condition.

Within C13's dictionary, the same lower bound by `κ>2/3` on the properly
normalized Gaussian dimension form excludes all remaining resonant
non-Gaussian vertices. At each fixed cell size, keep the finite
low-dimension mismatch list separated under dressing, as required by
C13. This gives a conditional stable `(N+1)`-boson target with total
Virasoro central charge `N+1`. It does not supply a microscopic mobile
cell, calculate its finite-coupling coefficients, or remove C13's
finite-ring charge/flux qualifications.

## 8. The finite-core construction and remaining composition

The sparse species target reduces the infrared task to a multicomponent
Luttinger neighborhood with a finite, explicitly spaced list of
backscattering mismatches. It is an alternative to the original target
when robustness matters more than maximum channel count.

[EXISTENTIAL_CORE.md](EXISTENTIAL_CORE.md) proves a finite initialized
Hamiltonian that prepares the equivariant sparse code, preserves all
`q` species counts exactly, exposes both ports at flagged times,
excludes unwanted inputs and controller states, and has local dimension
`Cq log q` with a constant independent of the program length. Its exact
branch rotations eliminate a gate-synthesis approximation parameter.
The two endpoints may have different positive readout time weights;
simultaneous flag activation is unnecessary.

[COMPILER_INTERFACE.md](COMPILER_INTERFACE.md) supplies the spectral
and readout arguments: charge-preserving decoration of a control graph,
a gap bound without commutation on malformed inputs, cyclic
compute/readout/uncompute, and clean endpoint observables. The
[oracle scheduler](../exact_route/ORACLE_SCHEDULER.md) and concrete sweep
interface meet its constructive controller requirements without a
second charged species register.

Applying that finite core to the full conjecture uses the C8 mobile
wrapper and the local finite-coupling and infrared arguments. The
Gaussian claims in this note retain their stated field-theory
assumptions; the finite-core theorem by itself does not assert a
thermodynamic central charge.

## Reproducing

Run `python3 checks.py` in this folder. It uses NumPy and Python's
standard library, takes seconds, and writes `checks_results.json`.
It checks primes `2,3,5,7`, the six-cycle compression at six logical
modes, and the momentum lists for all resulting channel counts.
No many-body thermodynamic central charge is estimated here.

## Relation to the rest of the repository

This refines [the torus construction](../chat_followup/torus_encoding.md)
and the target choice in [infrared control](../ir_control/README.md).
This note supplies the protected graph target and links its exact
finite-core realization. It also isolates a symmetry issue in
transplanting the [programmable-cell](../programmable_cell/README.md)
Bell-check menu. The full c–d conclusion belongs to the combined
microscopic and infrared argument, rather than these graph statements
alone.
