# Small-d conjectures from onsite stacking and a mobile vacancy

Date: 2026-11-16. Author / AI assistance: OpenAI Codex, GPT-6 family.
No independent agent or human review was performed.

## Claim and status

**[Working conjecture: attainability, not an upper bound.]** There are
one-site translation-invariant, nearest-neighbour Hermitian chains with

\[
 (d,c_{\rm tot})=(2,2),\quad(3,3),\quad(4,4).
\]

The d=2 entry is the repository's existing numerical candidate, not a new
confirmation. The d=4 entry follows **exactly conditional on that seed** by
putting two independent qubits inside each four-state site. The d=3 entry
is a new application of the repository's vacancy construction to the same
seed; it additionally requires a stable two-boson seed phase and
finite-coupling charge-payload matching. The explicit bond is below.

The stronger recursive hypothesis is attainability of c=d for every fixed
finite d>=2. It requires the vacancy conjecture to survive each iteration.
It is less established than the first three targets. **Neither statement
asserts c_max(d)<=d**, finite c_max(d), or optimality of these examples.

**[Proved algebra / kinematics in this note.]** Exact onsite stacking,
static-level embedding with a uniform defect gap, and an entropy identity
for factorized ordered-particle states. The finite-ring first-order vacancy
projection, already derived in the repository, is checked independently
here for the qubit seed and generic complex bonds. **[Conditional continuum
deduction.]** The existing resonance and Gaussian bounds give a sufficient
stability test for lifting a general stable Gaussian seed. No new
thermodynamic central charge is measured or proved.

## 1. Scope and the useful conjecture table

Define F(d) as the supremum of the total ordinary nonchiral Virasoro
central charge in the following class:

* H_L=sum_j h_(j,j+1) on the full tensor product (C^d)^L, with one fixed,
  finite Hermitian bond repeated at **every physical site**;
* strictly nearest-neighbour interactions, with uniform onsite terms
  allowed as halves of adjacent bonds;
* compact unitary conformal infrared sectors with finite positive speeds,
  whose central charges are added even if their speeds differ;
* all couplings fixed before L tends to infinity. A singular-coupling
  projection or a zero-speed degeneracy is not an example.

For a symmetry-broken state choose a bulk sector rather than adding the
central charges of different vacua. No spin-rotation symmetry, integrability,
reflection, or time-reversal symmetry is required.

| d | Established benchmark c | Lower bound if only the qubit c=2 seed is confirmed | Target using vacancy lifting as well |
| --- | ---: | ---: | ---: |
| 2 | 1 | 2 | 2 |
| 3 | 2 | 2 | 3 |
| 4 | 3 | 4 | 4 |
| 5 | 4 | 4 | 5 |
| 6 | 5 | 5 | 6 |

The first column uses the SU(d)_1 permutation chain; these benchmarks are
reviewed in [Latorre and Sierra](https://arxiv.org/html/2403.17242v2).
The middle column uses only these benchmarks, exact stacking, and static
embedding. The last column contains conjectured realizations, not upper
bounds or new numerical measurements. For d=4 the middle column already
reaches the target, without a vacancy argument.

Negative searches in irreducible spin-1 or spin-3/2 SU(2)-invariant families
do not bound this larger class. In particular, the two-qubit d=4 construction
need not have irreducible spin-3/2 SU(2) symmetry.

## 2. Exact restrictions on any proposed maximum

### Onsite tensor products

**[Proved.]** F(ab)>=F(a)+F(b). Given bonds h_A and h_B, identify a physical
site as C^a tensor C^b and use

\[
 h_{ab}=(h_A)_{A_jA_{j+1}}+\gamma(h_B)_{B_jB_{j+1}},\qquad\gamma>0,
\]

with identities on the other factors. A permutation of tensor factors
makes the full ring Hamiltonian H_A tensor I + gamma I tensor H_B.
Every energy is E_A+gamma E_B, a product ground state has
S_A(region)+S_B(region), and the conformal charges add. Both factors occur
at **the same physical site**; the spatial unit cell is still one site.
This is the familiar product construction, also discussed in the original
c-d paper, applied to the revised qubit seed.

Consequently a single verified qubit c=2 example would immediately give
a d=4, c=4 example. No robustness under coupling between the copies is
needed for this existence claim: those couplings are exactly zero. Their
absence is allowed even if a generic perturbation would gap a mode.

### Static embedding and its gap

**[Proved.]** F(D)>=F(d) for D>=d. Let P project onto the d payload levels,
Q=I-P, M=||h||, and define on C^D tensor C^D

\[
 h_D=(P\otimes P)h(P\otimes P)
       +\frac K2(Q\otimes I+I\otimes Q),\qquad K>2M.
\]

In a sector with r sites in Q the payload bond acts only on bonds joining
two P sites. Replace the Q sites by fixed payload states. At most 2r
missing bonds have been restored, each of norm at most M. The variational
principle therefore bounds that sector's energy by

\[
 E_{\rm sector}\ge E_0(H_d)+(K-2M)r.
\]

The same argument applies to mixed reduced states if the Q labels are
entangled. The all-P sector is exactly H_d, and every other sector has a
uniform positive defect gap. This preserves the original conformal limit.

These statements concern the true supremum F, or a conjectured *sharp*
formula for it. They need not hold for an arbitrary loose numerical upper
estimate. If unbounded central charge occurs at some d_0, it occurs at
every D>=d_0. They provide no universal obstruction to d_0=2 or d_0=3.

## 3. An explicit three-state candidate

Use the qubit bond, including the field halves,

\[
 h_2=\tfrac12(XX-YY)+\tfrac1{160}(XY-YX)+10ZZ
       -\tfrac{19.976047}{2}(ZI+IZ).
\]

Its finite-chain c approximately 2 evidence and remaining convergence
limitations are recorded in the
[qubit note](../d2_updated_conjecture_20261116/README.md).
Let e_min,e_max be the extremal eigenvalues of this **two-site** bond,

\[
 a=(e_{\max}+e_{\min})/2=9.99428065774773\ldots,\qquad
 M=(e_{\max}-e_{\min})/2=20.00678065774773\ldots,
 \qquad\widehat h_2=(h_2-aI)/M.
\]

Then ||hhat_2||=1 and its fully occupied chain has the same ground-state
physics as the seed. Add one vacancy |v> to the two payload states
|up>,|down>. Let iota embed the qubit in those two states and set

\[
 \boxed{\quad h_3(\eta)=
 -\sum_{s=\uparrow,\downarrow}
 (|vs\rangle\langle sv|+|sv\rangle\langle vs|)
 +\eta(\iota\otimes\iota)\widehat h_2
                 (\iota^\dagger\otimes\iota^\dagger),\quad\eta>0.\quad}
\]

This is one 9 by 9 Hermitian bond, repeated unchanged at each link. There
is no vacancy chemical potential. Occupied number is exactly conserved;
the payload qubit number need not be. All terms, including the field
halves, vanish unless both sites are occupied. Replacing these field
halves by unconditional onsite fields would produce a different effective
payload and is not the construction being claimed.

The constant a is specified deliberately: a constant shift of the payload
bond becomes an occupied-neighbour density interaction after lifting. It
cannot be silently discarded in the qutrit model. `checks.json` contains
the normalized 4 by 4 payload and a machine-readable 9 by 9 example at
eta=10^-6. That eta illustrates the literal bond; it is **not** a certified
c=3 parameter point or a proved matching threshold.

**[Conjecture.]** If the seed is in the interior of a stable two-boson phase
with the symmetry and resonance margins described below, h_3(eta) has a
three-boson infrared theory for all sufficiently small fixed eta>0.
The order is: fix the seed, then fix an acceptable positive eta, then take
L to infinity. Setting eta=0 does not give the desired CFT.

## 4. What can be established about this lift

### Exact ordered-particle representation and first-order coefficient

For N occupied sites, order their positions and retain their N payload
labels. Vacancy hopping moves the positions without interchanging the
cyclic order of the labels. The payload bond acts on consecutive labels
only when their positions touch. Across the periodic cut a hop rotates
the payload word. Its momentum sets the charge boundary twist, so the
finite-ring Hilbert space is not an unrestricted product with the same
charge boundary condition in every momentum sector.

**[Existing algebra, independently verified here.]** In the zero payload
momentum sector, the free charge ground state has energy

\[
 E_{\rm ch}(L,N)=-\frac{2\sin(\pi N/L)}{\sin(\pi/L)}.
\]

Writing S_0 for the isometry to cyclically invariant payload words and U
for its product with that charge state gives exactly the following identity
for the embedded interaction with its coefficient eta omitted:

\[
 U^\dagger H_{\rm lift,int}U
   =g_{L,N}\,S_0^\dagger H_{\widehat h_2,N}S_0,
 \qquad
 g_{L,N}=\frac LN\left[
 (N/L)^2-\frac{\sin^2(\pi N/L)}{L^2\sin^2(\pi/L)}\right].
\]

Individual contact probabilities of the rooted ordered particles need
not agree; their cyclic projection is what yields this identity. In the
thermodynamic limit at density rho,

\[
 g(\rho)=\rho-\frac{\sin^2(\pi\rho)}{\pi^2\rho},\qquad
 g(1/2)=\frac12-\frac2{\pi^2}>0.
\]

Thus the *leading* ordered-payload Hamiltonian is precisely
eta g(rho) H_hhat2. This does not integrate out a uniformly gapped charge
sector: that sector is gapless, and finite-size perturbation theory is
not uniform in L. See the existing
[vacancy projection and momentum audit](../counterexample_campaign/infrared/README.md).
The large-U Hubbard factorization of
[Ogata and Shiba](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.2326)
is a useful precedent, not a theorem for this nonintegrable payload.

### A rigorous density window without tuning

**[Inherited proved bound; its proof was read.]** The
[general density estimate](../vacancy_continuation/DENSITY_WINDOW.md)
applies to this normalized bond with t=1 and zero chemical potential:

\[
 |N/L-1/2|^2\le\frac{3\pi\eta}{16}+\frac{\pi^2}{32L^2}.
\]

For 0<eta<=10^-3, every thermodynamic minimizing density lies strictly
between 19/40 and 21/40. The estimate requires no internal conformal
assumption and rules out an empty or fully occupied ground state. It
does not establish positive compressibility or exclude coexistence of
two densities inside that interval.

### An exact entropy identity for the factorized reference state

**[Proved kinematics.]** Consider any fixed-N wavefunction of the form

\[
 \Psi(x_1,\ldots,x_N;a_1,\ldots,a_N)
       =\phi(x_1,\ldots,x_N)\,\chi(a_1,\ldots,a_N),
 \qquad x_1<\cdots<x_N.
\]

Cut the physical ring into sites 1,...,ell and the remainder. If k particles
lie on the left, their payload labels are exactly the first k labels.
Within that block the charge and payload Schmidt decompositions tensor
multiply. Different k blocks are orthogonal on both sides. If p_k is the
charge probability of k particles on the left, then

\[
 \boxed{S_{\rm physical}(\ell)=S_{\rm charge}(\ell)
                       +\sum_k p_k S_\chi(k;N).}
\]

Indeed the reduced density matrix is, up to a local change of basis,
the direct sum of rho_charge^(k) tensor rho_payload^(k), with
Tr rho_charge^(k)=p_k and Tr rho_payload^(k)=1. Its entropy is
H(p)+sum_k p_k[S_charge,conditional(k)+S_payload(k)], yielding the formula.
The particle-number Shannon entropy is already in S_charge and must not
be counted again.

If phi is a free charge sea and chi is a payload CFT ground state, their
usual logarithmic entropies give the coefficient (1+c_payload)/3 on a
periodic chain. For an interval of fixed fraction, k=rho ell has only
O(sqrt(log L)) fluctuations: this follows from
Var N_A=sum_(i in A,j outside A)|C_ij|^2 and the free sine kernel's
inverse-distance decay. Hence replacing k/N by ell/L changes no leading
logarithmic coefficient. This is an entropy statement about the reference
state, not a proof that it is the finite-eta lattice ground state.

At eta=0 the entire payload space has zero energy and zero velocity.
One can choose highly entangled payload states in that degenerate space.
Even an exactly additive entropy coefficient there is insufficient to
establish an additional conformal mode. Positive fixed eta and infrared
matching are essential.

## 5. A precise stability conjecture, and its remaining input

Here is a useful restricted formulation rather than an unrestricted
claim that adding a level always adds a CFT.

**[Stable Gaussian vacancy-lifting conjecture.]** Suppose a fixed q-state
payload bond has a compact k-boson phase with finite positive velocities.
Require that, under its retained exact onsite symmetries and translation,
all allowed uniform non-Gaussian vertices have dimension at least
Delta_pay,0>2, and all potentially relevant oscillating vertices have
nonzero momentum margins. Uniform linear currents may shift densities;
quadratic currents are retained as Gaussian couplings. After adding a
species-blind mobile vacancy, sufficiently small fixed positive coupling
produces a compact (k+1)-boson phase near half occupation.

This is a sufficient, restrictive set of hypotheses. It does not cover
every isolated CFT point or every marginally irrelevant interaction.
The strict dimension and momentum margins have **not** been extracted
for the current qubit seed. Its entropy data alone do not establish them.

**[Conditional field-theory deduction.]** Specialize the existing vacancy
operator dictionary to one payload site per material cell. A payload
vertex of momentum p and dimension Delta_pay, dressed by integer harmonic
ell, has reference physical momentum and dimension

\[
 k_{\rm phys}=\rho(p+2\pi\ell)\pmod{2\pi},\qquad
 \Delta_{\rm ref}=\Delta_{\rm pay}
                     +[(p+2\pi\ell)/(2\pi)]^2.
\]

Exact occupied-number conservation forbids dual-charge vertices. A
uniform physical resonance k_phys=2 pi m has

\[
 \Delta_{\rm ref,res}=\Delta_{\rm pay}+m^2/\rho^2.
\]

For m!=0 the density window gives a charge contribution greater than
1600/441>3.62. Thus no new such locking vertex is relevant or marginal
at the reference point, even if its payload factor has dimension zero.
For m=0 the momentum selection reduces to the original uniform payload
problem. An irrational filling is not needed for this conclusion.

Let \(\Delta_0=\min(\Delta_{\rm pay,0},1600/441)>2\). For the full Gaussian energy
matrix, including all unequal velocities and density-current blocks, the
[existing matrix certificate](../digit_compilation/audit/GAUSSIAN_CERTIFICATE.md)
implies that the sufficient relative bound

\[
 \epsilon=\|G_{\rm ref}^{-1/2}(G-G_{\rm ref})
                       G_{\rm ref}^{-1/2}\|
 <\frac{\Delta_0^2-4}{\Delta_0^2+4}
\]

keeps every resonant non-Gaussian vertex irrelevant. Indeed its dressed
dimension is at least sqrt((1-epsilon)/(1+epsilon)) Delta_ref>2.
The low-order oscillating vertices must also retain their momentum
margins. A common fast charge velocity alone is not a stability bound
for the slow payload modes.

The reason small coupling is promising is specific: the smooth charge
density response is bounded in the continuum by K/(pi v), including at
nonzero frequency. The [existing response audit](../vacancy_continuation/local_audit/README.md)
therefore gives corrections of order eta^2/t to slow payload interactions,
compared with their order-eta scale. Gaussian charge-payload cross terms
of order eta have relative size of order sqrt(eta/t) at fixed seed; they
also become small. Constants depend on the seed's slowest velocity and
operator amplitudes. They need not be uniform across different seeds.

**Open step:** justify this local low-frequency matching for the actual
qubit payload at fixed eta and retain its phase and momentum margins.
The continuum response bound is not a microscopic operator-norm bound;
the density estimate alone is not a compressibility theorem. This is
where the qutrit claim is less certain than the exact d=4 stacking claim.

The existing physical d=3 vacancy-plus-XX calculations at L=8,12,16 give
entropy slopes near two, consistent with the ordinary c=1 to c=2 lift
([data and limitations](../vacancy_continuation/README.md)). They are
calibration evidence for the mechanism, not evidence measured for the
new c=2 to c=3 candidate.

## 6. Consequences for further small-d work

For d=4, use the exact two-copy construction first. For d=3, the sharp
next test is the vacancy bond above, together with a determination of
the qubit seed's two-mode scaling dimensions and allowed perturbations.
If these pass, d=5 can be reached by lifting a four-state payload, and
successive lifts give the tentative c=d attainability hierarchy.

Each iteration chooses its next coupling only after fixing the previous
payload, including its smallest positive velocity and stability margins.
For any fixed finite number of levels this permits a finite hierarchy;
no coupling is taken to zero with L. Tensor-product seeds can retain
separate exact onsite symmetries when needed to protect the factors.
Whether every required strict stability condition holds is part of the
recursive conjecture, not a consequence of tensor-product existence.

The stronger guesses F(2)=2, F(3)=3 and F(4)=4 require independent upper
bounds. This work supplies none. Composite modes, extra neutral sectors,
and selected long emergent cells remain possible escape mechanisms.
Finite parameter counting, local entropy log d, and the two-valley
description of one effective Hamiltonian do not exclude those mechanisms.
The separate conditional large-c constructions at larger local dimension
also make a universal extrapolation F(d)=d unjustified.

## 7. Reproduction and provenance

Run from the repository root with Python 3 and numpy:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 contrib/small_d_lifting_20261116/checks.py
```

The checks cover the complete two-copy spectra on L=3,4; static embedding
on L=3,4,5 in all 53 nonempty defect patterns; six vacancy projection
identities on L=6,7,8 with q=2 or q=3 payloads, including oriented complex
bonds; four exact entropy identities for random complex factorized states;
and Hermiticity and exact occupied-number conservation of the explicit
qutrit bond. No entropy fit or thermodynamic extrapolation is performed.

[Results](checks.json): the largest identity error is 2.31e-14, every
embedding satisfies the proved bound, and all assertions pass. The recorded
runtime is below one second with Python 3.12.14 / numpy 2.3.5 on this
machine; the JSON records the exact installed versions and elapsed time.
These finite calculations supplement the stated proofs and do not verify
the infrared hypotheses. [Code](checks.py) uses only numpy and the Python
standard library.

The vacancy projection, density bound, harmonic dictionary, and Gaussian
certificate are inherited and explicitly linked above. New work here is
their application to the small-d qubit seed, the explicit normalized
qutrit candidate, the entropy-identity derivation and checks, and the
resulting distinction between attainable c and a conjectured maximum.
The primary c-d paper, Ogata-Shiba abstract, and the cited repository
proofs were opened. The recent
[qutrit XX-chain paper](https://arxiv.org/html/2509.20431) was also opened;
its uniform c=2 result supplies context, not a c=3 construction.
