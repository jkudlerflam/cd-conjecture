# Frozen selectors and incompatible ring lengths

2026-10-27, resumed C8. The frozen-marker theorem is proved under the explicit energy separation hypotheses below. The mobile-defect statement is proved only for the stated finite-band model. Neither is a theorem excluding every possible spontaneously blocked critical chain.

## 1. Frozen marker patterns force unbounded ground-state multiplicity

Let H_L be one-site translation invariant. Suppose every site has a binary marker projector M_i and

    [H_L,M_i]=0    for every i.                                (1)

Thus marker patterns label orthogonal invariant sectors. The marker can coexist with arbitrary quantum spectator degrees of freedom. Assume there are constants e in R, eta>0, C_0,C_1, independent of L, with the following properties for all sufficiently large L.

(a) Each constant marker pattern, all zeros or all ones, has energy at least (e+eta)L-C_0 in every state.

(b) Some state with a nonconstant marker pattern has energy at most eL+C_1.

Here e need not be positive. The notation eta>0 is the energy-density separation, not a conformal gap. These hypotheses formalize a stable marker crystal and a bounded cost for fitting it on an incompatible ring.

Theorem. For every sufficiently large prime ring length L, the ground-state multiplicity is at least L.

Proof. For L>(C_0+C_1)/eta, (b) lies below both constant-pattern sectors in (a). A ground state can be chosen in a joint eigenspace of all M_i because of (1). Its marker pattern is nonconstant. On a prime ring any nonconstant binary pattern has translation orbit of length L: its orbit length divides L and cannot be one. Translating the chosen ground state through this orbit produces L states of the same energy in mutually orthogonal marker sectors. They are linearly independent. QED.

This ground-state multiplicity is incompatible with convergence of the complete low-energy spectrum to a fixed finite collection of compact unitary CFT sectors with finite vacuum multiplicity. It does not become a larger central charge. Restricting by hand to one marker sector or to compatible lengths changes the spectral question.

The theorem is not a statement about CFTs with an explicitly allowed infinite vacuum multiplicity or a divergent torus trace. Such a convention would need to be specified separately and would also change the modular counting assumptions in C7.

### Application to the C7 selector

For a fixed selected period b, let e=e_b<0 be the optimal cell energy density from C7. A large ring of any length can be tiled by cells of size b except for one cell whose length lies between b and 2b-1. The excess energy over e_b L is bounded by

    C_b = max_(0<=r<b) [E_(b+r)-e_b(b+r)],                     (2)

which is finite for this fixed b. This proves (b) for the bare selector. The no-marker legal sectors have energy zero. C7's forbidden-bond estimate excludes lower-energy illegal sectors; for the bare selector the all-marker pattern is also uniformly costly. Thus (a) holds with a positive eta for the two constant patterns, and the theorem applies.

An encoded spectator construction inherits the same conclusion if it retains exact marker conservation and the two energy-density properties after adding its couplings. It is not necessary that the selected nonconstant ground pattern consist of only one defective cell. Primality alone gives the full orbit as soon as the constant patterns lose.

Consequently a proof that a frozen-marker encoding has the desired CFT on lengths divisible by b cannot meet the stronger all-ring requirement by appealing to a bounded defect energy. The defect position supplies an unbounded exact degeneracy. This is independent of whether the first-order code or second-order spectral gadget has been realized successfully.

## 2. Making an isolated defect mobile does not automatically fix this

Consider the following narrower effective model. A forced defect has a position x on a ring and a finite internal space of dimension s independent of L. Its translation-invariant Hamiltonian has finite-range hopping matrices,

    H_def=sum_x sum_(r=-R)^R |x+r><x| tensor A_r,
    A_(-r)=A_r^dagger.                                        (3)

Exponential hopping decay with an analytic Bloch matrix is equally sufficient. The Bloch matrix is

    h(k)=sum_r exp(ikr) A_r,

a Hermitian analytic matrix for real k. Let E_* be the minimum of its lowest band over the Brillouin zone.

Proposition. Such a defect has finite-size excitation energies above its ring ground state bounded by O(L^(-2)), or an increasing exact degeneracy. A nonzero linear cone cannot be the bottom of its finite-band spectrum.

Proof. Choose a normalized eigenvector u of h(k_*) at a minimizing momentum, with eigenvalue E_*. The scalar analytic function f(k)=u^dagger h(k)u obeys f(k)>=E_* for every real k and f(k_*)=E_*. Thus f'(k_*)=0 and f(k)<=E_*+C(k-k_*)^2 near k_*. For any fixed integer n, choose n distinct allowed ring momenta within O(n/L) of k_*. Their plane waves with internal vector u are orthogonal trial states, all with energy at most E_*+C_n/L^2. The Hamiltonian is diagonal in momentum, so the same bound holds throughout their span. The min-max principle gives at least n eigenvalues below that upper bound. The ring ground energy is at least E_*, proving the assertion. Since n is arbitrary, a fixed finite vacuum multiplicity cannot absorb all these levels after conformal rescaling by L. No assumption about analytic eigenvalue branches is needed. QED.

If this defect is decoupled from a critical spectator theory, the product spectrum contains extra quadratic or flatter levels alongside the conformal towers and is not purely a fixed finite collection of z=1 CFT spectra. This proves why an ordinary massive defect hopping term does not repair the frozen-selector construction.

The proposition does not apply without further work to a defect coupled to a gapless many-body bath. Its internal low-energy space may not have fixed finite dimension, and its effective dispersion need not follow from the analytic finite-band model (3). No universal theorem excluding a successful interacting repair is claimed. A complete counterexample must either control this coupling or avoid the forced-defect mechanism.
