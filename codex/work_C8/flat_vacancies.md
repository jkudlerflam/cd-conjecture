# A uniform-vacancy wrapper with explicit coefficients

Status: proved finite construction and gap bounds, with complete-sector
finite checks. This is the preferred ninth-run wrapper for the explicit
clean-readout core. It does not prove the microscopic infrared limit.

## 1. The modification

The earlier wrapper used a positive vacancy Hamiltonian with diagonal
domain-wall cost U and hopping w<U/(k+1). Its two endpoint droplets
produce a small even/odd tunneling splitting. Replace it by the
vacancy hopping Laplacian, w=U, and set the occupied-length filter
strength to one. Keep the same arrow alphabet, endpoint rewards
mu=2U, and violation penalty K=6U.

The local bond is still arrow_wrapper equation (4), with these changed
parameters. On each allowed vacancy/particle pair the charge block is

    U [[1,-1],[-1,1]].

On unmatched vacancy/particle pairs it has the positive diagonal U.
Consequently the charge term is positive on the full physical space.
All arrow violation projectors are conserved as before. The local
dimension remains 6Q+3, or 96Q+48 after the selector and equalizers.

The earlier proof using a large filter energy above nonzero charge
ground energies does not apply at w=U. Sections 2-4 give its replacement.
This parameter change is not being justified by extrapolating the
earlier small-w estimate to its boundary.

For all explicit constants below set U=1. Restoring a general U in
the exclusion gap is immediate; the payload and filter need not be
rescaled.

## 2. Uniform charge states and an elementary gap bound

[proved] In a legal N-site sector with both endpoint tags and r occupied
sites, write p=N-r. The ordered occupied sequence is L B^(r-2) R.
Vacancies have uniquely determined arrow labels. For 0<p<N the charge
Hamiltonian is the simple-exclusion Laplacian

    H_ch=U sum_(j=1)^(N-1) (I-S_(j,j+1))

on the p-hole configurations, where S exchanges binary vacancy labels.
Its unique zero vector is uniform,

    chi_(N,p)=binomial(N,p)^(-1/2) sum_(|S|=p) |S>,

and its gap is at least U/N^2. The p=0 sector is one dimensional,
with zero charge Hamiltonian and no charge excitation.

Proof. Connectedness of nearest-neighbor vacancy motion proves the
kernel assertion. For the gap, introduce the complete-graph exchange
Laplacian H_comp=sum_(i<j)(I-S_ij), at U=1. In the spin-half
representation,

    H_comp=N(N+2)/4-S_tot^2.

The largest total spin is N/2. In a fixed vacancy-number sector its
fully symmetric vector is unique. The next allowed total spin is at
most N/2-1, so H_comp>=N(I-P_chi).

Compare a transposition i<j with the adjacent exchange word that moves
the letter at i to j and backtracks to exchange the endpoints. Its
length is 2(j-i)-1<=2N and it traverses each path edge at most twice.
For any function f of the hole configuration, telescope f along this
word and apply Cauchy-Schwarz. Each intermediate configuration remains
uniformly distributed, so the mean-square increment at edge e is the
Dirichlet form of I-S_e on f. Summing over i<j gives

    <f,H_comp f> <= N^3 <f,H_ch f>.

Indeed an edge has at most N^2/4 crossing pairs, each contributing
at most 2N times two edge traversals. Combining the two inequalities
proves the gap U/N^2. This is deliberately weaker than the known exact
exclusion gap and needs no gap theorem as an external input.

[proved] Every adjacent pair of occupied labels in occupied order is
physically adjacent with weight r/N in chi_(N,p). One proof uses gaps:
each hole configuration corresponds to a weak composition of p into
r+1 gaps, and these compositions have equal weights. Setting a
specified internal gap to zero leaves binomial(N-1,r-1) configurations,
out of binomial(N,r). Their ratio is r/N.

## 3. A positive pinning estimate

[proved] Let A>=delta(I-P), W>=0, ||W||<=M, and let P0<=P satisfy
W P0=0. Suppose

    P W P >= a(P-P0),    a>0.

Then ker(A+W)=ran P0 when A vanishes on P, and the gap above this
kernel is at least a delta/(a+M+delta). If P0=0 this is a lower bound
on the entire spectrum.

Proof. For psi orthogonal to P0 write x=Ppsi and y=(I-P)psi.
The triangle inequality for W^(1/2) gives the energy lower bound

    delta ||y||^2+(sqrt(a)||x||-sqrt(M)||y||)_+^2.

In the region with a positive bracket, the corresponding two-by-two
matrix has determinant a delta and trace a+M+delta. Its smaller
eigenvalue is at least the determinant divided by the trace. In the
other region, ||y||^2>=a/(a+M) for a unit vector, giving the stronger
lower bound a delta/(a+M). Positivity proves the kernel statement.

This lemma is also recorded in the collaborator IR_AUDIT. The proof
is included here to make the new wrapper argument self-contained.

## 4. Complete occupied-length selection at both physical lengths

[proved] Let the positive payload core at occupied length k>=8 have
complete kernel V and gap gamma_pay. Use the inherited filter of
length_filter.md at strength one, with phase 2pi/(k+2), left pin
|+(0)> and right pin |+((k-1)phase)>. Put

    gamma=min{gamma_pay,1/[k(2k-1)]},
    epsilon_* = kappa sin^2(pi/(k+2)) /
                [(1+sqrt(kappa))^2+sin^2(pi/(k+2))],
    kappa=1/[(k+2)(2k+3)].

The filter has a unique zero state at r=k and rejection energy at
least epsilon_* at every 2<=r<=k+2 with r!=k. These inherited filter
facts are used without rederivation.

Fix N=k+1 or k+2 and a legal sector with both endpoint tags. In ordered
coordinates let W be the sum of adjacency-conditioned payload/filter
bonds and the two endpoint penalties. Let P project onto uniform
chi_(N,N-r), tensored with the entire occupied internal space. Section
2 gives

    PWP = (r/N) sum_(a=1)^(r-1) h0_(a,a+1)+a_left+b_right
         >= (r/N) H_core+filter(r).

All terms on the right are positive. At r=k its kernel is exactly
V tensor |helix_k>, with gap at least gamma. Every adjacency term
and both boundary penalties annihilate that kernel, so it is also
annihilated by the unaveraged W. Section 3 therefore gives the complete
ground code

    chi_(N,N-k) tensor V tensor |helix_k>.

At r!=k the filter alone makes PWP>=(r/N)epsilon_* P. Applying the
same lemma with P0=0 gives a strictly positive rejection energy.
At r=N there is no charge space to eliminate, and the filter bound
epsilon_* applies directly. No condition on the payload spectrum
at these wrong occupied lengths is needed beyond positivity.

For the explicit clean-readout core the local norm bounds at unit
scales are ||h_pay||<=17, ||a_pay||<=7, ||b_pay||<=6. Hence

    ||W|| <= (r-1)18+8+7 =18r-3.

With delta=1/N^2, the lemma gives the following simple conservative
bounds, using gamma<=1 and N>=k:

    selected-sector gap >= gamma/(40N^3),
    wrong-occupation rejection >= epsilon_* /(10N^4).

For the first inequality use a=k gamma/N, denominator at most 20N,
and k>=N/2. For the second use a>=2epsilon_*/N and the same denominator.
The p=0 rejection obeys this weaker bound as well.

Restoring the rewards gives energy -4 at r=k. A legal sector with
one endpoint tag has energy at least -2, and one with none at least
zero. With v>=1 arrow violations, the endpoint count is at most
2(v+1), so its energy is at least 6v-4(v+1)=2v-4. Thus all these
other sectors have separation at least two from -4. This includes
missing oriented rewards at physical endpoints, which only increase
their energies.

Consequently the complete open-cell code at each N=k+1,k+2 is the
stated code, with gap at least

    min{2, gamma/[40(k+2)^3], epsilon_* /[10(k+2)^4]}.

The complete ground energy is -4 at both lengths. The charge sector
no longer adds an exponentially small tunneling scale. The payload's
binary-clock gap remains in gamma.

The lower bound H_wrap(l)>=-4 for every open length l, and the
periodic-wrapper nonnegativity bound, still follow from the same
arrow count and positivity. They did not require w<U/(k+1).
Therefore the existing selector proof applies with e1=e2=0,

    mu_selector=t mu_c+4,    e_base=t e_c.

All wrong physical lengths, selector violations and marker-free
rings are excluded by the inherited proof. The new selected-cell
gap above replaces the old Delta_1,Delta_2 in its final gap minimum.
The positive equalizer tracks can then be attached unchanged.

## 5. Exact port factors and marker transport

[proved] For N=k+p, p=1,2, uniform charge weights give

    r_p=k/(k+p),
    s_p=k(k-1)/[(k+p)(k+p-1)].

These follow by counting hole subsets that avoid one or two specified
sites. They multiply the original core's endpoint and two-site density
maps. Core reflection symmetry is not needed.

The typed-vacancy cut conversion still has scalar payload transport.
Every term in its positive overlap sum has the same amplitude, so

    beta_wrap=[(k+1)/sqrt((k+1) binomial(k+2,2))]^2
              =2/(k+2).

Thus every wrapper normalization is elementary. There is no finite
charge eigenvalue or eigenvector left to compute. The full periodic
bond is periodic_clean_bond.PeriodicBond(flat_vacancies=True).

At k>=8 both r_p and s_p exceed 1/2. With phi=1/b the compensators
therefore obey

    ||M_L||,||M_R||,||M_D|| <=4pi^2 b^2.

The earlier elementary bounds A_b^2>1/2 and
cos(1/(2b))^(2b)>1/2 give

    zeta=A_b^2 beta_wrap beta_helix >1/[2(k+2)],
    u<2(k+2) J0.

This polynomial transport bound replaces the old wrapper's constant
1/64 bound. The constant bound must not be reused for this variant.

For J0=t_f=1, a full-space perturbation norm bound is consequently

    B_flat <= 1+2m+4(k+2)
              +409600 pi^4 b^4 m^3+320 pi^2 b^2 m.

The proof is the same triangle bound as run9_audit Section 5, using
the smaller compensators and u<2(k+2). The unperturbed wrapper bond
itself has norm at most 45 at unit scales: 6 for the arrow penalty,
2 for the charge block, 18 for the internal positive core/filter
bond, and 10 and 9 for the two oriented boundary terms including
their rewards. This bound is independent of k and q.

## 6. Numerical checks and scope

[numerical] Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 flat_vacancy_checks.py

The exclusion-gap test covers every 1<=p<N through N=9, 36 cases.
The minimum measured gap times N^2 is eight, above the proved one.
The comparison inequality H_comp<=N^3 H_ch is also checked through
N=6 in every such sector, with numerical negative eigenvalues below
1.2e-13.

Filter-only wrapper controls at k=3,4 include every occupied-number
sector and the complete legal arrow space at both selected lengths.
The legal dimensions are 315,961,961,2903. All ground energies are
-4 within 4.3e-14; all gaps exceed the bounds. The measured complete
gaps are approximately .0485006,.0381404,.0271197,.0219921.
Illegal arrows are excluded analytically by Section 4, not by a
full-Hilbert-space diagonalization.

For the actual q=4 clean core, the new periodic bond has the same
214272048 on-site states and is checked on 3232 base and perturbation
columns separately. All tested entries conserve microscopic charge,
marker number and the overlapping selector-legality projectors.
The wrapper agrees with the independent dense fixture below 9e-16.
The actual 5124-vertex histories reproduce the two-cell fermion matrix
with zero error, and compensation errors are below 3.1e-14. The new
marker visibility is .1869524555657523.

These are finite-dimensional checks. The flat wrapper preserves the
same first-order m+1-mode target and exact symmetries. It does not
supply a connected massive-sector integration theorem, measured
microscopic central charge, full-class bound or fixed-d construction.
The original conjecture remains nonconclusive.
