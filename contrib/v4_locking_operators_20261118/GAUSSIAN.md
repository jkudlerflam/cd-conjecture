# Luttinger matrix, drift and relevance

All phase criteria here are conditional on a stable two-component Gaussian
regime with the compactification in LEDGER.md. No values are fitted to V4.

## 1. General convention and symmetry reduction

Let Phi=(phi_c,phi_s,theta_c,theta_s)^T and

H_G=(1/2pi) integral dx (partial_x Phi)^T H (partial_x Phi).

H is real symmetric and positive definite. The canonical spatial algebra is
[Phi_i(x),Phi_j(y)]=(i pi/2) eta_ij sgn(y-x),
eta=[[0,I],[I,0]]. The overall chirality sign is conventional; **eta is
symmetric**, unlike a finite-dimensional oscillator symplectic matrix.
Define the real logarithmic covariance G by

< [Phi_i(x)-Phi_i(0)] [Phi_j(x)-Phi_j(0)] >_sym
=G_ij log(|x|/a0)+bounded terms.

Then x(exp(i l^T Phi))=l^T G l/4 and the two-point exponent is 2x.
For a single uncoupled liquid G=diag(K,1/K). This fixes all factors.

The exact unbroken symmetries allow

H = [[a_c, 0,   0,   lambda_1],
     [0,   a_s, lambda_2, 0],
     [0,   lambda_2, b_c, 0],
     [lambda_1, 0, 0, b_s]].

All four diagonal entries are positive and
lambda_1²<a_c b_s, lambda_2²<a_s b_c.
These are the exact positive-definiteness conditions for this Gaussian
matrix, not proof that a candidate ground state populates two valleys.
Without drift one can write a_nu=u_nu/K_nu, b_nu=u_nu K_nu. No equality of
u_c and u_s is required. With drift, these bare ratios do not themselves
specify the scaling dimensions. Linear partial phi_c and partial theta_s
are allowed zero-mode shifts; no derivative-only term pins a compact field.
The other cross terms are forbidden in a symmetry-preserving Gaussian state.

## 2. Correct covariance with drift

Set B=H^(1/2) eta H^(1/2). Since B is real symmetric and nonsingular,

\[
G=H^{-1/2}\,|H^{1/2}\eta H^{1/2}|\,H^{-1/2}.
\]

To derive it, diagonalize B=O Lambda O^T and put
Phi=H^(-1/2) O |Lambda|^(1/2) chi. The chi fields have diagonal signed
commutator metric and Hamiltonian velocities |Lambda|; their normalized
logarithmic covariance is the identity. Transforming back gives G above.
The eigenvalues of eta H are the signed velocities. This derivation neither
identifies the velocities nor assumes they are equal. G eta G=eta provides
an internal canonical check.

If density-current mixing vanishes but the density and phase blocks A,B
are general positive matrices, the usual Luttinger matrix K satisfies
K A K=B and G=diag(K,K^-1). This is consistent with the multichannel
normalization in [Kagalovsky, Lerner and Yurkevich, equations 8--10,
opened](https://publications.aston.ac.uk/id/eprint/30823/1/Local_impurity_in_multi_channel_Luttinger_liquid.pdf).
That paper studies impurity relevance; the threshold here is **bulk x=2**,
not an impurity threshold borrowed from it.

### Correction to C180, agreement with C191

[C180 optable.py](../../codex/work_C180/optable.py) instead uses iJ, with
J=[[0,I],[-I,0]], directly on the four spatial fields. That is incorrect
when drift is present. Choose all a_nu=b_nu=1 and lambda_1=lambda_2=0.6.
In the valley basis this is two free K=1 liquids with opposite drift;
right/left speeds are 1+0.6 and 1-0.6, while exponents remain unchanged.
The correct G is I and x_quartet=x_backscattering=2. C180 returns diagonal
entries 1.25, cross entries -0.75, and x_quartet=2.5.

This is not a dispute about the decoupled normalization: both formulas
agree at zero drift. It refutes the use of C180's drift-adjusted numbers
and the inference that separate bare K inequalities suffice with drift.
[C191 gauss191.py](../../codex/work_C191/gauss191.py) already uses the correct
signed chiral metric. The present implementation agrees with it for twenty
positive random Gaussian matrices. It also agrees with an independent
8-dimensional real-Fourier oscillator calculation, where the oscillator
symplectic matrix *is* appropriate after including cosine and sine modes.
Maximum discrepancy is below 3e-14. No prior files were changed.

## 3. Exact conditional inequalities on the generic operator lattice

Write Q as the submatrix of G for the commuting fields (theta_c,phi_s):
Q=[[a,h],[h,b]]. For the generic vertices of LEDGER.md,

\[
x_{r,s}=2(a r^2+b s^2-2h r s).
\]

The primitive quartet is relevant when **a<1**, relative backscattering
when **b<1**, and a mixed vertex when **a r²+b s²-2h rs<1**.
Equality is marginal: the amplitude and coupled RG flow decide its fate.
Strict irrelevance at a Gaussian point gives local stability only within
its basin; a finite large bare perturbation is not certified harmless.
Higher powers and all generated integer combinations are included by the
same formula.

For all generic permitted perturbations to be irrelevant, the necessary
Gaussian test is the shortest-vector condition

min_{(r,s) nonzero integers} (a r²+b s²-2h rs)>1.

Apply it to all nonzero amplitudes; a symmetry-allowed amplitude accidentally
zero at one point cannot protect an open phase without further evidence.
When h=0, the sufficient and necessary strict test reduces to a>1,b>1,
i.e. K_c<1 and K_s>1 in the uncoupled convention. Different velocities are
still permitted. A nonempty stable region exists in this abstract theory;
for example zero drift, K_c=1/2,K_s=2, away from other resonances. It is
not evidence that this microscopic bond realizes those parameters.

**Mixed terms matter even if the two primitive tests pass.** A valid
symmetry-allowed positive Gaussian example is
Q=[[1.1,0.9],[0.9,1.1]]. Its pure dimensions are 2.2 and 2.2, but
x_(1,1)=0.8. checks.py constructs the compatible full canonical G and H.
It is an algebraic counterexample, not an estimated V4 matrix.

The generic arguments contain only theta_c and phi_s, which commute.
One relevant independent cosine can lock one combination, leaving one
Gaussian mode (c=1) if no further instability occurs. Two independent
compatible relevant cosines can lock both (c=0). It is unjustified to
identify either remaining mode as pure charge or pure relative motion
when mixing is appreciable.

## 4. Resonances and non-null vertices

For any integer (n,m), rotate to charge/relative components and set
l=(-2m_c,-2m_s,n_c,n_s). Compute x=l^T G l/4, impose the actual momentum
and symmetry-orbit rules, and use x<2 for weak bulk relevance. This covers
commensurate charge locking, density-assisted pair creation and valley
transfer, including simultaneous resonances. One cannot decide resonance
using raw microscopic density or unresolved one-body maxima alone.

For straightforward simultaneous sine-Gordon pinning, the candidate
vectors must be mutually null under eta: l_i^T eta l_j=0, including i=j.
In (n,m) language the cross condition is
m_i dot n_j+n_i dot m_j=0. Non-null, nonzero-conformal-spin perturbations
can alter momenta or produce other instabilities; relevance alone is not
a proof that they gap one mode. The generic family is mutually null.

Relevant detuned vertices cease to add coherently beyond a length of order
1/|q|. Whether they reach strong coupling first depends on their amplitude,
velocities and running dimensions. This is why a generic incommensurate
selection rule is not a quantitative exclusion of a nearby locking tongue.

## 5. What the symmetries protect, and what they do not

Translation forbids generic two-particle pairing and direct valley transfer;
it does **not** forbid their momentum-zero products. A and Sigma relate the
valleys but permit the quartet, relative backscattering, mixed cosines and
drift. Thus they allow, but do not symmetry-protect, a c=2 phase. The quartet
is the principal additional number-breaking threat in passing from the
N-conserving projected model to the full spin chain. Relative and mixed
locking remain independent threats.

The relevance of an interaction depends on the complete Gaussian matrix,
not just the number of valleys. Published two-component zigzag field theory
also imposes a matrix-dependent backscattering criterion; see
[Hikihara et al., equations 60--61, opened](https://arxiv.org/pdf/1004.0622).
Those model-specific parameters are not V4 data.
