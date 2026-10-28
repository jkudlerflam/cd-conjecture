# An explicit Gaussian stability certificate for the species-graph route

**Status: proved matrix bound; conditional infrared application.** The
species-graph candidate's Gaussian-neighborhood requirement can be replaced
by the explicit sufficient condition

    ||G_ref^(-1/2) (G-G_ref) G_ref^(-1/2)|| < 5/13.

This condition includes arbitrary density-current blocks and noncommuting
velocity matrices. It guarantees that every vertex of reference dimension
at least three has dimension strictly greater than two. Translation
nonresonance and the correct local operator lattice are separate inputs.
The condition is a sufficient bound, not a necessary stability boundary.

Author: OpenAI Codex, GPT-6 family, 2026-10-28. This is an independent
derivation from the displayed Gaussian Hamiltonian in C13 and
`research_continuation/IR_AUDIT.md`. No human review or microscopic
thermodynamic calculation is claimed. All code in this folder is independent
of the candidate's implementation.

## 1. Covariance with a general density-current block

[assumed, Gaussian field theory] Let there be n compact bosons, use
X=(phi_1,...,phi_n,theta_1,...,theta_n), and write

    H = (1/(2 pi)) integral dx (partial_x X)^T G (partial_x X),
    G = [[A,C],[C^T,B]] > 0,
    J = [[0,I],[I,0]],
    [X_i(x), partial_y X_j(y)] = i pi J_ij delta(x-y).

An overall Fourier/sign convention for the commutator does not affect the
equal-time covariance. A vertex exp(i ell^T X) has dimension

    Delta_G(ell) = ell^T M(G) ell / 4.                 (1)

[proved within this Gaussian theory] The metric is

    M(G) = G^(-1/2) |G^(1/2) J G^(1/2)| G^(-1/2),    (2)

where |Q|=(Q^2)^(1/2) for a real symmetric Q. To obtain (2), define
Y=G^(1/2)X. The Hamiltonian has unit matrix in Y coordinates and its
commutation matrix is V=G^(1/2)JG^(1/2). Orthogonally diagonalize
V=O diag(v_i) O^T. Its inertia is (n,n), by congruence with J, so no
v_i vanishes. The fields

    chi_i = (O^T Y)_i / sqrt(|v_i|)

are independently normalized chiral fields. Their equal-time logarithmic
covariance is -delta_ij log|x|/2. Thus the logarithmic covariance of X is
-M(G) log|x|/2, which gives (1). The sign of v_i chooses chirality; its
magnitude sets velocity and drops out of the normalized chiral covariance.
This derivation does not require C=0 or simultaneous diagonalization of
A, B and C. It also explains why using only the eigenvalues of A or B
would be insufficient.

As a normalization check, for G=diag(v K^(-1),v K), (2) gives
M=diag(K,K^(-1)). The standard backscattering vertex exp(2i phi) has
dimension K, while a free complex fermion exp(i(phi +/- theta)) has
dimension 1/2.

For one mode with G=[[A,C],[C,B]] and |C|<sqrt(AB),

    M=diag(sqrt(B/A),sqrt(A/B)),
    v_right=C+sqrt(AB), v_left=C-sqrt(AB).              (3)

Thus a stable common drift C changes velocities but not dimensions. In
multiple modes a general C need not behave this way; (2) handles it.

## 2. A self-contained matrix proof of the bound

[proved] Define the positive matrix geometric mean by

    P # Q = P^(1/2) (P^(-1/2) Q P^(-1/2))^(1/2) P^(1/2).

It is the unique positive solution Z of Z P^(-1) Z=Q. Inverting this
identity gives Z Q^(-1) Z=P, so P#Q=Q#P. It is monotone in Q because
the positive matrix square root is operator monotone. For completeness,
the latter follows from the convergent resolvent integral

    sqrt(X) = (1/pi) integral_0^infinity
                t^(-1/2) X(X+tI)^(-1) dt.

If 0<X<=Y then (X+tI)^(-1)>=(Y+tI)^(-1), and therefore
X(X+tI)^(-1)=I-t(X+tI)^(-1)<=Y(Y+tI)^(-1). Integrating proves the
square-root statement. Symmetry of # now proves monotonicity in its
first argument as well. Direct substitution also gives

    (alpha P) # (beta Q) = sqrt(alpha beta) (P # Q).

[proved] Equation (2) is equivalently

    M(G) = G^(-1) # (J G J),                           (4)

because (G^(1/2)JG^(1/2))^2=G^(1/2)JGJG^(1/2).
In particular M(G) is the unique positive solution of

    M(G) G M(G) = J G J.                              (5)

Suppose

    (1-delta) G0 <= G <= (1+delta) G0,  0<=delta<1.

Inverse order reversal and congruence with J give

    G0^(-1)/(1+delta) <= G^(-1) <= G0^(-1)/(1-delta),
    (1-delta) JG0J <= JGJ <= (1+delta) JG0J.

Apply joint monotonicity and homogeneity to (4). The result is the
dimension-independent estimate

    kappa M(G0) <= M(G) <= M(G0)/kappa,
    kappa = sqrt((1-delta)/(1+delta)).                 (6)

Consequently kappa Delta_0(ell)<=Delta_G(ell)<=Delta_0(ell)/kappa
simultaneously for all real ell, hence for the entire allowed integer
vertex lattice. No finite enumeration or cutoff is needed for this step.

The bound is sharp even for one mode: G0=I and
G=diag(1+delta,1-delta) give M=diag(kappa,1/kappa).

[proved] If every surviving non-Gaussian vertex has Delta_0>=3, then

    delta < 5/13

implies kappa>2/3 and Delta_G>2 for every such vertex. A convenient
strict choice is delta<=1/4, which gives

    kappa=sqrt(3/5),   Delta_G>=3 sqrt(3/5)=2.3237900077...

The energy matrix is automatically positive under this relative bound.
This does not exclude exactly marginal quadratic current couplings; those
are precisely the entries already included in G.

## 3. Applying the bound to mobile cells

[proved, coordinate calculation] Let a=b+nu, f_alpha=k_alpha/pi, and

    T=[[a,0],[-f,I]],
    G0=diag(T^T V T, T^(-1) V T^(-T)),  V=diag(v0,...,vN).

The canonical transformation

    S=diag(T^(-1),T^T),   S J S^T=J

sends G0 to diag(V,V). Equation (2) therefore gives

    M0=diag(K,K^(-1)),  K=T^(-1)T^(-T).

In particular K_marker=1/a^2. The very small physical marker entry
does not invalidate (6): delta is measured in the full relative energy
norm, not an unweighted norm of microscopic or physical-coordinate
couplings. The same relative bound can be evaluated stably in the
squeezed canonical coordinates; its generalized eigenvalues are invariant
under this congruence.

[proved, conditional on C13's local-operator dictionary] A physical
marker harmonic has integer label D. For payload chiral charges q let
S_pay=(p.q)/(2 pi). The marker comb and payload translation phase give

    physical momentum = 2 pi(D+S_pay)/a mod 2 pi,
    Delta_0 = (D+S_pay)^2/a^2 + |q|^2/2.

At an exactly uniform resonance with integer reciprocal harmonic h,

    D+S_pay=h(b+nu),
    p.q-2 pi h nu=0 mod 2 pi,
    Delta_0,res=h^2+|q|^2/2.                          (7)

The integer hb is absorbed into D. In particular, large b produces many
small-dimension *oscillatory* physical harmonics, but not extra relevant
uniform resonances. Treating D as small independently of b at resonance
would give the wrong answer.

[proved, using graph girth six] Species neutrality means Wz=0, where
z is the total flavor-change vector. A nonzero integer circulation is
a sum of oriented cycles, so ||z||_1>=6. With chiral m=z modulo 2,

    Delta_pay=(||m||^2+||z||^2)/4 >= ||z||_1/2 >= 3.

When z=0, m=2u and Delta_pay=||u||^2. Thus the only potentially
dangerous non-Gaussian classes at a resonance are individual payload
backscatterings, pairs of them, pure length h=+/-1, and one length
backscattering times one payload backscattering. The target

    k_alpha=pi/4+alpha pi/[12(N+1)],
    nu in [5/12,7/12]

detunes these classes by the already stated cell-momentum bounds
pi/[6(N+1)] and pi/6. Every other resonant vertex has Delta_0>=3.
Equation (6) makes the required Gaussian margin quantitative.

**The filling caveat matters.** When equilibrium fillings move, both the
translation character and the free reference geometry move. One may
evaluate the certificate against the free reference at the actual dressed
fillings, together with the actual nonresonance margins, or retune the
fillings and then use the original reference. A small normalized stiffness
perturbation alone is not a bound on the filling shifts. Strict margins
and continuity give a neighborhood at fixed N and b, but that is a separate
input. The certificate should not be quoted without it.

## 4. Independent checks and what remains unresolved

[numerical] Run

    OPENBLAS_NUM_THREADS=1 python3 contrib/digit_compilation/audit/gaussian_bound_checks.py

with NumPy installed. The script writes `gaussian_bound_results.json` and
checks 60 full random G matrices with 1, 2, 7, 22 and 53 modes, relative
deformations from 0.05 to 0.6, and velocity ratios through 10,000. These
include nonzero generic C and noncommuting blocks. It checks the Riccati
identity (5), both inequalities (6), four exactly saturating one-mode
cases, and fifteen one-mode drift cases approaching the positivity
boundary. Nine physical-coordinate cases reach b=10,000 with 21 payload
modes. Exact rational arithmetic enumerates all low-dimensional neutral
payload classes at N=1,6,21,52,186 and checks the mixed margins and the
large-b cancellation in (7). All tests passed in the recorded run.

[audit verdict] This resolves an avoidable vagueness in the Gaussian
stability argument; it finds no infrared obstruction in the stated
selection rule. It does **not** supply a complete compiled high-channel
microscopic core, a local massive-sector elimination, or a measured
central-charge violation. The largest unresolved physics input remains
the local elimination in `research_continuation/IR_AUDIT.md` Section 9.
Its geometric connected-excursion bound is explicitly a controlled
physics assumption, not an established convergence estimate. The
defect-counting bound by itself does not establish the all-order linked
expansion with moving cell projectors. A complete microscopic composition
is also still needed before its hypotheses and constants can be evaluated.

The new result gives that eventual calculation an explicit target: certify
the effective G relative to the correct moving free reference, preserve
the finite momentum margins, and enter delta<5/13 (or the simpler 1/4
budget). It should not be recorded as a counterexample by itself.
