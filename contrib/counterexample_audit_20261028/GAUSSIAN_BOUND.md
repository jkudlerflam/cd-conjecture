# Full Gaussian stability bound, including tilted cones

Status: **[proved conditional continuum statement]**. This supplies the
matrix inequality used by the d=24 argument for all current couplings,
including density-current mixing and unequal left/right velocities. It
does not establish the microscopic matching coefficients or the presence
of the 24 modes; those require the separate microscopic and vacancy audits.

## 1. Scaling metric

Use canonically normalized chiral fields with

    K = diag(I_p,-I_p),
    S = (1/4pi) integral [i partial_tau phi^T K partial_x phi
                         + partial_x phi^T V partial_x phi],
    V = V^T > 0.

For the 24-boson construction p=24. Choose a reference V_0>0 that commutes
with K. It can have arbitrarily different positive velocities; no lower
bound uniform across constructions is needed. In these coordinates the
free scaling dimension of a vertex labelled by real vector l is
Delta_0(l)=l^T l/2. The allowed compactification lattice, including marker
gluing, restricts l separately.

Let A=V^(1/2) K V^(1/2), and diagonalize this real symmetric invertible
matrix as O Lambda O^T. The transformation

    phi = V^(-1/2) O |Lambda|^(1/2) chi

makes the kinetic signature diag(sign Lambda) and the velocity matrix
|Lambda|. Sylvester inertia gives exactly p right and p left chiral modes.
The equal-time covariance metric and vertex dimensions are therefore

    G = V^(-1/2) |A| V^(-1/2),
    Delta(l) = l^T G l/2.                              (1)

This argument allows independent speeds on both chiralities. Changing
those speeds or mixing the currents while V remains positive does not
change the oscillator count. Each pair contributes one to the usual
nonchiral central-charge convention, so c=p; this is not c_R+c_L=2p.

## 2. A relative bound independent of velocity ratios

Since |A|^2=A^2,

    G V G = K V K.

The unique positive solution of X V X=K V K is the matrix geometric mean

    G = V^(-1) # (K V K),
    A_1 # B_1 = A_1^(1/2)
                (A_1^(-1/2) B_1 A_1^(-1/2))^(1/2) A_1^(1/2).

For completeness, its monotonicity in each argument follows from the
positive integral representation

    A_1 # B_1 = (1/pi) integral_0^infinity
                  s^(-1/2) (A_1^(-1)+s B_1^(-1))^(-1) ds.

If either positive argument increases in Loewner order, its inverse
decreases, and the displayed integrand increases. Thus the mean is
jointly monotone. It is also homogeneous in scalar factors with weight
one half for each argument.

Assume the complete retained Gaussian perturbation satisfies

    (1-eta)V_0 <= V <= (1+eta)V_0,       0<=eta<1.       (2)

Inversion, conjugation by K, and [V_0,K]=0 give

    V_0^(-1)/(1+eta) <= V^(-1) <= V_0^(-1)/(1-eta),
    (1-eta)V_0 <= K V K <= (1+eta)V_0.

Apply monotonicity to these pairs and use V_0^(-1)#V_0=I:

    sqrt((1-eta)/(1+eta)) I <= G
       <= sqrt((1+eta)/(1-eta)) I.                     (3)

Consequently every vertex, regardless of charge or polynomial degree,
obeys

    Delta(l) >= sqrt((1-eta)/(1+eta)) Delta_0(l).        (4)

If every allowed non-Gaussian uniform vertex has Delta_0>=3, then
eta<5/13 keeps all such vertices strictly above dimension two. The
campaign's eta<=1/4 leaves Delta>=3 sqrt(3/5)>2.32379.
The bound does not require a common velocity, reflection symmetry, or
a Gaussian Hamiltonian without cross density-current terms.

## 3. Scope in the counterexample audit

The strict inequality controls the entire vertex lattice, not only a
finite list of six-fermion terms. Lower-degree nonforward terms still
need a nonzero momentum mismatch, and quadratic forward terms are
retained in V. Exact microscopic U(1) symmetries and locality still need
to determine the allowed lattice. A finite-order local perturbation
bound is not by itself (2): its coefficients must be compared with the
smallest retained velocity during matching. The nested parameter choices
are designed to make that normalized error arbitrarily small.

This proof is the coordinating auditor's independent derivation,
2026-10-28, OpenAI Codex GPT-6 family. No human review is claimed.
