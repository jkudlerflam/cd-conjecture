# A stronger translation-covariant normal form

Status: **[proved, under the finite-range Lieb-Robinson assumptions stated
in UV_NORMAL_FORM.md]**. This strengthens the finite-depth construction:
the remainder has cutoff power 10 instead of 18, and the change of basis
is translation covariant. The improvement was proposed by the coordinating
agent and verified and written out independently by the vacancy worker.

## Statement

Use exactly the bounded, finite-range hypotheses and definitions in
[UV_NORMAL_FORM.md](UV_NORMAL_FORM.md), Sections 1–3. In particular,

    H0=sum_x a_x,             V=sum_x v_x,
    s_x=integral_(-T)^T F_delta(u) tau_u^(I_x)(v_x)du,
    S_R=sum_x s_x,            T=delta^-2,
    R=O(delta^-2),            ||s_x||<=C/delta.

Every `s_x` is anti-Hermitian and supported in a radius-`R` interval.
The exact first-order identity already proved there is

    V+[S_R,H0]=W_delta+E_delta,
    ||E_delta||_1<=C delta^-5 exp(-c/delta^2).                 (1)

The Gaussian-filtered local operator `W_delta` is unchanged. Define the
single exponential

    U(s)=exp(s S_R).

There are constants `C0,C1,c>0`, independent of volume, such that, whenever

    q=C0 |s|R/delta <= 1/2,

the exact transformed Hamiltonian satisfies

    U(s)(H0+sV)U(s)*=H0+s W_delta+R_delta(s),

    ||R_delta(s)||_1
       <= C1 s^2 delta^-10
                   +C1 |s| delta^-5 exp(-c/delta^2).          (2)

Here the first spatial moment is the same explicitly decomposed norm

    ||Phi||_1=sup_y sum_(X contains y)(1+diam X)||Phi(X)||.

The nonlinear remainder has exponentially summable spatial tails on a
length scale `O(R)`. If `H0` and `V` are translation invariant, then `S_R`
and `U(s)` are translation invariant as well. More generally, translating
the coefficient array translates the normal-form interaction exactly.
There is no coloring choice.

In physical units, choose `s=J/t` and `delta=s^a`, `0<a<1/10`. Then
`q=O(s^(1-3a))` tends to zero, and

    ||t R_delta(s)||_1/J
       <=C s^(1-10a)+C s^(-5a) exp[-c s^(-2a)] -> 0.          (3)

This remains a Gaussian-filter ultraviolet lemma. It proves neither the
low-frequency ground-state comparison nor the c–d counterexample.

## Connected commutator count

Fix an original local term `B_x` supported within radius `r0` of `x`, with
`||B_x||<=M`. Take `R>=max(r0,1)`; changing the earlier fixed cutoff
constant ensures this.

Expand the `n`th iterated commutator as an ordered sum:

    ad_(S_R)^n(B_x)
       =sum_(y_1,...,y_n)
          [s_(y_n),[...,[s_(y_1),B_x]...]].                  (4)

Any nonzero summand has the following necessary condition: the interval
support of the next `s_(y_j)` intersects the union of the supports already
used. If it does not intersect, the next commutator is zero. This condition
may count additional summands that vanish for algebraic reasons, which
only makes the upper bound larger.

At step `j-1`, the union lies within radius

    r0+2(j-1)R

of the original anchor. A radius-`R` interval that intersects this union
has its center within radius `r0+(2j-1)R` of the anchor. Hence there are at
most `C j R` choices for `y_j`. This estimate is valid on finite intervals
and rings too, since a finite volume can only reduce the number of choices.
The total number of potentially nonzero ordered sequences is therefore

    product_(j=1)^n (C j R)=(C R)^n n!.                      (5)

Every commutator increases the norm by at most `2||s_y||<=C/delta`.
Combining this bound with (5),

    ||ad_(S_R)^n(B_x)||
        <= M n! (C R/delta)^n.                              (6)

All terms in (4), and therefore their sum, are supported in the fixed
anchored interval

    X_(x,n)={z:dist(z,x)<=r0+2nR}.                           (7)

This last fact matters: it permits a specified interaction decomposition
without converting an extensive operator bound into a local bound.

## The factorial cancellation and weighted norm

Assign the order-`n` BCH contribution from the original anchor `x` to
the support (7):

    b_(x,n)=s^n/n! ad_(S_R)^n(B_x).

Equation (6) gives

    ||b_(x,n)||<=M q^n,                                     (8)

after choosing `C0` in the definition of `q` sufficiently large.
For each fixed `n>=1`, a site `y` belongs to at most `C nR` supports
`X_(x,n)`, and each such support has diameter at most `C nR`. Consequently

    ||sum_x b_(x,n)||_1 <= C M n^2 R^2 q^n.                 (9)

If the supporting interval wraps all the way around a short ring, use
the whole ring as its support. Its size and diameter are still bounded
by the right side of (7), so (9) is uniform in all ring lengths.

For `q<=1/2`, elementary power-series estimates give

    sum_(n>=2) n^2 q^n <= C q^2,
    sum_(n>=1) n^2 q^n <= C q.                              (10)

Apply the first estimate to the transformed `H0` after subtracting its
zeroth and first orders. Apply the second to `sV` after subtracting its
zeroth order, retaining the prefactor `|s|`. The nonlinear remainder
therefore obeys

    ||R_nonlinear||_1
       <=C R^2(q^2+|s|q)
       <=C s^2(R^4 delta^-2+R^3 delta^-1)
       <=C s^2 delta^-10.                                  (11)

The last inequality uses `R=O(delta^-2)` and `delta<=1`. The term from
`sV` alone has the smaller power seven. Together with (1), (11) proves
(2).

An exponential spatial weight can also be included. The diameter in (7)
is at most `C_* nR` after increasing a fixed constant. Choose
`kappa=(log 2)/(2 C_* R)`. Then multiplication by
`exp(kappa diam X_(x,n))` replaces `q` by at most `sqrt(2)q` in (9).
At `q<=1/2` the series remains summable, with bounds by constants times
`q^2` and `q` as in (10). Thus (11) also holds, with changed constants,
in the norm with the extra factor `exp(kappa diam X)`.

## Equality with the actual unitary in every finite volume

In every fixed finite volume, the sums `S_R,H0,V` are bounded operators,
even when the on-site spaces are infinite dimensional, because they
contain finitely many bounded summands. The usual exponential series
therefore proves

    exp(s S_R) B exp(-s S_R)
       =sum_(n>=0) s^n/n! ad_(S_R)^n(B)

in ordinary operator norm. This finite-volume equality needs no
smallness assumption, although its crude extensive bound depends on
volume.

For `q<=1/2`, (9) supplies an additional absolute convergence estimate
in a local interaction norm that is independent of volume. For each
finite volume the anchored decomposition is consequently the same
operator as the exact exponential conjugation. No thermodynamic-limit
interchange is used to establish the equality. The uniform estimates
then define the associated infinite-volume local interactions by the
usual consistent local limits, if desired. The claimed finite-volume
normal form and defect estimates need only the finite-volume result.

## Covariance and the remaining obstruction

The cutoff intervals, restricted dynamics and kernel are chosen by the
same geometric prescription at every anchor. Under a lattice translation
`T`, therefore,

    s_x[T H0 T*,T V T*]=T s_(x-1)[H0,V] T*,

with the index shift adjusted to the translation distance. Thus the
entire anchored commutator construction transforms covariantly. In
particular, one-site translation invariance of the microscopic pair
`(H0,V)` implies one-site translation invariance of the transformed
Hamiltonian and its specified interaction decomposition.

This removes the artificial phase dependence of the earlier color
circuit. A cluster entirely inside a translated selected domain is the
translated reference cluster. Clusters that reach its interface can be
counted with the first spatial moment. This is a statement about the
constructed Hamiltonian interactions, not about connected ground-state
energies. Actual bulk cancellation in a sector-energy comparison still
needs the identification of the sector Hilbert spaces and the
low-frequency ground-state estimate.

The central obstacle is unchanged: `W_delta` is a low energy-transfer
filter, not a projection onto states of low absolute bath energy. Its
ground-state response with an arbitrary array of quantum cells has not
been bounded here. No claim of a completed vacancy defect inequality,
microscopic central charge, or counterexample follows from (2).
