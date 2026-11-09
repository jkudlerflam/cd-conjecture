# An exact binding criterion for the Fibonacci cusp

AI assistance: OpenAI Codex (GPT-6 family), 2026-11-09. Independent
parallel Codex derivations and numerical checks are identified below.
No human review is claimed.

## Claim and scope

**[proved]** In the number-conserving, second-order blockaded model of
C173, at the symmetric phase and fixed real-gauge total momentum K=0,
write the nearest hop as u, the next-nearest hop as +1, and the
interactions at separations 2 and 3 as V2=-2+v and V3=-g. For
0<u<4 and 0<g<2, a two-particle bound state lies strictly below the
global free two-particle threshold if and only if

    v < g u^2 / [2(2-g)].

At v=0 there is exactly one such eigenvalue in this momentum sector.
For the microscopic candidate's effective coefficients, g=2/3 and v=0,
the bound-state threshold is v_c=u^2/4. In particular, at u=1, raising
V2 from -2 to at least -7/4 is necessary and sufficient to remove this
K=0 bound state while the other effective coefficients stay fixed.
This does not exclude bound states in other total-momentum sectors.

This is a theorem about a two-body effective Hamiltonian. It invalidates
the unbound, dilute two-valley premise; it does **not** prove a
many-body central charge, rule out c=2 at density 0.17, or establish the
phase of the full nearest-neighbour qubit chain. The candidate's
microscopic bond has d=2 and one-site translation invariance; its
effective next-nearest hop is generated perturbatively, not added to
that microscopic bond.

## Relation to current work

While this calculation was in progress, upstream commit `1a4015f`
published [C175](../../codex/work_C175/RESULT.md), which independently
found the bound state and its zero-energy resonance condition. This
contribution credits that discovery. Its additional results are a
half-line square factorization, an exact necessary-and-sufficient
criterion with a proof by quadratic forms, and a finite integer-vector
certificate requiring only Python's standard library.

The independent [scattering calculation](../cusp_scattering_20261109/README.md)
supplies a complementary algebraic root, normalized wavefunction and
pair mass. The [research continuation](../counterexample_push_20261109/README.md)
records the resulting priorities.

## Proof

Let n=r-1>=1 index the relative separation r>=2. The real-gauge
two-body Hamiltonian is

    (H f)_n = 2u(f_{n-1}+f_{n+1}) + 2(f_{n-2}+f_{n+2})
              + (-2+v) delta_{n1} f_n - g delta_{n2} f_n,

with f_0=f_{-1}=0. The one-particle band is
epsilon(k)=2u cos k+2 cos(2k). For 0<u<4 its minimum is
-2-u^2/4, so E_th=-4-u^2/2.

Let A be the half-line adjacency operator, (Af)_n=f_{n-1}+f_{n+1},
f_0=0. Since A^2 has diagonal 1 at n=1 and 2 elsewhere,

    H-E_th = 2(A+u/2)^2 + v |1><1| - g |2><2|.                 (1)

Separating the first row of the square and completing its boundary
square gives, for every finitely supported f,

    <f,(H-E_th)f>
      = 2 sum_{n>=2} |f_{n-1}+(u/2)f_n+f_{n+1}|^2
        + (2-g)|f_2 + u/(2-g) f_1|^2
        + (v - g u^2/[2(2-g)]) |f_1|^2.                    (2)

Thus v>=v_c implies H>=E_th. This includes equality and establishes
nonbinding in this sector without finite-size extrapolation.

For v<v_c, choose f_1=1 and f_2=-u/(2-g), then solve

    f_{n+1} = -(u/2)f_n - f_{n-1}, n>=2.

Its characteristic roots lie on the unit circle and are distinct for
0<u<4, so f is bounded. Multiply it by a cutoff which is 1 through
n=R, decreases linearly to 0 over the next R sites, and vanishes
thereafter. The boundary square stays zero. The bulk residual is
O(1/R) at O(R) sites, so its squared norm is O(1/R).
Equation (2) is therefore v-v_c+O(1/R)<0 for sufficiently large R.
By the variational principle the spectrum extends below E_th.
The boundary terms have finite rank, and E_th is the bottom of the
essential spectrum of the unperturbed half-line band; this is a
discrete, normalizable bound state.

When v=0, (1) is a nonnegative operator minus a rank-one positive
operator, so it has at most one negative eigenvalue, counting
multiplicity. The preceding construction proves there is exactly one.

This also explains C175's resonance determinant:

    D = (V2+2)(4+2V3) + V3 u^2 = 2(2-g)(v-v_c).

The determinant's sign is now connected to the spectrum by (2), rather
than a continuity argument alone.

## Finite exact certificate and robustness

At u=1, g=2/3, v=0, `check_certificate.py` uses an explicitly listed
24-component integer vector and evaluates its Rayleigh quotient with
exact rational arithmetic. It proves

    E_pair <= E_trial < -4.51 < E_th=-4.5.

The vector is extended by zero on the infinite half-line. Its finite
support makes this an infinite-volume variational certificate, not a
claim that a truncated eigenvalue equals the infinite one.

For an effective-model perturbation with hopping changes delta_w and
delta_t2 and separation-potential changes delta_V2, delta_V3, operator
norm bounds give a two-body change at most

    4(|delta_w|+|delta_t2|) + max(|delta_V2|,|delta_V3|).

The free two-particle threshold moves by at most
4(|delta_w|+|delta_t2|). Consequently the same variational argument
certifies binding whenever

    8(|delta_w|+|delta_t2|) + max(|delta_V2|,|delta_V3|) < 1/100. (3)

Here the fixed test sector is K=pi in C173's original gauge (K=0 in
the real gauge). The threshold in (3) is the global free-pair
threshold, so unequal valley minima are allowed. Number-changing
perturbations and the remainder of the microscopic reduction are
outside this estimate.

## Checks and reproduction

Run from the repository root:

    python3 contrib/cusp_binding_certificate_20261109/check_certificate.py

No third-party dependencies are needed. The checker verifies (1) and
(2) with rational arithmetic for multiple parameters and finite
vectors, verifies the integer witness, and reports exact bounds.

The separate numerical check uses NumPy and SciPy:

    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
      python3 contrib/cusp_binding_certificate_20261109/numerical_check.py

It checks cutoff convergence against E_pair=-4.519526534666398 and
the threshold transition. Those numerical checks are supplementary;
the theorem and rational certificate do not depend on them.

Author and initial checker: root Codex agent, independent of C175.
Parallel algebraic and numerical cross-check: `cusp_scattering` Codex
agent. See the saved result files for actual rerun output and versions.
