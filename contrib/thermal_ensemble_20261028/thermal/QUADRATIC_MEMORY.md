# Exact thermal matching for a quadratic retarded action

Status: **[proved finite-orbital theorem; numerical independent control]**.
This is an exactly solvable control of the thermal matching mechanism. The
C44 and vacancy actions have higher-degree interactions, so this theorem
does not apply to them without an additional RG argument.

Provenance: OpenAI GPT-6-family Codex thermal worker, 2026-10-28. AI-assisted;
no human review. The derivation uses finite Hermitian matrix algebra and
the Fermi occupation trace; no interacting universality input is used.

## 1. Hypotheses and the effective kernel

Let a one-dimensional translation-invariant, number-conserving quadratic
fermion Hamiltonian have smooth periodic Hermitian one-particle symbol

    h(k) = [ H(k)    B(k)^dagger ],
           [ B(k)       A(k)    ]

with q logical and r auxiliary orbitals, and A(k)>=Delta I for a fixed
Delta>0. A finite-range microscopic symbol is one sufficient smoothness
condition. Work at Fermi energy zero. Define

    F(k)=H(k)-B(k)^dagger A(k)^(-1)B(k),
    Z(k)=I+B(k)^dagger A(k)^(-2)B(k).                   (1)

Assume F has finitely many zeroes k_i, one-dimensional kernel at each,
and for a normalized u_i in that kernel,

    d_i=u_i^dagger F'(k_i)u_i != 0.                    (2)

No smallness of B is needed for the theorem. These explicit hypotheses
are what rule out annihilation or creation of Fermi crossings.

The exact Schur-complement kernel after Gaussian elimination of the
auxiliaries is

    S(z,k)=H(k)-z-B(k)^dagger(A(k)-z)^(-1)B(k)
          =F(k)-z Z(k)
             -z^2 B(k)^dagger A(k)^(-2)
                               (A(k)-z)^(-1)B(k).     (3)

For |z|<Delta, the last term is bounded in norm by

    |z|^2 ||B(k)||^2 / [Delta^2(Delta-|z|)].            (4)

The resolvent identity proves both formulas. In imaginary time the
positive auxiliary matrix has propagation exp(-tA), so its retarded
kernel has memory controlled by the auxiliary gap. Equation (3) also
shows why the frequency derivative cannot simply be dropped: it changes
the physical velocities even though it does not change the Fermi count.

## 2. Exact dressed velocities

**[proved]** The full Hamiltonian has exactly the simple zero crossings
k_i, with signed speeds

    v_i = d_i / [u_i^dagger Z(k_i)u_i].                (5)

Proof. Block Gaussian elimination gives

    det h(k)=det A(k) det F(k).

Since A is invertible, ker h(k) consists exactly of vectors

    w=(u,-A(k)^(-1)B(k)u),    u in ker F(k).           (6)

Their squared norm is u^dagger Z u. At a zero with one-dimensional
kernel, the full eigenvalue is differentiable in k. Hellmann-Feynman
gives its derivative as w^dagger h' w/(w^dagger w). Differentiating
F in (1), using (A^(-1))'=-A^(-1)A'A^(-1), makes the numerator exactly
u^dagger F' u. Equations (2) and Z>=I give (5) and transversality.
All other eigenvalues are separated from zero in neighborhoods of these
points; compactness supplies a nonzero lower absolute bound away from
the neighborhoods. This also excludes additional massless modes.

## 3. Thermal coefficient including the scalar factor

**[proved]** The thermodynamic free energy of the full chain satisfies

    f(T)-e = -(pi/12) sum_i |v_i|^(-1) T^2+o(T^2).     (7)

For completeness, diagonalize h(k) and write its eigenvalues E_a(k).
The Fock trace and the momentum Riemann sum give

    f(T)-e = -T integral_(BZ) dk/(2pi)
                          sum_a log(1+exp(-|E_a(k)|/T)).

Outside the crossing patches the integrand is exponentially small.
On a patch set x=E_a(k)/T. Its Jacobian tends to 1/|v_i| and is bounded,
so dominated convergence reduces the patch integral to

    -(T^2/(2pi|v_i|))
                     integral_R log(1+exp(-|x|))dx.

Expanding the logarithm and integrating its absolutely integrable series
gives the integral pi^2/6. Summing proves (7).

Gaussian auxiliary integration at finite temperature is an exact
determinant factorization with the original trace normalization retained:

    Z_full = Z_A Z_S.                                  (8)

One may define Z_S as the ratio on the right before removing a common
coherent-state time regulator; this avoids any independent choice of
normalization for an infinite Matsubara product. Equivalently, Gaussian
block elimination at each finite regulator proves the factorization,
and the physical trace fixes its continuum limit. The scalar factor
has free energy

    f_A(T)=-T integral_(BZ) dk/(2pi)
                            log det(I+exp(-A(k)/T)),
    |f_A(T)|<=r T exp(-Delta/T),       e_A=0.           (9)

Therefore the exact retarded action has the same T^2 coefficient (7),
and the eliminated positive modes contribute exponentially little.
Dropping the scalar without (9) would not be a justified general
procedure; here its size is proved from the complete trace.

The result is independent of the residue u^dagger Z u **after** that
residue is included in the speeds (5). Treating F alone as the logical
Hamiltonian gives the wrong thermal coefficient unless Z=I on its roots.

This is a theorem about an entire temperature-independent microscopic
Hamiltonian, not an arbitrary multiplication of a generating functional
by a temperature-dependent scalar.

## 4. One light band and one massive orbital

Take H(k)=-2 cos k-mu, A=Delta and B=eta real. Choose parameters with

    |mu+eta^2/Delta|<2.

There are two roots, cos k_F=-(mu+eta^2/Delta)/2, with speeds

    |v_F|=2 sin k_F/(1+eta^2/Delta^2),
    -(f(T)-e)/T^2 -> pi/(6|v_F|).                      (10)

The action is retarded for every eta!=0, yet the two crossings are one
nonchiral free mode. Their standard fermion oscillator algebra has c=1.
The residue renormalizes its speed without creating an additional
central-charge contribution.

`checks.py` independently diagonalizes the two-by-two symbol and performs
quadrature concentrated at its two exact roots. It compares the thermal
coefficient with (10), checks the full-matrix eigenvalue slope against
(5), and reports the incorrect coefficient obtained by discarding Z.
It also checks the scalar auxiliary pressure bound (9). Run

    python3 -B contrib/thermal_ensemble_20261028/thermal/checks.py

The generated `checks.json` is an algebra/thermodynamic control. It is
not evidence of a c-d violation or a multichannel interacting theorem.

## 5. Exact boundary of this control

The matrix calculation pays no inverse logical bandwidth in its auxiliary
elimination. Every denominator in (3)-(4) uses Delta. The low-temperature
coefficient, as it must, contains the inverse *dressed* light speed.

For an interacting action there is no single finite matrix determinant
whose zeroes exhaust the spectrum, and (8) becomes a connected vacuum
functional. Its scalar RG terms, including finite-temperature flow
differences, must obey a bound such as THERMAL_CRITERION.md (2)-(3).
The free calculation validates the matching interpretation and supplies
an exact test case; it does not prove those interacting bounds.
