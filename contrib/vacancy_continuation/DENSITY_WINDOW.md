# A density window without chemical-potential tuning

Status: **[proved]** finite-volume energy and density bounds;
**[proved conditional on the bosonization dictionary]** resonance bounds;
**[physics inference]** their application to the proposed d=24 critical phase.
No central charge is measured or proved by the density estimate.

## 1. Statement for an arbitrary internal bond

Let h be any Hermitian nearest-neighbor bond on q internal states, with
M=||h||. Add a vacancy state 0 and define the periodic physical chain

    H = -t sum_(x,a) (|0,a><a,0| + |a,0><0,a|)
        + J sum_x h_lift(x,x+1),
    a=1,...,q; t>0, J>=0.

Here h_lift vanishes unless both sites are occupied. There is **no occupied
particle chemical potential**. All coefficients are independent of L>=3.
For any ground-state occupied-number sector N, put rho=N/L. Then

    |rho-1/2|^2 <= 3 pi J M/(16t) + pi^2/(32 L^2).       (1)

In particular, if J M/t <= 1/1000, every thermodynamic ground-state density
satisfies

    19/40 < rho < 21/40 < 1/sqrt(3).                     (2)

The conclusion holds for arbitrary q and arbitrary internal interactions,
including frustrated, asymmetric, or very large-norm bonds. It does not
assume a conformal internal chain, a gap, a product spin ground state,
or a perturbative expansion. Equation (1) controls a density window; it
does not prove differentiability of the equation of state or exclude
coexistence of two densities within that window.

## 2. Finite-volume proof

The minimum vacancy-hopping energy at occupied number N is

    E_kin(L,N) = -2t sin(pi N/L)/sin(pi/L).               (3)

To obtain the lower bound, order the occupied positions. Internal cyclic
translation decomposes the boundary hop into twisted hard-core-particle
sectors. The lowest N free-fermion levels on such a ring form a consecutive
sea. Their sum is bounded below by (3), irrespective of the twist.
The completely polarized internal word realizes the untwisted hard-core
sector and attains (3), including the usual fermion parity in the boundary
condition. Thus (3) is the exact minimum over the full internal space.
It also gives zero at N=0,L.

Since -M n_x n_(x+1) <= h_lift(x,x+1) <= M n_x n_(x+1)
and sum_x n_x n_(x+1) <= N as operators in that sector,

    E_kin(L,N)-JMN <= E(L,N) <= E_kin(L,N)+JMN.          (4)

For the upper bound any kinetic ground state is a valid trial state;
it need not minimize the internal interaction. Set N0=floor(L/2),
rho0=N0/L, and A_L=2t/[L sin(pi/L)]. Comparing a minimizing sector to N0
gives

    A_L[sin(pi rho0)-sin(pi rho)] <= JM(rho+rho0)
                                  <= (3/2) JM.          (5)

Write delta=rho-1/2. We have A_L>=2t/pi,

    1-sin(pi rho0) <= pi^2/(8L^2),
    1-cos(pi delta) >= 4 delta^2,   |delta|<=1/2.

For the second inequality use the chord bound
sin x >= (2 sqrt(2)/pi)x on 0<=x<=pi/4 and
1-cos(pi delta)=2 sin^2(pi delta/2). Substitution into (5)
proves (1). To prove (2), use pi<22/7:

    3 pi/16000 < 1/1600 = (1/40)^2,
    3(21/40)^2 = 1323/1600 < 1.

If a ground state superposes degenerate occupied-number sectors, (1)
holds for each sector in its support, since occupied number is conserved.
The construction does not require exact uniqueness of the finite-ring
ground state.

## 3. Safe resonance margins throughout the window

Use the material-cell operator dictionary established in the campaign's
[infrared audit](../counterexample_campaign/infrared/README.md). For a
neutral payload vertex of cell momentum q0 and free dimension Delta_pay,
and any integer cell harmonic ell,

    k_phys = rho(q0+2pi ell)/b mod 2pi,
    Delta_ref = Delta_pay + [(q0+2pi ell)/(2pi b)]^2.      (6)

At a uniform physical resonance k_phys=2pi m,

    Delta_ref,res = Delta_pay + m^2/rho^2.               (7)

For every rho in (2), m!=0 contributes at least
(40/21)^2=1600/441>3. For m=0 only the original zero-momentum
payload selection problem remains. Thus a nonresonant fermionic payload
whose only dimension-at-most-two uniform operators are Gaussian currents
still has no new relevant or marginal non-Gaussian uniform vertex.
All remaining non-Gaussian resonances have reference dimension at least
three. This includes fractional marker harmonics.

Let delta_int>0 be the smaller distinct-root and distinct-pair-sum
momentum margin of the payload, with delta_int<=2pi. For every
non-Gaussian reference vertex of dimension below three, (6) implies

    dist(k_phys,2pi Z) >=
       min[19 delta_int/(40b), 2pi(1-21 sqrt(3)/40)].     (8)

Indeed its unwrapped |k_phys| is less than 2pi rho sqrt(3).
This bounds the distance to any nonzero reciprocal vector. Its distance
to zero is bounded by rho delta_int/b for a nonforward payload term,
or 2pi rho/b for a nonzero pure marker harmonic. Forward Gaussian terms
are retained. Derivatives do not reduce dimensions. The inequality
controls the reference dictionary; dressed root shifts must still fit
a reserved fraction of the margin.

For the inherited n=23 C33 certificate, b=24 and
delta_int>3.1558e-7 give

    delta_phys > 6.24585e-9.                             (9)

The reference charge Luttinger parameter in gap coordinates satisfies
K_gap=rho^(-2)>1600/441. Its contact harmonics remain irrelevant for
both bulk locking (K_gap>2) and a static impurity (K_gap>1).
The Gaussian matrix condition eta<=1/4 used in the campaign keeps all
reference-dimension-three vertices above 3 sqrt(3/5)>2, uniformly over
this density window. The large material period changes coefficients and
the momentum cutoff, but not these strict dimension inequalities.

## 4. Consequence and limits for the proposed counterexample

First fix the internal h23 with its small-g, small-lambda hierarchy and
positive material-defect costs. Then choose J positive and sufficiently
small relative to t, including J||h23||/t<=1/1000 and the stronger
payload-matching tolerances. This is compatible with all preceding
choices because ||h23|| is finite and independent of L. Set the vacancy
chemical potential to zero. The density automatically remains in (2).

**[physics inference]** Under the same local charge-payload matching,
defect stability, and Gaussian phase-stability assumptions as the
[d=24 synthesis](../counterexample_campaign/COUNTEREXAMPLE_D24.md),
the expected c=24 conclusion therefore does not need an exact
half-filling adjustment. Half filling remains a convenient reference
point for numerical checks. The exact statement established here is
the density bound, not the matching or a new lattice-to-CFT theorem.
Bounded rounding to N=b times an integer remains a separate zero-mode
step and preserves the limiting window.

## 5. Checks and provenance

`density_checks.py` independently assembles sparse physical d=3 chains
with arbitrary Hermitian internal bonds, compares every occupied-number
sector with (3)-(4), and verifies (1). It also checks the rational
window/momentum constants and full bond symmetries. These finite checks
are diagnostics for the proof, not thermodynamic central-charge data.

Author: OpenAI Codex, GPT-6-family coordinator, 2026-10-28. The general
density estimate and tuning-free application were derived in this
continuation; the harmonic dictionary and C33 certificate are inherited
from the linked campaign. Independent reviews and other checks are in
this folder's README. No human review is claimed.
