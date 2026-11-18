# Moving only the conserved charge field outside the contact embedding

2026-11-18. [Proved algebraically, with literal-bond checks; not independently
reviewed.] This is a proposed physical modification of the d=10 model. It
does not replace the head chemical potential or charged masses by genuine
on-site terms.

## 1. Charge conservation is exact

On the ten physical letters define

    q(L_a)=q(R_a)=2a, q(W)=q(V)=0,
    q(A)=q(XA)=-1, q(B)=q(XB)=+1,
    Q=sum_i q_i.                                        (1)

Each literal off-diagonal rule preserves the two-site total charge:

| Process | Charge on both sides |
| --- | --- |
| L_a H <-> H R_a | 2a+q(H) |
| y s <-> y X_s; s y <-> X_s y | q(y)+q(s) |
| A R_1 <-> B R_0 | -1+2=+1+0=1 |
| W W <-> X_A B; W W <-> X_B A | 0=-1+1 |
| a V <-> V a | q(a) |

The diagonal terms commute with Q. Therefore

    [h10, q tensor 1+1 tensor q]=0,
    [H10,Q]=0                                          (2)

for every parameter value, including nonzero port and conversion. The
occupied number and head number are conserved by the same literal rules.
The charge field is not an approximation to a symmetry of the effective
band: it is an exact symmetry of the microscopic Hamiltonian.

The ORIGINAL contact field is NOT itself the global charge operator. Its
many-site expression is

    H_field,contact = (J f/4)sum_i(q_i n_(i+1)+n_i q_(i+1)). (3)

It commutes with Q but generally not with H10 at f=0. For example, on a
four-site ring, the vacancy hop

    |L1,V,L0,V> -> |V,L1,L0,V>

preserves Q=2 but changes (3) from zero to Jf/2. Its commutator with the
vacancy hopping therefore has a nonzero matrix element of magnitude
J|f|t/2. Commutation with total charge does not make (3) a uniform field.

## 2. The proposed bond and its exact effect

Let b9^(0) be the parent's frozen bond with f=0. Replace ONLY the charge
field by

    h10,external = J embed(b9^(0)) - t sum_a(|a V><V a|+h.c.)
        -(mu/2)(n tensor 1+1 tensor n)
        +(J f_ext/4)(q tensor 1+1 tensor q).             (4)

On a periodic chain the last line sums to J f_ext Q/2. The factor 1/4 in
the bond is essential: each site belongs to two bonds. This Hamiltonian
still has physical d=10, range one, Hermiticity, and one-site translation
invariance. The squeezed grammar, cell lengths, and head number remain
conserved. No parameter depends on physical chain length.

Within every fixed-Q sector, (4) differs from its zero-field Hamiltonian
by the EXACT scalar J f_ext Q/2. In the vacuum logical space Q=2N_1,
where N_1 counts occupied data bits. Hence the logical field is exactly

    J f_ext N_1.                                       (5)

It has no dependence on bit position, vacancy placement, head position,
or J/t. Since all virtual microscopic processes preserve Q, a
charge-respecting effective isometry, if constructed, carries (5) exactly,
including finite g and lambda. This algebra does not construct a dressed
code or prove that such a low-energy isometry exists. A virtual denominator
inside one fixed-Q block
receives the same field on its initial and intermediate states and is
unchanged. Thus the referee's order-(J/t)f position-dependent field error
is removed identically, not merely reduced at leading order.

A literal n=3, N=4,V=3, J=1,t=25 check also reproduces the original
contact-field variation across the three one-bit positions; the values
are stored in `checks.json`. For the external field the three shifts
are exactly the same scalar (5). The finite free contact factor and
the thermodynamic factor used in (7) are reported separately in that
check, rather than adjusted to conceal their finite-size difference.

This does NOT make the full effective Hamiltonian free or uniform in
flavour: the field-independent dressed hoppings, self-energies, and
many-body interactions still need control.

## 3. Correct normalization of the intended band field

At fixed thermodynamic occupied density rho, the free contact factor is

    alpha_infinity(rho)=rho-sin^2(pi rho)/(pi^2 rho).      (6)

To leading order in J/t, all embedded material terms have the common
scale J alpha_infinity. If the frozen certificate calls for the
dimensionless physical field f9 in b9, choose the SIZE-INDEPENDENT value

    f_ext=alpha_infinity(rho) f9.                       (7)

Simply setting f_ext=f9 would overstate the normalized logical field by
1/alpha_infinity. At half occupancy alpha_infinity=0.29735763271532445...
At finite N,V, the leading compression has alpha_(N,V) in PROOF.md (9),
so the field in those finite normalized units is
alpha_infinity f9/alpha_(N,V). That O(L^-2) difference must not be hidden
by choosing a new microscopic coefficient separately at every size.
At nonzero fixed J/t, matching the actual dressed hopping scale is a
further task; (7) matches the leading frozen target only.

For the n=9 certified schedule, write

    E_* = lambda^2 alpha_clock^4 beta^2 g^2,
    alpha_clock^4=(2/11)^2 sin^4(pi/11),
    beta=1/d_10, d_0=1, d_1=1004, d_(k+1)=1004 d_k-d_(k-1),
    g=2^-163, lambda=2^-359, f9=E_* 2^-83.              (8)

The logical unit in the mobile leading Hamiltonian is J alpha_infinity
E_*. The band's allowed 2^-35 is a tolerance in its certified normalized
matrix/metric convention. It is not a selector energy density and not
an absolute microscopic allowance of J*2^-35. Removing (3) eliminates
one identified field error but supplies no bound on the other dressed
matrix errors in this normalization.

## 4. Selector and charged-sector consequences

For vacuum words and f_ext>=0, all-zero bits minimize (5). Their
Hamiltonian and ground-sector comparisons are exactly the f=0 problem.
For comparisons at the same fixed Q, the field cancels exactly even
when head counts differ. Thus the finite tests and the relative-response
problem in PROOF.md are not damaged by (4).

This assertion alone does not exclude negatively charged token sectors.
Their charge can be negative, so their energy can decrease in a positive
field. PROOF.md §6 proves zero-field charged domination under its stated
mass condition and supplies a finite-word port/field allowance. A uniform
charged gap and a finite-density, nonzero-port defect theorem remain
separate. Nor does the modification establish a compact unitary CFT or
its Virasoro central charge.

## 5. Literal complexity and coefficient hierarchy

There is one replaced diagonal interaction, no new physical state, no
extra independent field coefficient, and no new off-diagonal process.
It would nevertheless be misleading to say that the physical bond has
fewer literal terms.

At generic nonzero coefficients in the specified ten-letter matrix-unit
basis, with K>0 and mu removed in the canonical ensemble:

| Bond | Hermitian off-diagonal pairs | Nonzero diagonal matrix units |
| --- | ---: | ---: |
| Original contact-field bond | 50 | 78 |
| External-charge-field bond | 50 | 90 |

The 50 pairs are 10 material head moves, 14 left rail switches, 14 right
closing switches, 1 conversion, 2 ports, and 9 material/vacancy swaps.
Each counts an operator plus its Hermitian conjugate; counting directed
matrix units gives 100 off-diagonal entries. The field adds support on
the 12 ordered charged-material/vacancy diagonal pairs. This basis count
includes tied coefficients; it is not the number of independent real
controls. With the phases fixed to the schedule, sw=1, fixed N, and t
as common scale, the general amplitudes are J/t, nu, Delta, D, K, t_d,
g, lambda, and f_ext: nine real controls before imposing their schedules.
Allowing the three README phases to vary adds three controls. Restoring
the occupied chemical potential adds one in a grand-canonical problem.

The small microscopic alphabet does not remove the hierarchy. At n=9,

    log10 beta = -30.0173332505...,
    log10 f9 = -402.9768488466....                       (9)

At half occupancy the largest material diagonal is at least JD, whereas
a charged/vacancy field diagonal has magnitude J f_ext/4. Consequently
the ratio of nonzero literal bond entries is at least

    4D/f_ext > 10^407.10.                               (10)

This lower bound has already removed the occupied scalar and a common
positive energy scale. The zero VV and uncharged-vacancy entries fix a
convention with no irrelevant identity term; the field is in the
symmetric bond convention with its antisymmetric telescoping freedom
removed. Removing a gradient component from K Pi_bad does not affect
the charged-head/vacancy entries used in (10). Counting the field
as six one-site projectors instead changes factors of two, not this
enormous spread. This is a microscopic coefficient hierarchy, distinct
from small emergent velocities. The rational selector nu=1.77 and
moderate selector test J/t=1/25 do not make the certified band schedule
a simple Hamiltonian.
