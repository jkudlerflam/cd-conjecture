# Noncommuting contact-source diagnostic

Status: **numerical diagnostic passed**, including a two-impurity extension.
The tested finite systems show bounded O(J^2) charge dressing when the
source is an operator that does not commute with a machinery Hamiltonian
whose own gap is O(J). This checks a mechanism relevant to the campaign;
it is not a proof for the full d=24 bond or for a thermodynamic phase.

Independent calculation by the microscopic-audit GPT-6-family Codex worker,
2026-10-28. No human review is claimed.

## Hamiltonian and baseline

There are N periodic material-gap sites, N vacancies, and physical length
L=2N. The gap-chain hopping has unit amplitude, t=1, without square-root
occupation factors. Its basis dimension is binomial(2N-1,N-1). Let A_j be
the projector onto zero vacancies at material bond j. The single impurity is

    H = H_gap + J Delta/2 (1-sigma_x)
              + J u (A_0-alpha_L) sigma_z,
    Delta=0.2, u=1, J=0.04,0.02,0.01.

Here alpha_L is the finite-ring ground-state expectation of A_0. The
unperturbed impurity ground state is |+x> with energy zero and its excited
state costs J Delta. Thus the correct baseline is E_gap, with no added
J Delta/2. The source F=A_0-alpha_L satisfies <F>=0. In the sigma_x basis,
the exact coupled matrix after subtracting E_gap is

    [ H_gap-E_gap            J u F                 ]
    [ J u F                  H_gap-E_gap+J Delta   ].

The independent second-order expression is

    delta E_2 = -J^2 u^2 <psi|F
                  (H_gap-E_gap+J Delta)^-1 F|psi>.

It retains the impurity's excitation energy in the denominator. The source
is orthogonal to the charge ground state, so there is no artificial
1/(J Delta) contribution from that state. At N=4 and N=6, complete spectral
decomposition and an independent conjugate-gradient resolvent calculation
agree within 2.2e-16 in the susceptibility. For N=8,10 the same resolvent
formula is evaluated by sparse linear solves.

The independently diagonalized charge ground states also reproduce

    E_gap = -2/sin(pi/(2N)),
    alpha_L = 1/2 - 1/[2N^2 sin^2(pi/(2N))].

The finite-size centering is intentional. These tests do not claim that
their N-dependent centering defines a single fixed microscopic bond for
all lengths.

## Single-impurity results

At J=0.01:

| N | Bath dimension | Exact delta E/J^2 | Second-order delta E_2/J^2 | Relative error |
| --- | ---: | ---: | ---: | ---: |
| 4 | 35 | -0.07740672279 | -0.07740679660 | 9.54e-7 |
| 6 | 462 | -0.08474251008 | -0.08474261664 | 1.26e-6 |
| 8 | 6435 | -0.08834740386 | -0.08834752511 | 1.37e-6 |
| 10 | 92378 | -0.09047699771 | -0.09047712658 | 1.42e-6 |

The largest coupled space has 184756 states. At N=10 the coupling sweep is:

| J | Exact delta E/J^2 | Relative error versus delta E_2 | Dressed impurity gap | Bare gap |
| --- | ---: | ---: | ---: | ---: |
| 0.04 | -0.09012514611 | 2.17e-5 | 0.00799849244 | 0.008 |
| 0.02 | -0.09035949497 | 5.61e-6 | 0.00399981152 | 0.004 |
| 0.01 | -0.09047699771 | 1.42e-6 | 0.00199997644 | 0.002 |

All single-impurity eigenpair residuals are below 2.7e-14. The discrepancy
from the retained-gap second-order expression scales approximately as J^4
in this finite-size sweep. The coefficient delta E/J^2 increases moderately
in magnitude with N; these four sizes establish neither convergence nor a
uniform infinite-volume bound. No finite-size instability is observed.

## Two quantum impurities reveal an additional operator effect

To avoid relying solely on a single source, the script also couples
independent impurities at adjacent sites 0 and 1:

    H = H_gap + J Delta (Q_0+Q_1)
              + J u [F_0 sigma_z,0 + F_1 sigma_z,1],
    Q_a=(1-sigma_x,a)/2.

N=6,8,10 are tested; the largest coupled space has 369512 states. Bare
second-order ground-energy shifts add, because the two intermediate spin
flips are orthogonal. The exact answer also contains a charge-mediated
pair-flip interaction. This is an operator effect absent from a scalar
source test.

At fixed finite N define

    chi_01(0)=<psi|F_0 (H_gap-E_gap)^-1_perp F_1|psi>.

Eliminating charge excitations to second order gives a pair-flip matrix
element -2J^2 u^2 chi_01 between |+x,+x> and |-x,-x>. The latter state costs
2J Delta. Its leading extra energy correction is therefore

    -2 J^3 u^4 chi_01(0)^2 / Delta.

The measured correction beyond the sum of the two retained-gap
second-order energies has precisely this scale:

| N | Predicted coefficient of J^3 | (delta E-delta E_2)/J^3 at J=0.01 |
| --- | ---: | ---: |
| 6 | -0.0011456530 | -0.0011377589 |
| 8 | -0.0005527341 | -0.0005494544 |
| 10 | -0.0002979007 | -0.0002977493 |

For N=10 the measured coefficient is -0.0003012082, -0.0002982441 and
-0.0002977493 at J=0.04,0.02,0.01. The exact pair delta E/J^2 at J=0.01 is
-0.18095723066, versus -0.18095425316 at second order. Eigenpair residuals
remain below 5.1e-14.

This extension confirms why one should retain the quantum machinery
operators during matching: their induced interactions affect the answer.
The J^3 effect remains smaller than the leading J^2 dressing, and hence
does not contradict the proposed scale separation.

## Reproduction and limitations

Run from the repository root:

```sh
PYTHONPATH=/private/tmp/cd-campaign-python /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 contrib/counterexample_audit_20261028/microscopic/quantum_impurity/check_quantum_impurity.py
```

Use `--max-N 8` to skip the largest systems. The complete sweep takes about
five seconds in the audit environment and writes `results.json` and
`pair_results.json`. The code constructs the gap-chain basis and hopping
independently; it does not import the earlier scalar-defect toy.

These calculations test one and two fixed material-coordinate impurities
in a finite periodic charge chain. They do not establish bounds uniform in
the number of impurities, control an extensive payload with many slow
velocities, prove summability of a full operator-valued dynamical kernel,
or measure a central charge. Their positive result is that the simplest
noncommuting O(J) machinery/contact mechanism behaves as the claimed
O(J^2) dressing expects, with the expected additional quantum correction.
