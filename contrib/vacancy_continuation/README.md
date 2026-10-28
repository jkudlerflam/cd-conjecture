# Independent continuation of the vacancy counterexample

## Claim and status

**[proved; checked]** The vacancy construction does not require tuning
its particle chemical potential to exact half filling. For any fixed
internal Hermitian bond h, vacancy hopping t>0, internal coefficient
J>=0, and zero particle chemical potential, every thermodynamic ground
density lies in (19/40,21/40) if J||h||/t<=1/1000. This entire interval
has the required strict reference charge-harmonic margins.

**[checked physics inference]** An independent primary-source review
finds no concrete failure of the proposed one-site invariant,
nearest-neighbour d=24 chain's c=24 argument at the project's stated
theoretical-physics standard. A new bounded dynamical charge-response
calculation strengthens its matching argument. The conclusion retains
the inherited finite certificates and conventional Gaussian RG and
defect-matching inputs. It is not a rigorous lattice-to-CFT theorem or
a microscopic measurement of c=24.

**[checked finite construction]** C39's all-n reference-band argument
passes an independent analytical review. Exact root and pair-sum checks
now extend to n=4. This helps the separate fixed-d, unbounded-c effort;
it does not supply its full spontaneous selector or mobile completion.

## New results

1. [Density window](DENSITY_WINDOW.md). A volume-independent energy
   comparison gives

       |N/L-1/2|^2 <= 3 pi J||h||/(16t) + pi^2/(32L^2).

   It covers arbitrary internal bonds and odd physical rings. The
   inherited n=23 reference momentum margin remains above 6.24585e-9
   throughout the stated window. The
   [independent review](local_audit/DENSITY_REVIEW.md) checks the proof,
   10,710 twisted kinetic sectors, and 89,532 harmonic instances.

2. [Local perturbation and dynamical response](local_audit/README.md).
   The finite-order local Schrieffer-Wolff argument and auxiliary
   threshold pass primary-source inspection. In the Gaussian charge
   reference, G_rr(k,omega)<=1/e_gap'' uniformly in momentum and
   frequency, equal to 4/(pi t) at half filling. This controls the
   smooth influence-action coefficient at order J^2/t even for quantum
   internal sources. It is a continuum bound, not a lattice norm bound.

3. [All-n finite-input audit](all_n_audit/README.md). An independently
   derived positive-matrix-pencil argument validates the order of
   C39's parameter choices. A separate exact determinant and rational
   Cayley verifier checks the n=1,2,3 certificates and a new n=4
   reference: eight simple roots and all 630 pair comparisons, with
   delta4>1.0437e-54. The rational Gram evaluator is inherited.

The [full d=24 synthesis](../counterexample_campaign/COUNTEREXAMPLE_D24.md)
remains the main construction. These results strengthen specific steps;
they do not change its stated standard of evidence.

## Reproducing

From the repository root:

    python3 contrib/vacancy_continuation/density_checks.py
    python3 contrib/vacancy_continuation/local_audit/checks.py
    python3 contrib/vacancy_continuation/local_audit/density_review_checks.py
    python3 contrib/vacancy_continuation/all_n_audit/regenerate.py
    python3 contrib/vacancy_continuation/all_n_audit/verify.py

Only the first command requires numpy and scipy; the others use the
Python standard library. The density script checks all occupied-number
sectors of physical d=3 rings through L=7, using generic complex and
scalar internal bonds and four interaction strengths. Its largest
kinetic formula discrepancy is 5.4e-15 and eigenvector residual is
8.1e-12. It takes less than one second in this environment. The exact
certificate replay takes roughly ten seconds. Runtime varies by machine.
All generated outputs remain beside their scripts.

## Coordination and provenance

The user requested an independent continuation alongside the existing
agents under the same guidance. A GPT-6-family Codex coordinator and
three parallel research workers produced this folder on 2026-10-28.
The density estimate received a separate analytical review. The
coordinator reran the exact certificate verifier, the dynamical response
checks, and the physical density diagnostics. Existing collaborators'
sources were read without modification; each new worker used a separate
subfolder. No human review is claimed.

The ongoing C43 task elsewhere handles the all-n mobile completion.
This continuation deliberately supplies independent input checks for it.
