# C85 result

[Proved] The eleven-letter grammar needs three vacancy colors, not two. The boundary pin distinguishes L_1 from L_0 and the markers, while R,D,X retain a third outgoing state. Three colors are sufficient and minimal for a nonnegative nearest-neighbor pair-penalty repair admitting every legal material word at arbitrary vacancy placements. This restricted minimum does not exclude redesigned grammars.

[Proved finite ingredients; assumed infrared transfer] The lift has d=11+3=14. A regenerated C84 certificate at n=13 active flavors, p=14 data sites and cell length b=15 meets n+1=14>d-1=13. Its carrier retains 29 vertices, including the inactive final data endpoint, and the product resolvent has dimension 841. This matches C83's local dimension rather than improving it.

| Certified finite quantity | Result |
| :-- | :-- |
| Simple Fermi points | 26 |
| Distinct unordered pair sums, including repeats | 351 |
| Root separation | greater than 0.0156742281 |
| Pair-sum separation | greater than 3.2513114e-6 |
| Additional normalized Bloch allowance | 2^-35 |
| Finite-g error | less than 1.727e-37 |
| Physical finite-K error, including zero-port splitting | less than 3.024e-83 |
| Local two-cell port remainder | less than 2.396e-54 |
| Length-15 selector density margin | greater than 2.8116274e-5 |

[Definition] The material schedule is D=1000, Delta=4-2cos(pi/16), phi_A=pi/28, phi_B=0, psi_A=pi, psi_B=pi/2, t_d=beta*2^-66, g=2^-272, lambda=2^-608, K>=2^2560, and physical field lambda^2 alpha^4 beta^2 g^2 2^-120. Here alpha=sqrt(2/16)sin(pi/16) and beta=1/d_15 with d_0=1,d_1=1004,d_(j+1)=1004d_j-d_(j-1). The midpoint selector and additional mobile K conditions are specified in [PROOF.md](PROOF.md). [MEMBER.json](MEMBER.json) records the exact definitions.

[Proved] Legal coloring is unique and N_material_bad<=N_bad. Projected swaps have amplitude -t without a multiplicity factor. The exact charges are occupied number, marker number and payload Q. The lift adds no independent color mode.

[Proved] The complete leading grammar correction on the legal space is

    -(t^2 C+J^2 A+2J^2 lambda^2 Z)/(JK).

C counts color-ambiguous occupied-vacancy contacts. Its coefficients are 1 on L_0,L_1,D,X,E, 2 on R_0 and 0 on R_1,M. It produces a t^2/(JK) tail-bit field and its coupling to vacancy and clock fluctuations, with no lambda^2 suppression. A is the inherited vacuum-tail bit self-energy, of scale J/K. Z is the wrong-end port self-energy, of scale 2J lambda^2/K. The mixed t/K operator vanishes because its two types of illegal intermediate are disjoint. All higher color-dependent fields, hopping and interactions are retained in the explicit Feshbach formula and the fixed-order local remainder. PROOF Section 4 lists the channels and their scales.

[Numerical] The leading correction formulas hold on every legal nonempty ring of lengths 3,4,5, up to 14497 legal states, with maximum discrepancy below 4.5e-16. ED with all four two-bit strings independently resolves the color-induced energy shifts.

[Numerical] A nonzero-bit two-cell closure of dimension 3348 has illegal material weights 8.30676e-5, 8.39127e-7, 8.39976e-9 at K=100,1000,10000. The uncolored closure instead approaches 0.0229959. Nonzero-port Q=2 tests with one and two vacancies, dimensions 368 and 1945, also show K^-2 suppression. All these eigenpair residuals are below 2.7e-11. A J=.01 test uses the actual gap JK and shows the same suppression.

[Proved] F3 transfers to full frozen marker sectors at nonzero port. The new last-L_1 penalty belongs to the left cell, while the marker-to-first-letter penalty belongs to the right cell. The dressed frozen defect bound is k(gamma-5sqrt(2)|lambda|)>0. F4 transfers with kappa=JK and a volume-independent fixed-order local remainder. Both have separate small-matrix checks.

[Assumed; heuristic conclusion] The candidate count is c=14>13. A microscopic central-charge result still needs finite-J matching, a uniform mobile dressed-defect comparison, multichannel stability and conformal-observable control. No central charge was measured. The material value K=2^2560 alone does not control the mobile model at arbitrary small J; the t^2/(JK) bit field must be below the much smaller J times logical scale. The required extra inequalities are explicit, but no full numerical infrared window is claimed.

Reproduce every check with `python3 -B run_all.py` from this folder. All dependencies and outputs are local. No files outside work_C85 were modified.
