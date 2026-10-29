# A twelve-state mobile candidate, and why it is harder than it looks

## Claim

A single uncoloured vacancy does suffice for the eleven-letter C84/C88 grammar, provided the vacancy is allowed only where the grammar automaton is in state 0. This gives an explicit one-site translation-invariant nearest-neighbour bond on C^12 (eleven material letters plus one vacancy V) whose legal sector at zero hopping is, cell by cell and bond by bond, the frozen C84 chain with the port removed on occupied bonds. With C84's certified n = 11 member the conditional count is

    c = 11 + 1 = 12 > 11 = d - 1.

The price is new physics that the two-colour chain C88 does not have. Because V cannot sit to the right of the walker, where a vacancy may stand depends on the clock state, so the vacancy fluid is no longer spin-charge separated from the material. Two consequences are proved or computed here. (i) The chain is mobile only in a slow-vacancy regime t << J, where a vacancy crosses a cell only when that cell's walker is parked at its right end, with amplitude t_eff = c_p t^(p+1) J^(-p). (ii) The walker mediates an interaction between the two gaps that flank a cell, and in the bare design it is attractive, X = -14.2 t^4/J^3 at p = 12, far larger than t_eff. Left alone this phase-separates the vacancies and would remove the extra boson, leaving c = n = d - 1. A contact energy w < 0 on the pair (V, R_0) makes the interaction repulsive at order |w| t^2, and with it there appears to be a consistent parameter window.

On the obstruction side, I prove that d = 12 is impossible for this grammar inside the class used by C82, C83 and C88 (every bond of every legal word admits a vacancy), and that of the three maximal single-vacancy grammars only the state-0 one gives a vacancy that is both mobile and blind to the payload bits.

The smallest mobile d achieved is therefore 12, as a conditional candidate whose mobile reduction is heuristic, one level weaker than C88's d = 13.

## Status

- [proved] The squeezing lemma. If aV and Vb are both legal then ab is legal, so a legal physical word has a legal squeezed material word, N_matbad <= N_bad. Colouring is trivially unique (one colour). Exhaustive on all 12^L rings with L <= 6 as well.
- [proved, exhaustive L <= 6] Vacancy hopping maps legal words to legal words except in one context, o(x) = 1 followed by E V, which the hop sends to x V E. This leak is bit blind (it depends only on whether the walker is at the right end). The vacancy-context switch never leaves the legal space. The pin count F of C88 is conserved by every literal transition, including illegal ones (0 violations in 7.7 million transitions at L = 6).
- [proved, exact identity checked numerically to 7e-15] At t = 0 and w = 0, the block in which every vacancy sits in an inter-cell gap equals the frozen eleven-letter ring on the squeezed word with the port deleted on each occupied bond, entry by entry, penalties included, with nonzero bits, tokens and port. Without the vacancy-context switch the identity fails at order one, because a vacancy in a gap cuts the folded carrier of the next cell.
- [proved, by inspection] On the F = 0 sector of the pinless grammar, H_C84 = H_pinless + K F gives H_pinless = H_C84, so C84's n = 11 certificate is the one that applies (not C85's n = 13).
- [numerical] Illegal material weight falls as K^-2 or faster with nonzero bits, vacancies and port; the transparent uncoloured control plateaus (1.29e-5 and 6.0e-3 in two of the closures).
- [numerical, exact linear algebra] Vacancy-fluid coefficients from Rayleigh-Schroedinger perturbation theory at p = 2 to 12, cross-checked against exact diagonalization at p <= 6 (agreement to 0.2 percent at the smallest t).
- [numerical] Payload dictionary at small size. The port-induced energy of a three-cell, one-vacancy ring is 2/3 of the frozen one (one bond in three is occupied) up to O(t^2), and the payload-specific one-fermion energy approaches 2/3 with an O(t^2) deviation of large coefficient.
- [heuristic] The slow-vacancy fluid is a gas of hard rods on the gap lattice (on-site and nearest-neighbour repulsion much larger than t_eff), hence a Luttinger liquid at vacancy filling 0 < nu < 1/2, with all payload-scale couplings far below t_eff. The parameter windows in MEMBER.json.
- [assumed] Everything C88 needs at the infrared level, (I1) to (I4) of contrib/small_d_note, transposed to this chain, plus the frozen pin selection (S_pin) at C84's schedule.

Nothing here measures a central charge. No probabilities are stated anywhere in this folder.

## How established and checked

Author: Claude (Opus 5.5) construction agent for Jonah Kudler-Flam, 29 October 2026. No Codex call; the work was repository-bound and numerical. Python 3 with numpy and scipy only. Every computation ran locally in minutes and under 1 GB, so no Typhon job was needed. Read before starting: contrib/small_d_note/note.tex, codex/work_C88 (RESULT, bond13.py, check_lift.py, MEMBER.json), codex/work_C84 (RESULT, MEMBER.json) and contrib/redteam_d13/README.md. Nothing outside this folder was modified; check_parent.py imports codex/work_C88/bond13.py read-only for comparison.

Tags. [proved] is an algebraic argument, an exhaustive finite enumeration, or an exact identity verified to rounding. [numerical] is floating-point evidence at the stated small sizes and diagnostic couplings, far from the certified schedule. [heuristic] is a physics argument not made controlled. [assumed] is an input used but not established.

## Construction

### Letters and grammar

The eleven material letters are C84's, L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E, M_A, M_B, with C88's automaton (outgoing state o, incoming set I) and no pin. The twelfth letter is a vacancy V with o(V) = 0 and I(V) = {0}. Every pair ab with o(a) not in I(b) costs JK. So V may follow only L_a, E, M_s or V, and may precede anything except R_1.

The squeezing lemma is one line. A legal aV forces o(a) = 0, a legal Vb forces 0 in I(b), and o(a) in I(b) is exactly legality of ab. The hidden pair R_0 V L_0 that killed the uncoloured vacancy of C75 and C80 is excluded because R_0 V itself is illegal. The cost is that vacancies are excluded from every bond to the right of a walker or token.

### The bond

With b11 the C88/C84 material hopping terms (walker, rails, conversions, t_d link, marker switch, port), o the material on-site energies of C84 and o(V) = 0,

    h12 = J b11 (material-material block only)
          - J sum_s ( |V X_s><V D_s| + h.c. )                  vacancy-context switch
          - t sum_{a in S} ( |a V><V a| + h.c. ),  S = {L_0, L_1, E, M_A, M_B}
          + J w |V R_0><V R_0|                                 vacancy-walker contact
          + (J/2)(o x 1 + 1 x o) + JK p_ext .

Three choices differ from C88 and each matters.

1. The vacancy-context switch. In C84 the D-to-X rail switch at the first data site is identified by the marker on its left. A vacancy in the inter-cell gap would otherwise cut that edge and change the carrier graph of the next cell at order one (frozen_identity.json, control entry error 1.0).
2. Genuine on-site energies. In C88 every material diagonal is a contact term, so a vacancy beside a token or marker removes half of D, Delta or mu_sel. Here a vacancy in a gap changes no diagonal energy, which is what makes the t = 0 identity exact.
3. The contact w < 0 on (V, R_0). V R_0 occurs only when a vacancy sits immediately left of the walker, and never while a token is in that cell (position 0 is then the token or an L). It is bit blind and invisible to every virtual carrier path, and it shifts the clock vacuum of a cell whose left gap is occupied. Its only purpose is the sign of the gap-gap interaction below.

h12 is Hermitian and conserves the occupied number N, the marker number and the C84 charge Q. At a generic point it has 199 nonzero entries, 62 off-diagonal (w = 0). The material block equals C88's bond13 material block to 1.4e-14 at eight random points.

### Why the regime must be slow

A vacancy may stand only on state-0 bonds. In a cell L^j T R^(p-1-j) E these are the j + 1 bonds from the previous marker up to the walker; the p - j bonds from the walker to E are closed. The only cell state that is open throughout is the parked state L^(p-1) L_0 E. Hence [proved, combinatorial] a vacancy crosses a cell only through the parked state, and a vacancy inside the L region confines the walker to its right. If t were comparable to or larger than J, the vacancy kinetic energy would park every walker and destroy the clock vacuum on which the certificate rests [heuristic]. So the mobile member lives at t << J, where vacancies sit in the inter-cell gaps (after a marker) and tunnel through cells.

### The slow-vacancy fluid

perturb.py computes, at zero port, K = infinity and bits zero (exact for every bit string by bit blindness), the Rayleigh-Schroedinger coefficients of a cell flanked by two frozen gaps, and the lowest-order tunnelling amplitude through a cell. At p = 12 (C84's member):

| quantity | value |
| :-- | :-- |
| first entry cost into a cell, 2cos(pi/14) - 2cos(pi/13) | 0.0080 J |
| binding of a vacancy to its gap | 123.1 t^2/J (+ 0.0071 abs(w)) |
| extra cost of a doubly occupied gap U | 123.1 t^2/J |
| gap-gap interaction across a cell, w = 0 | X = -14.2 t^4/J^3 (attractive) |
| gap-gap interaction, linear in w | +7.7e-4 abs(w) t^2/J^2 for w < 0 |
| tunnelling through one cell | t_eff = 1.92e9 t^13/J^12 |

The attraction grows with p (-1.55, -2.02, -3.62, -7.24, -10.4, -14.2 at p = 2, 3, 5, 8, 10, 12). It has a clear mechanism: a vacancy entering from the left pushes the walker to the right, which raises the parked amplitude that the right-hand vacancy needs. Since t_eff is of order t^13, an attraction of order t^4 dwarfs the hopping, so with w = 0 the vacancies would bind into clusters of consecutive occupied gaps, an incompressible region with no contacts, and the extra boson would be lost [heuristic, standard for hard-core bosons with nearest-neighbour attraction much larger than hopping]. With abs(w) > 1.84e4 t^2/J the interaction is repulsive at leading order. The fluid is then hard-core (U >> t_eff) with nearest-neighbour repulsion X >> t_eff, which at filling 0 < nu < 1/2 per gap is a gas of hard rods of length two and a Luttinger liquid with K = (1 - nu)^2 [heuristic]. Interactions between gaps two apart need a port process across an empty gap and are of order lambda^2 alpha^4 abs(w), which the windows keep below t_eff.

### The member

[definition] The material schedule is C84's MEMBER.json unchanged (p = 12, b = 13, n = 11, D = 1000, K = 2^2048, lambda = 2^-480, g = 2^-216, t_d = beta 2^-55, mu_sel at the 13/12 midpoint), J = 1, on the F = 0 sector. MEMBER.json records two windows, each satisfying every [heuristic] constraint of member12.py (perturbative, repulsive, clock distortion below the certificate tolerance, t_eff above the second-neighbour and payload scales, hard rods).

| window | t/J | -w/J | log2 X | log2 t_eff | log2 payload W | relative clock distortion |
| :-- | --: | --: | --: | --: | --: | --: |
| normalized 2^-35 allowance | 2^-30 | 2^-40 | -110.3 | -359.2 | -1665.5 | 2^-47.1 |
| raw tolerance 2^-163 | 2^-88 | 2^-158 | -344.3 | -1113.2 | -1665.5 | 2^-165.1 |

The second window is tight. There the repulsion exceeds the residual attraction by a factor of about 14 and t_eff exceeds the second-neighbour estimate by about 2^22. The "relative clock distortion" is the size of the vacancy-induced shift of the clock vacuum (abs(w) alpha^2 + 123 t^2), which I use as a proxy for the non-scaling deformation of the certified Bloch matrix; it is not a certified bound.

## Checks

| script | output | what it shows | tag |
| :-- | :-- | :-- | :-- |
| check_parent.py | parent_checks.json | material block equals C88's to 1.4e-14; N, M, Q commute exactly; squeezing lemma on all 121 pairs | proved |
| check_rings.py | ring_checks.json | all 12^L rings, L = 2..6: 0 squeezing violations, 0 F violations, 0 switch leaks, hop leaks only in the E-context (12000 of 7.7e6 transitions at L = 6) | proved |
| check_frozen_identity.py | frozen_identity.json | t = 0 gap block = frozen ring minus cut ports, six closures up to dimension 8764, max entry error 7.1e-15; control without switch fails | proved |
| check_lift.py | lift_checks.json | illegal material weight at K = 10, 100, 1000, 10^4 (two cells, nonzero bits, one vacancy, port) 5.72e-3, 6.10e-5, 6.14e-7, 6.15e-9; transparent control plateaus | numerical |
| check_bitblind.py | bitblind_checks.json | at K = infinity, zero port, t != 0, w != 0, all 16 bit strings of two cells give the same spectrum to 1e-13 (equal cells give a subset) | numerical, exact to rounding |
| obstruction.py | obstruction.json | maximal rectangles, minimal colour cover, crossing analysis | proved |
| perturb.py | perturbative_scales.json | X, U, binding, t_eff coefficients for p = 2..12 | numerical (exact linear algebra) |
| vacancy_scales.py | vacancy_scales.json | ED cross-check of perturb.py: X/t^4 = -3.613 against -3.621 (p = 5, t = 0.0025); tunnelling splitting within 0.1 percent at p = 2 | numerical |
| check_dictionary.py | dictionary_checks.json | three-cell, one-vacancy ring against the frozen ring | numerical |
| member12.py | MEMBER.json | the member and the scale ledger | definition, heuristic |

The dictionary test needs a word of interpretation. In the frozen three-cell ring the port shifts the ground energy by -2.874e-5, and the one-fermion level by -4.862e-12. In the mobile ring with one vacancy, whose band is mobile (t_eff = 0.5 t^3 at p = 2), the expected factor is 2/3 because one of the three bonds is occupied at any time. The port shift comes out at 0.6859, 0.6763, 0.6716 of the frozen value at t = 0.1, 0.07, 0.05, deviating from 2/3 as t^2 (coefficient about 2). The payload-specific level comes out at 0.974, 0.827, 0.754 at the same t, also approaching 2/3 as t^2 but with a coefficient near 35. A contact w = -0.02 at t = 0.07 changes these to 0.6727 and 0.823, so w itself is harmless here. The payload hopping at these small cells is itself a strongly cancelled quantity (about 1e-6 of the naive lambda^2 g^2 alpha^4 estimate), which seems to amplify the relative effect of the vacancy-induced clock distortion. This is the reason the windows above put the relative distortion far below the tolerance rather than merely below it.

## Obstructions and what does not work

[proved] Material-independent placement needs two colours. Suppose every bond of every legal word must admit a vacancy, so that the vacancy fluid is independent of the material word (the Ogata-Shiba structure on which C88's dictionary rests). With k colours this asks for k rectangles A_c x B_c inside the legal relation whose union covers every ordered pair USED in legal words. For this grammar USED has 56 pairs and no single rectangle covers it (R_0 R_0 and L_0 L_0 are used, but R_0 L_0 is illegal), while C88's two rectangles do. So within that class d >= 11 + 2 = 13, and C88 is optimal. Combined with C86's restricted eleven-letter bound this gives d >= 13 for every single-token product encoding with material-independent vacancies.

[proved] Only one single-vacancy grammar works. The legal relation has exactly three maximal rectangles. State 0 (A = {L_a, E, M_s}) is the design above. State 1 (A = {R_a, tokens}, B = {R_a, E}) never admits the bond after a marker, so a vacancy can never leave its cell and the chain stays frozen. The third (A = everything, B = {R_0, E}) admits a vacancy only before R_0 or E, so it can cross a cell only when every site is R_0, that is, only for the zero bit string. Its kinetic energy is then bit dependent at order t, which splits the 2^n code at the vacancy scale.

[heuristic] d = 11 by reusing a letter as the vacancy, or by marker hopping. A letter that moves through the L region by swaps must be allowed on both sides of an L, which leaves L_a, E and M_s. L_a carries a bit and M_s a token. An E moving through a cell splits it, and the split fragment is a legal end-state cell only if the letter before the moving E is L_0. Moving E across L_1 changes the conserved pin count F, so the motion is bit dependent. Marker hopping without vacancies (idea (iii)) moves a data letter between cells and has the same defect, besides changing p and with it the number of flavours per cell, which the certificate fixes. I found no repair that keeps the code degenerate.

[proved, as a statement about this design] Two features of the single-colour chain are unavoidable rather than artefacts of my choices. Vacancy placement depends on the walker, since only the parked state is open, and the walker therefore mediates a gap-gap interaction at order t^4 against a hopping of order t^(p+1). Any single-colour mobile chain built on this grammar inherits a sign question for that interaction. Here it is settled by w, but it is a new place where the construction can fail.

## Reproducing

From this folder, `python3 -B run_all.py`, or the scripts one at a time in the order of the Checks table. check_parent.py needs the repository copy of codex/work_C88/bond13.py (read-only). Times on a laptop are a few seconds for the exhaustive and exact checks, about two minutes for perturb.py, and of order a quarter of an hour each for check_lift.py and check_dictionary.py. The `*_run.log` files are the raw logs of the last runs.
