# Parameter and term audit

Date: 2026-11-17. The user's preference is a genuinely simpler Hamiltonian:
fewer physical terms, moderate coefficient ratios and less delicate tuning.
This audit does not count repackaging a sum as deleting its terms.

## Independent controls versus physical coefficients

At fixed n, and fixing the material express/junction/clock normalization
used by the parent, its continuous controls can be taken as

    q, w, v, eta, D, g, lambda, f, K.                    (1)

There are nine. Delta=4+E_n(q), mu_sel=-(n+2)E_n(q), the terminal
coefficient q^2 and the phase pi/(2n) are prescribed, not additional
independent knobs. In the band schedule u replaces q through q=(wu)^2
and nu replaces f; that change of variables does not add controls.
The mobile model adds J/t and mu_occ/t after removing one common energy
scale, giving eleven continuous controls in this written family.
At fixed occupied number mu_occ N is a scalar; in the full density
selection problem it is a physical control. K has no effect inside the
exact legal block, but is physical on the full Hilbert space.

Nine controls do not mean nine matrix entries. The literal material
bond has 64 undirected off-diagonal transitions (128 directed entries):
4 vacuum moves, 8 data moves, 16 express moves, 28 marker-conditioned
junctions, 4 conversions and 4 ports. The mobile lift adds 20 vacancy
swaps. Symmetry ties many coefficients together, but does not erase the
corresponding processes. Onsite and grammar diagonal terms remain.

With t=1, the physically displayed nonzero scales include

    J, Jq, Jq^2, Jw, Jwv, Jg,
    Jlambda, Jlambda sqrt(eta),
    JDelta, J(Delta+D), Jmu_sel, Jf, JK, mu_occ.           (2)

Signs and the prescribed complex phases matter in addition to their
magnitudes. The vacuum entries 1,-q,q^2 are tied by a positive rank-one
projector; they are not three independent parameters.

Only a constant on the full 22x22 bond and a full-space coboundary
A tensor I_22-I_22 tensor A are automatically irrelevant on a ring.
A material identity before lifting becomes an occupied-contact
interaction. A material coboundary becomes
sum_j (A_j-A_(j-1)) a_j after contact weighting, and is generally
physical. In particular, redistributing the material onsite halves
without compensating the vacancy bonds changes this mobile model.
The marker contact energy cannot be removed as a constant even at
fixed marker number.

## What is structural, and what is a proof allowance

| Choice | What it does | What is actually established |
| --- | --- | --- |
| 0<q<=1/4 | Biased clock and geometrically decreasing cell energy | q=1/4 already proves the all-length selector and the new clock estimate. A much smaller q is not required for those two tasks. |
| Terminal q^2 projector | Gives E_k>0 and permits finite-length selection | Removing it destroys the selector. Its mobile clock leakage is real, not a missing large penalty. |
| mu_sel=-(n+2)E_n | Selects n among every cell length | Its exponentially small scale with n is structural for this selector, although its exact value can have a stability interval. It is not a removable global chemical potential in the mobile model. |
| Distinct w and wv, eta>0, relative phase i | Two independent hopping Gram matrices | The present arbitrary-n proof uses this distinction. Setting v=1 or eliminating the second bus invalidates that argument. A simpler alternative is not established here. |
| Small u, v, nu | Separates slopes, centers and low-order momentum sums | These are sufficient hierarchical choices in the existing band proof. No necessity theorem or moderate-ratio arbitrary-n certificate is known. |
| Large D | Justifies the express-rail Gram limit | Needed by the written limiting argument, not by grammar, the selector or the new clock theorem. Whether a fixed moderate D gives the desired bands for all n remains unproved. |
| Small nonzero g | Makes the Gaussian lambda^2 g^2 term dominant | g=0 removes that payload hopping mechanism. The extreme smallness is a matching allowance, not a demonstrated physical threshold. |
| Small nonzero lambda | Couples cells while controlling defects and higher orders | The frozen sufficient bound B<=E_(n+1)/5 is conservative. Nonzero coupling is essential for the proposed dispersing payload. The clock-pin theorem itself adds no lambda restriction. |
| Positive tiny f | Removes leading-band resonances | The current all-n argument needs a nonzero readout. The new clock estimate only asks nf<=1/64, a vastly weaker allowance. |
| J/t small | Fast vacancy regime | The old charged-rail estimate uses a full-bond norm involving D. The new clock proof uses only c+q^2+|mu_sel|/2+nf and permits J/t=1/80 independently of D and n. Full fast-vacancy matching may still impose stronger conditions. |
| Finite K | Penalizes forbidden arrows | K=1 suffices in the completed frozen comparison. Exact invariance removes any need to suppress legal leakage by taking K huge. A mobile illegal-sector comparison is still missing. |

The new theorem removes a **proof-induced** D and cell-length factor
from the vacuum-clock estimate. It does not remove the hierarchy of the
whole band construction or certify J/t=1/80 as a high-c microscopic point.
The existing architecture remains an elaborate existence construction.

## Two literal simplifications, with sharply different scope

### Set Delta=4 exactly

This removes the tiny calibration +E_n from one coefficient. It deletes
no physical operator. Grammar is unchanged, the vacuum selector is
unchanged, and every excited frozen cell still has energy before its
marker shift at least 4-2-|g|>=3/2>E_k. Both new clock theorems continue
to hold; they use no Delta-dependent error bound.

In the pair resolvent the scalar rail diagonal becomes a=4-E_n instead
of 4. It remains positive. The leading Gram expression is still positive,
with limiting entries

    C_a(0)_(ij)=(2i+2j+2)! /
       [(2i+1)!(2j+1)!(2a)^(2i+2j+3)].

The two-charge Gaussian cancellation and geometric valuations use a>0,
not exactly a=4. Thus the qualitative finite-order band construction
can be run with this simpler coefficient. No numerical band margins,
full mobile matching or CFT certificate are promoted by that observation.
The local and finite-pair checks include this literal variant.

### Set D=0

This genuinely removes the additional onsite projector on X and excited
M letters. With Delta=4, grammar and the all-length zero-port selector
survive by the same bound 3/2>E_k. The vacuum-clock Hamiltonian and its
proof are completely unchanged, since this projector vanishes on the
vacuum space. This is independently checked through n=2.

However, D=0 removes the large-D traversal regime used to prove the
parent's arbitrary-n Gram formula. It is a simpler **selector and clock
subsystem**, not a certified replacement for the high-c construction.
Deleting g, lambda, or the second bus would likewise preserve a selector
while losing the stated payload mechanism; that is not progress toward
a simpler counterexample without a replacement argument.

The parent files are left intact. The reproducible checks construct both
variants explicitly, and distinguish the full mechanism from these
limited simplifications.
