# Proposed coordinator summary

**[Proved; internal numerical checks, not independently reviewed]** The
modified markerless d=9 family now has a connected local normal form through
fourth port order with an explicit exponentially weighted interaction-norm
remainder independent of the number of cells. The unprojected residual is
O(lambda^5); only its code projection improves to O(lambda^6). This avoids
the invalid global Kato argument for an extensively broad code band, but
does not assert an exact low-energy block or a CFT.

The calculation retains the first conversion-generated quartic at
lambda^2 g^4, the complete connected three-cell lambda^4 port operator,
readout-port terms, and the readout's own clock distortion. An exact
all-n current formula proves C_ij>0, so the L1 field generates attractive
density interactions -2 h_L^2 C_ij n_i n_j. Literal n=1,2 Fock checks,
an independent wave-operator recursion and finite-lambda diagonalization
verify the coefficients and normalization.

**A separate exact theorem closes the nonzero-port frozen ground-sector
comparison on compatible rings for this modified d=9 bond.** Conserved
two-site boundaries partition all sectors into complete one-head cells
or inert headless cells; zero-data cells and altered head counts are
included. With gamma=E_(n+1), B=sqrt(2)|lambda| and gamma>5B, each
nonselected cell costs at least gamma-5B relative to the exact dressed
selected ground energy. K=0 is sufficient. Choosing
lambda<=(1-q^2)^2 q^(2n+4)/(10sqrt(2)) is compatible with subsequently
smaller matching tolerances and the actual nonnegative L1 readout.

The interaction/remainder budgets can be made small relative to the
smallest physical velocity and the pair-momentum margins with fixed
nonzero g,lambda. The finite-order forward Gaussian matrix stays positive
and close to free under the stated bound; forward scattering is retained.
An all-scale integration/RG theorem controlling residual auxiliary
coupling and dressed cutoffs is still missing. No unbounded-c microscopic
CFT theorem follows. The coefficient hierarchy remains severe and is
not improved by these deliberately conservative bounds. The original
mobile d=10 model and incompatible frozen rings remain separate issues.

All bounded checks passed on one worker in about 0.22 seconds; no new
root scan, production calculation or subagent was launched. Source:
[README.md](README.md).
