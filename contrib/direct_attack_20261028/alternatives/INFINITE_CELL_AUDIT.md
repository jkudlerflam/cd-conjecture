# Independent audit of operator-valued auxiliary matching

**Verdict: the revised abstract vacuum and thermal matching arguments pass this independent analytical review under their stated hypotheses.** The reviewed proof is [INFINITE_CELL_PROGRESS.md](../infrared/INFINITE_CELL_PROGRESS.md), including its thermal Section 10, continuous-time translation argument, strengthened source norm, and explicit application restrictions. No microscopic central charge or completed C44/C45 application is certified by this verdict.

The review was independent of the author's derivation. It included concrete objections and repairs exchanged with the author and root coordinator. It is an AI mathematical review, not formal verification or human refereeing. The separate finite-cell determinant argument has its own [audit](IR_AUDIT.md).

## Representation and coefficient norm

The embedding

\[
\mathcal K\longrightarrow\mathcal F_n\widehat\otimes
(\mathbb C\Omega\oplus Q\mathcal K)
\]

is an even isometry, taking the full code to logical states with auxiliary vacuum and excited states to logical vacuum with their original auxiliary state. The auxiliary on-site operator has a unique vacuum and preserves the original positive excited energies. No auxiliary eigenbasis or finite excited dimension is required.

Expanding only the finite logical matrix blocks gives the stated bound \((8H^2)^{nr}\|O_X\|\) on a support of \(r\) cells. Auxiliary coefficients remain bounded operators. This removes the dedicated-fermion-mode counting loss; it does not secretly assert a tensor factorization of the original cell. The strong input norm sums all coefficient operators and species per support, so it also pays the rooted-source sums that a pinned-field seminorm alone would miss.

The baseline condition is two-sided annihilation of the product auxiliary vacuum on every support. Logical matrix-unit independence then implies the same property for each auxiliary coefficient. This is stronger than zero vacuum expectation and is sufficient for the rightmost-vertex cancellation with arbitrary external logical fields.

## Block estimate and exact factorization

The proof uses the semigroup Duhamel formula. Every \(e^{-sa}\), \(s\ge0\), is a contraction. It never estimates an unrestricted interaction-picture expression \(e^{sa}Ve^{-sa}\), which could be unbounded.

For fixed active interval labels, summing the label multiplicities and the ordered time simplexes gives the product \(\prod_X(e^{\tau v_X}-1)\). A passive excited site undergoes a full block of free evolution and pays the positive auxiliary-gap factor. These estimates use operator norms and remain independent of excited-space dimension.

I checked the proposed exact polymer factorization as a tensor-network statement. At an inactive vacuum node, write \(P=|\Omega\rangle\langle\Omega|\) and cut that wire. Different resulting connected contractions factor. The author's occupied square-grid components may merge several otherwise independent contractions, but cannot split a genuine contraction: an active interval fills its full spatial support in its block, and successive occupied blocks sharing a cell are vertically connected. This establishes factorization between distinct coarse components.

No additional binary projector decoration is needed inside active blocks: their exact resummed operators remain unsplit. The inactive excited/vacuum choices are exactly the recorded passive set. My initial suggestion to add an active-node projector factor was unnecessary for this representation.

The graded tensor product is essential. Each input vertex has even combined parity. A completed vacuum matrix element or thermal trace vanishes for odd total auxiliary parity; hence every surviving component has even external Grassmann parity. Completed activities therefore commute, although open auxiliary operators need not. This is the required algebraic premise for the hard-core polymer logarithm.

## Arbitrarily long interval supports and convergence

The revised proof sums active interval decorations with their activities, rather than counting all interval subsets indiscriminately. With \(\epsilon=e^{-L}\),

\[
\eta=\sup_j\sum_{X\ni j}\epsilon^{-|X|}
(e^{\tau v_X}-1)<1,
\]

the decoration sum on a specified \(N\)-node support is at most \(\epsilon^N e^{\eta N}\). The extra thermal heat-trace factor changes this only by \(e^N\). This is sufficient for arbitrary exponentially summable support lengths and avoids the invalid limit of a finite-range \(2^{RN}\) count.

The depth-first-walk animal count \(16^N\), nearest-neighbor incompatibility neighborhood \(5N\), and stated choices of \(L,\tau,\kappa\) leave ample room in the Kotecký–Preiss criterion. I independently opened [Fernández–Procacci, Section 2, (2.15)](https://webspace.science.uu.nl/~ferna107/papers/ferpro07.pdf), which states this criterion and its pinned connected expansion bound. Applying its scalar majorant to commuting coefficient-algebra activities uses the submultiplicative norm and exact factorization just checked. No positivity of the quantum activities is needed.

The output bound and its Cauchy interpolation give a correction proportional to the residual interaction. The auxiliary interaction can be larger than the logical bandwidth; no logical propagator or inverse logical gap appears in this integration.

## Continuous time and pinned kernels

The finite logical operators are normal ordered. Point splitting fits all fields of each vertex into one interval of uniformly bounded width, with creation fields later than annihilation fields. This prevents an unbounded-degree loss when spatial support sizes grow. The coefficient norm and polymer majorant control the limiting operation.

The initial cluster bound is an integrated coefficient density. The revised proof explicitly establishes continuous time translation invariance before using that density as a pinned-coordinate bound. Arbitrary block origins produce outer partial blocks of width at most \(\tau\). An inactive excited projector there directly meets a vacuum boundary and gives zero; active-set estimates improve when a block is shortened. Thus the majorant is uniform in the block origin. Exact Duhamel partition independence and exponentially suppressed boundary-reaching clusters give the same infinite-slab kernels for every origin.

For a translation-invariant coefficient kernel, integrating all coordinates and dividing by spacetime volume equals the integral over relative coordinates with one coordinate pinned. Summing anchor species is stronger than taking their supremum. The diameter weight is invariant under this change of variables. This justifies the proposed pinned norm in the homogeneous setting. A general inhomogeneous theorem would require the separately marked local estimate that the author does not claim here.

## Thermal trace and winding/source control

The thermal extension requires an actual heat trace, not only a spectral gap. Its assumption and choice are

\[
\operatorname{Tr}_E e^{-\tau a}\le e,
\qquad \beta\ge2\tau.
\]

Using \(\lfloor\beta/\tau\rfloor\) blocks gives widths in \([\tau,2\tau]\), so the argument covers arbitrary sufficiently large \(\beta\), with the enlarged active-set bound and reduced small-coupling threshold stated by the author.

I independently checked the half-semigroup extraction. At a passive site, no insertion in that block touches the site, so one half of its free factors can be collected with its excited projector without crossing a noncommuting insertion. It pays \(e^{-g\tau_t/2}\). The other half of every global free segment is retained. Generalized Schatten Hölder with exponents \(\beta/s_i\) bounds the completed component trace by

\[
\operatorname{Tr}e^{-\beta A_S/2}\prod_i\|B_i\|
\le e^{|S|}\prod_i\|B_i\|.
\]

The component's spatial support satisfies \(|S|\le N\), giving the advertised \(e^N\) loss. Inactive rank-one projectors split thermal wires as well, including components using the same site in separated time intervals. Graded cyclicity and the antiperiodic identification of odd logical sources preserve the thermal functional; absolute coefficient kernels have the required continuous translation invariance.

The thermal baseline estimate controls nonconstant kernels as well as pressure. A nonzero baseline-only component cannot leave an entirely empty time layer. Cutting at such a layer would produce a vacuum amplitude, killed by its first baseline insertion. Every surviving component therefore has at least \(\lfloor\beta/\tau\rfloor\) nodes. The additional exponential cluster reserve suppresses its complete weighted coefficient norm, including arbitrary external logical degrees. This supplies the stated exponentially small baseline term without inferring a source bound from a pressure bound.

## Invalid states and the available C44/C46 input

For the abstract embedding the invalid heat trace is exactly

\[
z_{\rm inv}(\beta)=(2^n-1)\operatorname{Tr}_Q e^{-\beta A}.
\]

Embedded interactions touching an invalid site vanish. Cutting them costs at most the stated incidence norm times the number of invalid sites. The resulting full-trace comparison with rate \(g-w\) is valid when \(w<g\). This replaces the earlier vacuum-only concern about sectors orthogonal to the global auxiliary vacuum, while leaving the physical-length ensemble transfer separate.

[C46 §2](../../../codex/work_C46/PROOF.md) gives a concrete heat-trace input **on its specified legal/well-formed rooted space**. With

\[
\tau_0=\frac{\log(2d_0)}\gamma,
\quad C_0=D_{\rm fin}+b d_0^{b-1}+2d_0^{b+2},
\]

its counting argument at zero intercell norm yields

\[
\operatorname{Tr}_Q e^{-\beta A}\le C_0e^{-\beta\Delta}
\quad(\beta\ge\tau_0).
\]

Thus, taking \(g=\Delta\), a sufficient finite choice for the thermal theorem is

\[
\tau\ge\max\left\{\frac{2L}{g},\tau_0,
\frac1g\log\frac{C_0}{e-1}\right\}.
\]

This demonstrates that the heat-trace hypothesis has an explicit candidate input. It does not extend that input to every malformed selector string without the original invariant-sector qualifications. In particular, C45's unbounded positive multi-token readout terms must remain in their invariant onsite blocks or be treated by the existing exclusions, not inserted into the bounded interaction norm.

## What this verdict does not establish

The abstract theorem still needs its complete hypotheses instantiated for the chosen microscopic invariant space, with all constants fixed before the final coupling. The physical total-length constraint, cyclic quotient, minimizing marker density, grading/charge dictionary, and competing microscopic sectors need their own transfer arguments.

The resulting action has general retarded, high-degree logical interactions. Its multichannel infrared flow and a Virasoro central-charge or current/stress scaling limit remain separate mathematical obligations. An exponentially small auxiliary thermal baseline is useful for that analysis; it does not itself calculate the \(T^2\) coefficient of the interacting logical pressure.

Reviewed 2026-10-28 by a separate GPT-6-based Codex agent. The revised vacuum and thermal steps were checked after substantive corrections and explicit qualifications were added. No independent human review is claimed.
