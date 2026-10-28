# C32 local state

[Proved record] Final handoff, 2026-10-28. Read RESULT.md for status,
PROOFS.md for all derivations, and PROGRESS.md for the audit history.
The parent project STATE was not modified because the user restricted
all writes to this folder.

[Proved] The positive modified folded bus has 28 core letters. Its core
length is k=n+2. The C8 wrapper costs d0=96*28+48=2736 and selects physical
lengths n+3,n+4 with the same n-mode payload. Exact U(1) and parity survive.
All finite-cell equalization and crossing constants are given in PROOFS.

[Proved] The leading hopping is
T=-g^2 [2/(n+2)^2] D^(-2n-4) Gamma_n + higher terms,
when the neutral express phase is pi/(2n). Gamma_n is positive and given
by a finite Lyapunov equation. PROOFS gives a finite express-energy bound
and a finite small-g/small-loop schedule preserving all n band crossings.

[Assumed, unresolved] No all-n distinct-center/Sidon proof was obtained.
The stable infrared hierarchy and c=n+1 remain conditional. Positivity
of the hopping is insufficient, and the independent-momentum C20 volume
argument cannot be applied without a nonidentity lemma for this family.

[Numerical] All four two-/three-cell, n=1,2 checks completed in the
fast-wrapper model with active lengths, full legal core excitations,
and total charges 0 and 1. Largest sector 91872. They show lambda^4 error
against the independently computed second-order prediction. This is not
full d0=2736 microscopic ED and does not test every payload charge sector.

[Proved record] Reproduction commands and all numerical limitations are
in RESULT.md. ED resumes saved sectors and couplings. The six Python
sources use only numpy/scipy and Python standard-library support.
