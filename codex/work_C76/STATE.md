# C76 working memory

Completed 2026-10-28. Scope is the small-d search outside route A, with all writes confined to this folder.

Verdict: no admissible counterexample established at d=3 through 12. Best c-(d-1)=0, attained by the fundamental SU(d) permutation and finite-density one-hole q-flavor t-J family. Largest surveyed c/d is 11/12 at d=12. No general upper-bound theorem was proved.

The hard-core question is resolved for the proposed mechanism: q complete density liquids have c=q, with spin-charge separation giving 1+(q-1). Hole density is constrained. Pure hopping preserves the ordered color word, even with species-dependent hopping. Equal hopping has one propagating charge liquid plus frozen color degeneracy. An arbitrary color-conserving nearest-neighbor chain could in principle have additional neutral critical degrees of freedom; this survey does not exclude them.

Important examples and caveats:
- Interacting two-orbital Fock-space swap: d=4,c=3, beating C1's quadratic bound c<=2 but saturating d-1. Explicit parity-correct fermionic formula verified. Three-orbital version has d=8,c=7.
- Doped spin-S integrable family uses d=4S+1, not a spin multiplet plus a single scalar hole. Low-doping charge/WZW/minimal-model count is c=5-6/(2S+1). d=5 gives c=3. Generic-doping magnetic thermal coefficients should not each be identified with fixed independent minimal-model central charges.
- Mixed rational OSp(n|2m) standard bond is non-Hermitian in the ordinary site metric. E²=(n-2m)E and ||h-h†||_F=8sqrt(nm)/|2-n+2m|. Nonunitary/twisted effective charges are excluded.
- Sp(2r) with positive rank-one K_J uses h=P+K_J/(r+1). The sign was audited by [H,Q3] on L=4. A minus sign with the positive K convention fails that check. Published multivelocity interpretation gives c=r, and can give a misleading one-velocity coefficient exceeding d-1.
- Ordinary q-component Hubbard has d=2^q. The literature projected degenerate model allowing at most two particles has d=1+q+C(q,2), with q modes.

Files:
- RESULT.md: survey, explicit representative bonds, citations, verdict, finite-size table.
- PROOF.md: density count and its limitations, order conservation, t-J identities, fermionic swap, OSp matrix audit, finite-size methodology.
- ed_checks.py and ed_results.json: SU(3), SU(4), t-J, and frozen-order ED.
- extra_checks.py and extra_results.json: Hubbard spectra and OSp bonds.
- audit.py and audit_results.json: fermionic swap identity, half-chain estimators, free correlation-matrix check.
- symplectic_check.py and symplectic_results.json: Sp(4) sign and Casimir checks.

Reproduction commands are in RESULT.md. Main ED used SU(3) L=6,9,12; SU(4) L=4,8; t-J L=6,9,12. Degenerate t-J L=9 entropy is explicitly excluded. Half-chain c values are 1.982403, 3.089079, and 1.981850 respectively. No new thermodynamic c is inferred from these small rings. Hubbard entropies must be blocked by the two orbital charges, not by four local-configuration counts. Free Hubbard L=8 agrees with the independent correlation-matrix entropy to 2.23e-15. No project-level memory or other folders were modified.
