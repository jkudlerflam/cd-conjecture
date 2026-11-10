# Archived physical-chain lead audit

Read-only audit on 2026-11-10. `audit_archived_leads.py` recomputes the current
central-three-cut and fractional-window estimates from all physical-chain
DMRG JSONL records in C184's `typhon_out/out`, `local_out`, and `c152_j2_copy`.
It retains 194 distinct source/parameter groups with R≤10 and at least two
lengths; it does not rank single-size profile estimates as candidates.
`archived_lead_audit.json` retains coefficients, exact source lines, sizes,
chi, discarded weights, profile fits and all three entropy-window slopes.

R is the ratio of the largest to smallest nonzero displayed Pauli/field
coefficient. It is coordinate dependent, so it is a screening preference,
not a physical invariant. Existing C184 records often round couplings to
10 decimals and entropies to 6 decimals. Their discarded weight is not a
substitute for a chi-convergence test. Most old runs did not enforce parity,
so comparisons to the new parity-resolved states require caution.

No stronger overlooked preferred-ratio (R≤3) lead was found than the
phi=1.2, Delta=0.8, V=-2 reference now under active continuation.

## Ranked useful checks

For these Peotta-family rows, t=1 and hz=0. Their underlying parameterization
is Jx=(cos(phi)+Delta)/2, Jy=(cos(phi)-Delta)/2, D=sin(phi)/2, Jz=V/4.
The exact stored, rounded coefficient vectors are in the JSON audit.

| phi, Delta, V | R | Sizes; final chi | Fractional slopes | Latest profile estimates | Assessment |
|---|---:|---|---|---|---|
| 1.2, 0.8, -2.0 | 2.656 | 48/96; 128/160 | 1.4215 | 1.316 / 1.527 | Best preferred ratio; now `pref_ref` and rationalized `pref_rat` |
| 1.2, 0.8, -2.8 | 3.199 | 48/96; 128/160 | 1.3351 | 1.220 / 1.369 | Closest untested broader alternative; only one size pair |
| 1.2, 0.6, -1.8 | 4.050 | 48/96; 128/160 | 1.2891 | 1.286 / 1.668 | Broader alternative; fit-window dependence substantial |
| 1.0, 0.7, -2.3 | 7.767 | 48/96; 128/160 | 1.3140 | 1.168 / 1.254 | Adjacent to existing reference; no old 192 row |
| 1.0, 0.7, -2.4 | 7.767 | 48/96/192; 128/160/224 | 1.1333, 1.3108 | 1.179 / 1.200 | Best old third-size broader lead; already `attr_ref` |
| 1.2, 0.6, -2.2 | 4.629 | 48/96; 128/160 | 1.2269 | 1.229 / 1.550 | Lower-priority one-pair lead |

The profile columns use the two current windows, not C184's original profile
convention. For example, the unsmoothed central-cut slope of the first row
is 1.6078, whereas the registered fractional estimator is 1.4215. Ranking
on the unsmoothed number alone would exaggerate this lead.

The remaining preferred-ratio rows have weaker evidence. At phi=1.2,
Delta=0.8,V=-2.6, R=2.970, the fractional 48→96 slope is 1.1137; at V=-1.8
it is 1.1126 but the latest profile already falls below one. These are useful
neighbors for a coupling scan, not stronger alternatives.

## Other archived regions and current coverage

The C152 copy contains preferred-ratio XZ-neighborhood points. At phi=1.45,
s=1.0, V=+0.8,hz=0 (Delta=s sin(phi), R=2.783), fractional slopes fall from
1.2158 (32→64) to 1.1073 (64→128), chi96/128/160. At s=0.9 the corresponding
slopes fall from 1.1606 to 1.0734. Negative V=-0.8 produces large early slopes
but then drops to 0.35–0.41 or saturates; these are not retained as increasing
leads. The broader worst XZ row phi=1.45,s=1.1,V=+0.4 has R=6.062 and falls
1.4502→1.2453→1.1319 through L192, consistent with the already recorded
crossover in REPORT §9.143. This audit does not exclude a small excess above
one at unattained sizes, but it supplies no reason to prioritize those rows.

C185's only saved full microscopic completed points in `typhon_out_bf/oc`
are F8ref L64 and F8dm_p06 L64, with R=133.679 and 135.474. The latter adds
transverse DM=0.3 and has profile estimates 1.192/1.227. They are single-size
points far outside R≤10; C185's block/effective-model numbers are not full
physical-qubit-chain evidence and are not mixed into this shortlist.

C190's named physical families (`models.py`) have Jz=10 or20 and D around
1/160 or1/320 with fields near -20 or-40. Their R is approximately3196 or
12796. Their stronger c≈2 finite-size readings are already described in
REPORT §§9.175–9.178 and fail the recorded L256 chi-convergence gate.
They are valuable bridge references, but are outside this ratio screen.

REPORT §9.159 already records the broader Peotta third-size flags and the
first-order-boundary alternative. It does not elevate the phi=1.2 preferred
point explicitly. The original `d2_order_one_20261110` search instead tested
the rational attractive R=2 point and the phi=1,Delta=0.7,V=-2.4 reference.
The current continuation addresses this coverage gap with `pref_ref` and
`pref_rat`. The rational point previously tested is a different Hamiltonian.

The larger repository's old XZ flags (C163 hi_T_40, dist_T_50 and C155
neighbors) should not be silently revived based on their old c≈1.3 fits:
REPORT §9.126 documents Néel order and pinning-field entropy saturation, and
§9.136 closes hi_K_1 as a gapped period-two state. Those are useful negative
controls for the present ordered-cat concern.

One distinct lower-priority broader-ratio point remains locally unresolved:
C163/C166 KI26 has diagonal J=(-0.5601821693,-0.4603422503,-0.4179566666),
Dy=-0.4825839437 and fields hx=-0.1574668648,hz=0.2047195604, giving R=3.5575.
The old L128 profile of about1.64 sits on a jump in entropy along the path.
The saved independent two-site VUMPS data in
`codex/work_C166/local/vc/KI26_n2_quick.json`, chi8/16/32, give recomputed
finite-entanglement slopes 1.0053 and0.9831 (6 Delta mean(S)/Delta ln(xi)).
Only these low-chi data are present locally; they give no extra evidence
for c>1. This bond breaks the parity used by the current C190 screen, so it
would require a more general solver. KI22's last analogous slope is0.9808
but its displayed coefficient ratio is26.95. Neither should displace the
current higher-priority preferred-ratio test.

This is an audit of saved local outputs, not a new search over all parameters
or a check for uncollected remote jobs. No old conclusion is treated as a
general no-go theorem.
