# C86 result

[Result] No ten-letter design with a certified member n>=a was found.
The fallback result is a stronger restricted lower bound and analytic
obstructions to the proposed replacement mechanisms. Eleven remains
the smallest certified alphabet among the supplied constructions.
No universal ten-letter impossibility theorem is claimed.

[Proved] Reflection does not supply the desired positive hopping in
the retained common-carrier reduction. Its hopping is F=GJ, where
G is a positive semidefinite Gram matrix and J reflects the flavors.
For independent filters, Herm(F) has one positive direction per even
flavor and one negative direction per odd flavor. A common port phase
cannot make it definite. The result holds for complex reflected carrier
graphs and extended endpoint sources, not just real paths.

[Proved] The same reflection also forces

    h(-k)=J h(k) J.

With four or more distinct Fermi points, two distinct unordered pairs
have the same sum. C84's pair-sum margin is therefore zero. Internal
flux compatible with this reflection does not help. Equal scalar
sources at both ends have the same reciprocal-spectrum problem.

[Proved, restricted] At least eleven letters are required in the
single-token product class with four binary L/R backgrounds, frozen
on-site markers, matched common carrier species, payload-independent
pair ports, and nearest-neighbor token motion. With at most ten
letters there is only one transport type per species. Locality makes
its scalar graph a path and its sources point vectors, which forces
reciprocal spectra. This extends C84's faithful-fold lower bound to
changed transport mechanisms in that class. The assumptions and
loopholes are explicit in [PROOF.md](PROOF.md), Section 7.

[Proved, restricted] A single fixed-charge species cannot perform a
number-preserving conversion that flips one charged data bit. Letting
conversion create or destroy the carrier can instead give a quadratic
hybridization model, but an intercell hopping block of rank r then
allows at most 2r isolated Fermi points. Obtaining n Dirac pairs by
that mechanism requires r>=n. This bound does not apply to general
interacting, data-dependent ports.

[Proved, computer-assisted] C84's integer resolvent enclosure method
rejects 118 proposed scalar instances before positive-metric whitening.

| Family | Exact checks | Outcome |
| :-- | --: | :-- |
| Single rail, n=2,3,6,9,10, four terminal masses | 20 | Negative principal determinant |
| Two-vertex endpoint sources, n=9,10, seven complex amplitudes per end | 98 | Negative principal determinant |

The certificates include rational sign-stability bounds and reject
every common hopping phase. The arbitrary-mass two-vertex obstruction
is also proved analytically and checked at twelve rational points.
Four direct Fraction inverses independently check selected integer
enclosures. Numerical heat-kernel checks verify the reflection and
Gram identities on complex graphs of sizes 4,5,6.

[Rejected check, corrected] The first extended-source checker assumed
positive entries and stopped with an assertion failure. The corrected
checker uses signed interval determinant bounds. The failure is
retained in [rejected_checks.json](rejected_checks.json).

[Scope] No candidate reached the root-inertia or finite-K transport
stages, so there is no new root certificate, microscopic coupling
schedule, or central-charge prediction. The asymmetric extended-source
scan is finite, not an all-parameter exclusion. Unequal carriers with
controlled onsite terms, contextual encoding, and nonquadratic
data-dependent ports remain open. No independent referee was used.

Reproduce from this folder with Python, NumPy and SciPy.

    python3 certify_obstructions.py
    python3 independent_checks.py
    python3 consistency_checks.py
    python3 explore.py

[Completed] All scripts, outputs, and notes are in work_C86. No network
command, subagent, cluster job, or modification outside this folder was
used. [PROOF.md](PROOF.md) contains the derivations,
[ALPHABET.md](ALPHABET.md) separates the lower-bound classes, and
[STATE.md](STATE.md) carries the handoff.
