# C86 alphabet bounds

[Result] No alphabet reduction was certified. The new eleven-letter
bound below is restricted, as are C84's bounds.

| Assumptions | Bound | Status |
| :-- | --: | :-- |
| Retained on-site charge roles and neutral two-channel port equations | 7 | Inherited from C84 |
| Also retain independent binary L/R backgrounds | 9 | Inherited from C84 |
| Faithful shifted-bit express-fold encoding | 11 | Inherited from C84 |
| Single-token product cells, matched common species, local point ports, C84 leading nonresonance | 11 | Proved here, without assuming the express fold |

[Proved, restricted] The last count is 4+3+2r. Four backgrounds store
the independent bits on both sides of the head. Three markers carry
the neutral input and opposite-charge outputs. Every transport type
costs two charge-distinct carrier letters. With a<=10 there is only
one such type. Its single-token nearest-neighbor graph is a path and
its local intercell source is supported at one vertex per cell.
All scalar path phases can be gauged away. The resulting Bloch
spectrum is reciprocal and cannot pass the pair-sum certificate for
n>=2. See PROOF Section 7 for the representation assumptions.

[Proved] Granting extended source vectors does not rescue a reflected
common carrier. Reflection gives F=GJ and h(-k)=J h(k) J, including
complex scalar graphs. The hopping is indefinite when both parity
sectors have independent filters, and the pair-sum margin vanishes
whenever there are at least four distinct roots. One common phase
only shifts the reciprocal center.

[Proved, restricted] Encoding conversion in a charged carrier bit
still requires two orthogonal carrier charge states. If conversion
instead interchanges a charged carrier and the vacuum, a wholly
quadratic enlarged model requires at least n independent intercell
quadratic channels to have n simple Dirac pairs. An interacting
data-dependent port is outside both statements.

[Proved finite exclusions] Twenty asymmetric-terminal single rails
and 98 extended-source instances have certified negative principal
determinants. Their rejection is robust below the rational sign gaps
stored in obstruction_certificate.json. These tests do not establish
that every asymmetric extended-source design fails.

[Open] The remaining ways below eleven must leave at least one of
the common-carrier, product-token, local point-source, or linear
conversion assumptions. Possible changes include unequal carriers
with a controlled onsite cancellation, a new data-dependent port,
or contextual storage of transport information. No complete bond for
any of these was obtained. There is no universal a>=11 theorem here.
