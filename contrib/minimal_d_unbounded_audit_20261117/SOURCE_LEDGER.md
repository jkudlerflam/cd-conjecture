# Smallest-d source and scope ledger

Date: 2026-11-17. Status: source audit, not a new CFT certificate.
Author and AI assistance: OpenAI Codex, GPT-6 family. Repository sources
were opened; their old large certificates have not been independently
regenerated in this round. Physical states, including markers and holes,
are counted before encoding or blocking.

## Status of the reported ten-state arbitrary-c example

**Not yet located.** The coordinator searched main, the fetched alternate
branch and public results and requested the collaborator's source. This
round independently inspected the concrete local d=10 records below and
searched the public phrases `"central charge" "d=10" "unbounded"`,
`"central charge" "eleven" "Latorre" Sierra`, and
`"c-d conjecture" counterexample dimension ten`. No matching literal
ten-state arbitrary-c nearest-neighbor Hamiltonian was identified.

This is a search status, not a disproof. A supplied source becomes the
primary target. The original [Latorre–Sierra paper](https://arxiv.org/abs/2403.17242)
states the nearest-neighbor local-dimension conjecture; its abstract
does not supply the reported construction. The scope below is based on
the actual local files, not on interpreting a title or central-charge
number as a microscopic state count.

## Literal inventory

| Record | Physical states and d | Range; microscopic translation | n dependence and certified members | Physical ground sector and infrared status |
| --- | --- | --- | --- | --- |
| [Single-rail seed](../../codex/work_C86/RESULT.md), now being audited by the complementary resonance chat | L0,L1,R0,R1,DA,DB,E,MA,MB: **9** | Nearest neighbor; one-site repeated rules in the retained seed | Arbitrary cell lengths can encode bits; no arbitrary-n band or CFT certificate | Dispersion and resonant interaction audit in progress. C86's failed sufficient certificate is not a physical no-go. A one-hole d=10 lift has not been established. |
| [Frozen demonstrator](../routeA_demonstrator/cell.py) | P0,P1 and B(sL,c,sR), with three binary indices: 2+8=**10** | Nearest neighbor; one-site invariant bond | One logical qubit per PBP cell; no arbitrary-flavour family at fixed d=10 in this source | Exact local checks and second-order mechanism; numerical period-three formation and c approximately 1. Full nonzero-port sector/IR transfer is not established by those numerics. |
| [Two-flavour mobile toy](../mobile_cells/ed2.py) | L,R and #(a,b),X(a,b), a,b binary: 2+4+4=**10** | Nearest neighbor; one-site invariant bond | Two payload flavours; replacing the two bits by n bits would give d=2+2^(n+1), not fixed ten | Finite-ring checks and conditional payload/fluid reasoning. Proposed total c=3 is itself unestablished, and is not arbitrary c. |
| [C53 exact completion](../../codex/work_C53/RESULT.md) | Same P2+B8 alphabet: **10** | **Four-site** BPPB interaction; one-site invariant rules | One logical flavour, no arbitrary-n family | Exact logical reduction and noncode gap within the frozen PBP sector. At the free point c=1 is proved there. Unrestricted type-sector selection and interacting thermal c remain separate; not a nearest-neighbor counterexample. |
| [C84 shifted-bit folded bond](../../codex/work_C84/PROOF.md) | L0,L1,R0,R1,DA,DB,XA,XB,E,MA,MB: **11** | Nearest neighbor; one-site invariant | d stays 11. p=n+1 data sites plus one marker encode n bits. Certified **n=11**, physical cell length 13; not an all-n nonresonance result | Zero-port all-length selection at finite K is proved in the source. Finite-g/finite-K **quadratic** band certificate has 22 simple roots and separated pair sums. Nonzero-port dressed defects, volume-uniform matching, interacting stability and conformal-observable transfer remain open. Conditional frozen c=11. |
| [C80 shifted-bit predecessor](../../codex/work_C80/RESULT.md) | C84 plus E1: **12** | Nearest neighbor; one-site invariant | d stays 12; certified n=12 member, not arbitrary-n nonresonance | Similar finite ingredients and infrared/selection qualifications; not preferable to C84 for this task. |
| [Restricted one-vacancy route](../mobile_below_13/README.md) | C84's eleven material letters plus V: **12** | Nearest neighbor; one-site invariant, with a vacancy-context switch | Inherits C84 n=11 frozen member only in the specified sector; no all-n mobile certificate | V may occupy only automaton-state-0 gaps. Slow hopping, clock-dependent access, and mediated interactions are essential. Bare attraction threatens phase separation; a contact correction is proposed. Total c=12 is heuristic and requires more assumptions than C88. It cannot be dismissed merely because unrestricted vacancies need two colours. |
| [C88 pinless mobile bond](../../codex/work_C88/RESULT.md) and [small-d note](../small_d_note/note.tex) | C84's eleven material letters plus V0,V1: **13** | Nearest neighbor; one-site invariant | Fixed d; inherits the regenerated **n=13** material certificate, not an arbitrary-n theorem | Unique legal colouring and finite grammar corrections are explicit. Frozen pin-free selection is established at zero port; mobile wrong-length/pin selection, finite-J matching, interacting stability and conformal transfer remain assumptions. Proposed total c=14. |
| [C85 pinned mobile lift](../../codex/work_C85/RESULT.md) | Eleven material letters plus three vacancy colours: **14** | Nearest neighbor; one-site invariant | Certified n=13 material member | Three colours are minimal only for its pinned pair-grammar repair class. C88 removes the pin and uses two; the d=12 constrained-vacancy route uses different dynamics. |

All c labels mean the intended total nonchiral Virasoro charge of a
unitary compact infrared theory, not an effective encoded dimension.
No local or quadratic certificate alone identifies that microscopic CFT.
The known numerical d=6 one-flavour mobile toy has states L,R,#0,#1,X0,X1
and evidence for c approximately 2; it also supplies no fixed-d
arbitrary-flavour mechanism.

## What the lower bounds actually exclude

[C84](../../codex/work_C84/ALPHABET.md) needs eleven states for its
faithful shifted-bit two-rail representation: four backgrounds, four
charged carriers and three markers. Its literal charge-preserving merges
fail. [C86](../../codex/work_C86/RESULT.md) excludes the retained
positive-Hermitian/nonresonant certificate in specified single-token
product classes. [C87](../../codex/work_C87/RESULT.md) extends several
obstructions to homogeneous exchange and restricted port dressings,
while explicitly leaving contextual/interacting ports open.

None proves that d<=10 cannot have large or unbounded physical c.
The nine-state seed is being studied separately, including whether its
resonances represent actual instability or only failure of a sufficient
proof route. This round does not duplicate that mechanism analysis.

## What failed at larger sizes

There are three different records, which must not be conflated:

1. **C84 bracket resolution:** the inherited loop exponent 66 gave a
   numerical pair margin about 4.46e-10 but root half-width 1e-9. The
   exact pair-enclosure lower bound was negative, about -6.15e-9.
   Exponent 55 passes, with margin greater than 9.11916e-7. This is a
   rejected certificate prescription, not a demonstrated band collision.
   The failed file is [rejected_parameters.json](../../codex/work_C84/rejected_parameters.json).
2. **Fixed express-energy prescription:** [C33](../../codex/work_C33/RESULT.md)
   records an indefinite Hermitian part already at n=3 for D=10,g=.05,
   and failure of D=100 by n=12 at g^2 order. Its increasing-D(n) Gram
   bound proves positivity for every fixed n. Fixed D=1000 positivity
   certificates reached n=32; a larger computation timed out. These
   facts neither prove all-n centre separation nor an architectural limit.
3. **Unproved all-n nonidentity:** [C37](../../codex/work_C37/RESULT.md)
   proves pair-product noncollision through n=23 for the larger original
   alphabet using exact finite-field certificates. It explicitly leaves
   existence of such certificates at every n open. Principal-block
   positivity transfers to C84; a pair-product certificate does not
   automatically transfer after deleting an active conversion.

The concrete smaller benchmark is therefore C84's eleven-state literal
bond. Its inactive D_n endpoint and X partner must remain in every
propagator. The 2^n code and the finite-K split band must both be retained.
The current finite certificate is evidence of feasibility, not a proof
of arbitrary-c physics.
