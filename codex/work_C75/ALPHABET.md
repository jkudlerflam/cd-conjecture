# Alphabet audit

[Result] The smallest design obtained in this run has 15 letters. There
is no proof that 15 is minimal among all route-A constructions.

| Group in the 23-letter parent | Parent count | New count | Decision |
| :-- | --: | --: | :-- |
| L/R data letters | 4 | 4 | Retained as the domain-wall selector |
| Vacuum data token with a bit | 2 | 0 | Absorbed into R_a |
| Excited data token with a bit | 6 | 4 | Two species replace three |
| Express token with a bit | 6 | 4 | Same two-species reduction |
| Markers | 5 | 3 | Two excited markers plus one merged vacuum marker |
| Total | 23 | 15 | Eight letters removed |

## Vacuum letters and shared data/token encodings

[Proved] Separate vacuum-token letters are unnecessary. Their positions
are encoded by L^j R^(n-j), including the all-L endpoint. The data bit
remains on every data site. The empty and vacuum-occupied markers can
also be identified. This changes exact grammar invariance into a large,
finite penalty construction, which is accounted for in PROOF Sections
6 and 8. Treating the identification as a projection without those
corrections would not define the finite bond being claimed.

[Open] Identifying excited letters such as D_(A,1) and D_(B,0) would
save another letter while preserving their common charge. It also
identifies the endpoints of the conversion. The g=0 propagator would
then mix the two data/species channels at order one, so the present
small-g Gaussian expansion no longer applies. No controlled 14-letter
variant from that identification was found.

## Charged species

[Proved] Three excited species are unnecessary. Two species whose
charges differ by one data charge suffice, with their average charge
equal to the vacuum-token charge. Their conversion flips one data bit,
and the port creates one of each. In the integer charge convention used
here the excited charges are -1,+1 and the data charge is 2.

[Proved, restricted] Keeping only one charged excited species and the
original neutral vacuum cannot implement the old pair-creation port
without a compensating data flip. Such a pair changes Q. Deleting one
of the original two conversions also makes the old g_+ g_- hopping
coefficient vanish. These observations rule out literal deletions,
not every redesigned single-species mechanism. In the new model there
is one conversion term, but two excited species.

## Express track and folding

[Proved, restricted] The express track is still used. Simply identifying
it with the data track changes the folded graph to a path whose sources
see the data sites in opposite orders. The positive Gram argument no
longer applies. There is also no independent scalar loop flux in the
single-path equal-rail version, so its sea centers coincide and its
quartic pair-transfer margin vanishes.

For a concrete second-order obstruction take n=2 and identical real
scalar carriers. Reflection gives the hopping, up to its overall factor,

    F = [[a,b],[c,a]],
    a=<f_0,f_1>, b=||f_1||^2, c=||f_0||^2.

The endpoint filters are linearly independent, so bc>a^2. They have
different short-time powers, and for the endpoint path c>b. The Bloch
eigenvalues are

    2a cos k +/- sqrt((b+c)^2 cos^2 k+(b-c)^2 sin^2 k).

The upper branch is strictly positive and the lower strictly negative
for every k. Their spectral intervals are disjoint. No common chemical
potential crosses both branches. This is the simple deletion's failure,
not a lower bound for all possible folding mechanisms.

[Exploratory, not completed] A fast nearly uniform second carrier could
make the two source orders more alike without its express rail. The
different scalar propagators then lose the cancellation of onsite
quadratic terms, which can be much larger than the hopping. No finite
certificate or local counterterm prescription was obtained for that
alternative.

## L/R letters, markers, and two-token cells

[Proved, within this encoding] L and R now store the vacuum-head position
as well as the data bit. Identifying them removes the n+1-position
selector path. M_A and M_B have different exact charges and cannot be
identified while retaining this U(1). The unexcited marker supplies the
cell boundary and the neutral input of the pair port. These are reasons
to retain the seven L/R and marker letters here, not universal lower
bounds.

[Open] A two-token cell may implement a different fold or encode the
species in the relative token position. It needs a number selector,
the correct boundary source operators and a proof of the Gaussian
second-order structure. No such smaller complete bond was obtained or
numerically certified in this run. The present reduction makes a second
vacuum token unnecessary, but does not exclude that possibility as a
route below 15.
