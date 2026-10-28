# A spanning-tree norm for the matched initial kernels

Status: **[proved geometric extension of the abstract polymer theorem;
independently AI-audited]**. This strengthens the output norm of
[INFINITE_CELL_PROGRESS.md](INFINITE_CELL_PROGRESS.md). It does not prove
the subsequent multiscale RG.

Author: coordinating GPT-6-family Codex agent, 2026-10-28. The construction
below uses the already-reviewed polymer representation. A separate
alternatives agent checked the network geometry, cluster weights,
continuous-time norm and periodization, with one counting clarification
incorporated. No human review.

## Statement

For a nonconstant kernel with external spacetime points
`P={(x_i,t_i):1<=i<=p}`, define `T_(a,kappa)(P)` as the minimum spanning-tree
length in the metric

    d((x,t),(y,s))=a|x-y|+kappa|t-s|.

At finite inverse temperature use circular time distance. Repeated
external points contribute no distance. Define the pinned coefficient
norm exactly as in the matching theorem, replacing its diameter weight by

    exp[T_(a,kappa)(P)].

This weight controls all field degrees simultaneously. It is stronger
than a diameter weight, and no estimate of the form `T<=p*diameter` is
used.

Fix `a>=0`, choose

    L>=8(a+log(256)+4),
    tau>=2L/g,             Tr_E exp(-tau a_aux)<=e,
    kappa=L/(32tau),       v_*=1/(8tau).

Here `a_aux` is the auxiliary onsite operator; it is distinct from the
spatial exponent `a`. At zero temperature the heat-trace hypothesis is
unnecessary. The same analytic smallness condition on the input
coefficient norm as in the thermal matching theorem gives a finite
bound in this stronger output norm, uniformly in volume, temperature
`beta>=2tau`, auxiliary cutoff and field degree. Its constants may
depend on the fixed code, gap, heat-trace time and input weight.

The leakage-linear Cauchy estimate, exponential thermal baseline, and
thermal-to-periodized-vacuum comparison consequently hold in this norm
as well, after enlarging their fixed constants. Scalar coefficients
retain their separate density convention.

## Geometry of one connected cluster

In the thermal proof, time blocks have widths in `[tau,2tau]`. A cluster
is a connected family of polymers under overlap-or-neighbor
incompatibility. Let `N` be the sum of their occupied-node counts,
including multiplicities. The union of their occupied square-grid
nodes is connected and has at most `N` nodes.

Draw a vertical segment at each occupied cell position, covering that
node's full time block. To accommodate normal-order point splitting,
extend each end by at most `tau/2`; its length is at most `3tau`.
Draw a horizontal unit segment between every pair of horizontally
adjacent occupied nodes, at the common block midpoint. There are at
most `N` such horizontal segments: each row with `m` nodes has at most
`m` adjacent pairs, including on a periodic spatial ring. Vertically
adjacent bars meet at the block
boundary. Thus the union is a connected network containing every
external point of the cluster. Its weighted length is at most

    (a+3kappa*tau)N.                                      (1)

An active interaction supported on an interval marks every cell of that
interval in its block. Its external fields therefore lie on these bars,
irrespective of how many fields or insertion times it has. Passive nodes
and incompatibility links supply the remaining connections. This is
why (1) has no field-degree factor.

For a cluster on the line, take a finite connected spanning subnetwork,
remove loops and traverse its tree twice. The resulting walk visits all
external points and has length at most twice the network length.
Shortcutting between successive external points in the metric gives a
connected terminal graph; removing its cycles gives a spanning tree no
longer than the walk. The same argument uses shortest circular-time
distances for a cluster on the cylinder. Consequently

    T_(a,kappa)(P) <= (2a+6kappa*tau)N.                    (2)

No assumption that the external points themselves connect the polymer
is required: the network may have auxiliary Steiner vertices. The factor
two is what converts that network into the terminal minimum spanning
tree.

## Cluster summation with the stronger weight

The thermal decoration bound for a fixed occupied support is
`exp[-(L-2)N]`. Combining (2), the rooted animal count `16^N`, the usual
cluster reserve `exp(N)` and one additional reserve `exp(N)` leaves the
geometric ratio

    16 exp[-L+4+2a+6kappa*tau].                           (3)

Since `6kappa*tau=3L/16` and
`2a<=L/4-2log(256)-8`, this ratio is smaller than `1/16`.
The same `5N` incompatibility-neighborhood estimate therefore verifies
the Kotecky--Preiss criterion with spare exponential decay.

Strict submultiplicativity of the terminal-tree weight for disconnected
coefficient products is not asserted. Instead, apply the coefficient
norm majorant before adding this weight, and charge the weight of each
completed connected cluster to its total occupied size using (2).
The resulting scalar polymer majorant is exactly (3). Thus the proof
does not silently assume that two disconnected external sets can be
joined at zero cost.

The connected coefficient density converts to the pinned norm by the
same translation-invariance argument as in the original theorem. The
tree metric is invariant under common spacetime translations. The
point-splitting limit is dominated uniformly by (3).

For periodization, the circular distance between any two points is at
most the line distance of any chosen lifts. Applying this to a minimum
tree of those lifts gives the necessary norm comparison. Clusters that
cannot be lifted through an empty time layer have total size at least
`floor(beta/tau)` and are suppressed by the additional exponential
reserve. This proves the stronger-norm periodization estimate, including
the same antiperiodic signs as in the original theorem.

Analyticity in the input interactions and Cauchy's formula now give the
leakage-linear estimate in the tree-weighted output norm. Bounded source
derivatives use a smaller analytic ball. Complex length tilts may be
included using [TILTED_MATCHING.md](TILTED_MATCHING.md), with a common
weaker gap and heat-trace time.

## Why the refinement matters

Some constructive RG formulations impose exponential decay in a shortest
connecting-tree length for initial retarded vertices. A diameter bound
alone would lose a factor growing with degree when converted to that
norm. The occupied spacetime network supplies the stronger bound directly.
It controls endpoint locality; the normalized multiscale determinant,
localization, running-coupling and source estimates remain separate work.
