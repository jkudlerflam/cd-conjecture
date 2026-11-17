# Sequential decisions

1. Saved L=96,chi=256 checkpoints are present (about 79 MB each).
   The prior raw outputs identify exactly V=4,kappa=a=1, full endpoint
   fields, both parity sectors. They will be read without modification.
2. The actual V=4 fixed-work L=144,chi=256 benchmark uses 0.793 GB RSS
   and about 0.55 seconds for a central 100-matvec update. Prior L=96
   chi=256 full-sweep timings, scaled by 144/96, remain below the
   registered 10-minute cap even at the observed maximum. This supports
   L=144,chi=256 production; it does not pre-approve chi=320.
3. All new density/bond/pair contractions, including adjacent overlapping
   bond products, agree with literal L=8 ED in both parities to below
   2e-9. The validation uses the physical bond in the computational gauge;
   density and entropy are invariant under both onsite changes of basis.
