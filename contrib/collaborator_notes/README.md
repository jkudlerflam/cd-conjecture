# Collaborator notes

Informal write-ups of the two main constructions, for sharing with people who want the mechanism and the evidence without the project record. Each folder has the LaTeX source, the compiled PDF, the figures and the script that draws them (`python3 make_figs.py`, needs numpy and matplotlib).

- [qubit_c2/notes_qubit_c2.pdf](qubit_c2/notes_qubit_c2.pdf): a nearest-neighbour qubit chain with c = 2. A blockade and a virtual pair process generate a next-nearest-neighbour hop, the band has two valleys, and DMRG to L = 512 reads c close to 2. Open points include convergence at L = 512 and whether a number-changing operator locks the charge mode beyond several thousand sites. By an order-of-magnitude estimate this locking would affect the symmetric-line bond but not the phase-1.53 bond.
- [unbounded_c/notes_unbounded_c.pdf](unbounded_c/notes_unbounded_c.pdf): ten states per site with predicted c = n + 1. A soliton lattice of period n + 1 carries n flat bands, and charge moves between cells only through virtual fractionally charged kinks, giving an n-band metal with 2n Fermi points. Worked out for n = 9 to 20, with the infrared step conditional. The technical record is [contrib/markerless_d9_d10](../markerless_d9_d10/README.md).

Written 2026-11-18. Every number traces to REPORT.md, STATUS.md or `contrib/markerless_d9_d10/`, cited in comments in the LaTeX source. Both notes were checked against those sources by an independent AI reviewer; no human has reviewed them.
