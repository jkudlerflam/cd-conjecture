# Phase diagram assignment for the simple chirality chain

Date: 2026-11-18. OpenAI Codex coordination. **[Research assignment;
no new phase claim.]** The user requested one new project chat to
investigate the phase diagram as a function of h of

    H(h) = sum_j [(9/16) X_j Y_(j+1) + (7/16) Y_j X_(j+1)
                  + Z_j Z_(j+1)] - h sum_j Z_j.

X, Y and Z are Pauli matrices. The field is counted once per site,
including both ends of an open chain. The three bond coefficients stay
fixed; only h varies. The microscopic dimension is d=2.

## Deployment

One new chat, **Chirality chain phase diagram versus field**, was created
in the c-d conjecture project using the user's requested GPT-6 Astra,
Max reasoning (`gpt-6-astra`, `max`). Its ID is
`01a1f0f8-ad15-77a2-8005-ae48e42307b6`, on the local host. A progress
snapshot confirmed the assigned turn in progress, with an initial plan
to check exact mappings, audit saved ladders and perform a bounded
one-worker spectrum survey. No additional chat or subagent is assigned.

The chat owns `contrib/qubit_chirality_field_diagram_20261118/`.
The preceding four stability chats completed their earlier assignments;
this is a new investigation of the simpler Hamiltonian.

## Scientific questions

1. Establish the exact symmetries, any h to -h relation, and useful
   anchors at h=0 and large absolute h. Do not assume an exactly
   saturated product state at finite field in the presence of pairing.
2. Survey the full field range with bounded calculations, refining
   candidate boundaries. Distinguish gaps from splitting among ordered
   vacua; examine magnetization, Neel and bond order, staggered chirality
   and incommensurate correlations without imposing an ordered unit cell.
3. Reanalyze the existing ladders around h=1.738 to 1.755 with explicit
   fit and convergence criteria. Determine what distinguishes a genuine
   c about 2 region or point from crossover, metastability and finite
   entanglement effects. Account for unequal velocities and multiple
   correlation lengths.
4. Independently test the proposed boson-times-Ising transition near
   h=1.7466 through operator content, connected chirality correlations,
   spectra, exponents or velocities. A transition located by imposing
   the Ising exponent must not use the same fitted exponent as independent
   confirmation. Weak first-order behavior remains an alternative when
   available scales cannot exclude it.

The expected output is a phase table, a static diagram with uncertain
intervals explicitly marked, reproducible small result files and a
concrete analytic or computational advance. Unsupported intervals must
remain unresolved. Candidate central charges are not established CFTs.

## Sources and distinctions

The main input is the collaborator's
[chirality-Ising contribution](../d2_xy_yx_chirality_imps_20261118/README.md),
including its target and control ladders. It reports c fits 1.53 to 1.59
near h=1.7466, apparent flow toward c=1 at h=7/4, and unresolved c about
1.9 to 2.2 near h=1.739 to 1.741. These are numerical candidates. The
underlying MPS states are not stored in Git.

Earlier sources include the
[V=4 bulk study](../v4_bulk_diagnosis_20261117/README.md),
[simple-bond interaction study](../qubit_stable_two_modes_20261118/README.md)
and [broader floating-Ising searches](../d2_floating_ising/README.md).
Their conventions and parameter points must be checked before comparing
results. Negative searches on other parameter cuts do not exclude this
model. The phase-1.53 p153c Hamiltonian is a different target, so its
matched quantities and exponents cannot be transferred here.

## Work and publication constraints

The earlier battery-saving stop on long DMRG and background production
continues. Begin with analytic work, saved data and bounded one-worker
calculations. No cluster submissions, detached jobs, self-resubmitting
scripts or long production are authorized by this deployment. If decisive
evidence requires production, complete independent work and specify a
minimal run protocol with resources and stopping criteria first.

Follow AGENTS.md and CONTRIBUTING.md. Preserve prior contributions and
stage only this chat's output. Serialize Git mutations with the shared
atomic directory lock and release only the owning chat's lock. Publish
completed checkpoints by commit, pull with rebase and push. The coordinator
owns shared summaries during deployment; proposed shared text belongs in
the new contribution folder.
