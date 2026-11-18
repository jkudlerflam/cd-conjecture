# Project instructions for every chat

Follow [CONTRIBUTING.md](CONTRIBUTING.md) for research and contribution conventions.

## Research priorities

The user's current highest-priority directions are large central charge at
fixed microscopic d=2, and unbounded central charge at fixed small d.
Prefer concrete progress on these directions over a general upper-bound
program. For the latter direction, the markerless d=9 frozen / d=10 mobile
family is the primary current construction; its finite certificates do not
yet establish arbitrary n or a microscopic unbounded-c theorem.

The user explicitly values simple Hamiltonians: fewer terms, ratios of
nonzero coefficients closer to one, and less delicate tuning. When ranking
candidate models, weigh these alongside local dimension and strength of
evidence. A small on-site dimension alone does not make a construction simple.

The user also explicitly prioritizes the smallest possible physical on-site
dimension. Search and audit smaller-d constructions before extending a larger-d
version when they can address the same research objective. Distinguish a
certified finite member from an arbitrary-c family at fixed d, and retain all
infrared and sector-selection assumptions when comparing them. A smaller
alphabet with only a finite-member certificate is not automatically an
unbounded-c construction. Count every physical state, including markers and
vacancies; do not substitute an encoded or blocked dimension for microscopic d.

For comparisons, distinguish independent real coupling parameters from
literal operator monomials in a stated basis. Remove irrelevant constants,
telescoping bond terms, and a common positive energy scale before reporting
coefficient ratios. Grouping terms or changing notation does not by itself
simplify the physical model. Distinguish microscopic coefficient hierarchy
from emergent small velocities, and prefer robust parameter regions over
long-decimal tuned points when the scientific evidence is comparable.

## Keep GitHub synchronized

The user has authorized routine pulls, commits, and pushes for work in this
project. Apply this workflow in every chat without asking for permission again
for ordinary synchronization.

- At the start of work and when resuming after a break, inspect `git status`
  and the current branch, then run `git pull --rebase` against its upstream
  when the working tree is ready. Preserve existing edits first.
- During longer work, synchronize at each meaningful checkpoint and roughly
  every 30 minutes at the next safe stopping point. Commit completed work in
  small, descriptive commits, pull with rebase, and push. Do not wait until
  the end of a long chat to share completed work.
- Before the final response, commit the task's completed changes, run the
  relevant checks, pull with rebase, and push. Verify that the local branch
  and its upstream agree, and report anything still uncommitted or unpushed.
- This repository's normal collaboration branch is `main`. Respect an
  explicitly requested branch or an existing isolated worktree; push that
  branch to its matching upstream instead of changing branches unexpectedly.
- If a push is rejected because the remote advanced, pull with rebase,
  resolve conflicts while preserving both contributors' work, repeat any
  affected checks, and push again. Never force-push shared history.
- Multiple chats may share a checkout. Stage only the files or hunks belonging
  to the current task unless the user requests a broader sync. Do not discard,
  reset, or silently stash another chat's edits, and do not run overlapping
  Git mutations. If another chat's work prevents a safe pull, fetch the latest
  changes and complete the integration at the next safe checkpoint.
- Keep credentials, ignored checkpoints, and oversized raw data out of Git.
  If access, network failures, or an ambiguous conflict prevents a safe sync,
  preserve the work and explain the specific blocker rather than claiming
  the repository is synchronized.
