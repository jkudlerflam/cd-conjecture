# Project instructions for every chat

Follow [CONTRIBUTING.md](CONTRIBUTING.md) for research and contribution conventions.

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
