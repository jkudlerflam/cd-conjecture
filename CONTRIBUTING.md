# Contributing

## Day-to-day workflow

If you have write access, work on `main` and keep the history linear:

```bash
git pull --rebase
```

Then make your changes, commit them in small pieces, and push:

```bash
git push
```

If someone pushed first, `git pull --rebase` again and resolve any conflicts before pushing. Without write access, fork the repository and open a pull request.

Use issues for open problems, proposed attacks and claims that need checking, so that discussion stays attached to the question.

## Adding a result

Put each new piece of work in its own folder, `contrib/<short-name>/`, with a `README.md` that states the following.

1. The claim, with every hypothesis written out.
2. Its status: proved, numerical, heuristic, or refuted.
3. How it was checked and by whom, including any AI assistance and the model used.
4. How to reproduce any numerics: the script, its parameters and its runtime.

If the result changes the overall picture, add a line to "Where things stand" in `README.md` and a short section to `REPORT.md`.

## Conventions

- Tag claims by how they were established. Use [proved] for a complete proof in the repository, [checked] when someone other than the author has verified it, [numerical] for evidence from finite-size data, [read] when taken from a paper that was opened, and [recalled] when taken from memory.
- Name the central charge being used. The conjecture concerns the Virasoro c of unitary compact CFTs. For non-compact or non-unitary limits, also report c_eff = c − 24 h_min.
- Keep files under about 50 MB (GitHub refuses files over 100 MB). Put larger raw data elsewhere and commit a pointer plus the script that regenerates it.
- Python code should run with numpy and scipy unless a `requirements.txt` says otherwise.
