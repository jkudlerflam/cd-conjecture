"""Finite checks of the C92 entropy channels. All outputs stay here."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from itertools import product
from pathlib import Path
import math
import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent


def gap_basis(n, vmax):
    return [g for g in product(range(vmax + 1), repeat=n) if sum(g) <= vmax]


def bath(gaps, t, tilt):
    n = len(gaps[0])
    ix = {g: i for i, g in enumerate(gaps)}
    h = np.diag([tilt * sum(g) for g in gaps]).astype(float)
    for i, g in enumerate(gaps):
        for j in range(n):
            if g[j]:
                for k in ((j - 1) % n, (j + 1) % n):
                    out = list(g)
                    out[j] -= 1
                    out[k] += 1
                    h[ix[tuple(out)], i] -= t
    return h


def entropy(r):
    vals = np.linalg.eigvalsh((r + r.conj().T) / 2)
    vals = vals[vals > 1e-14]
    return float(-np.sum(vals * np.log(vals)))


def gibbs(h, temp):
    vals, vecs = eigh(h)
    w = np.exp(-(vals - vals[0]) / temp)
    free = float(vals[0] - temp * np.log(w.sum()))
    r = (vecs * (w / w.sum())) @ vecs.conj().T
    return free, r


def expectation(r, h):
    return float(np.trace(r @ h).real)


def contact_channel_checks():
    x = np.array([[0, 1], [1, 0]], dtype=float)
    z = np.diag([1., -1.])
    sources = [2 * z + x, -z + 1.7 * x, .4 * z - 1.2 * x]
    m = max(np.linalg.norm(h, 2) for h in sources)
    rows = []
    max_entropy_violation = 0.
    min_variational = math.inf
    max_upper_violation = 0.
    for vmax, tilt, coupling, temp in product((3, 5), (.55, .7), (.01, .04), (.02, .08)):
        gaps = gap_basis(3, vmax)
        ix = {g: i for i, g in enumerate(gaps)}
        h = np.kron(bath(gaps, 1., tilt), np.eye(2))
        for j in range(3):
            a = np.diag([float(g[j] == 0) for g in gaps])
            h += coupling * np.kron(a, sources[j])
        free, r = gibbs(h, temp)
        c = tilt - coupling * m
        assert c > 0 and temp <= c / (4 * math.log(2))
        lower = min(.5, c * c / 256)
        for j in range(3):
            a = np.diag([float(g[j] == 0) for g in gaps])
            p = np.kron(a, np.eye(2))
            down = np.zeros((len(gaps), len(gaps)))
            for i, g in enumerate(gaps):
                if g[j]:
                    out = list(g)
                    out[j] -= 1
                    down[ix[tuple(out)], i] = 1
            d = np.kron(down, np.eye(2))
            trial = p @ r @ p + d @ r @ d.T
            contact = expectation(r, p)
            q = 1 - contact
            h2 = -contact * math.log(contact) - q * math.log(q)
            entropy_loss = entropy(r) - entropy(trial)
            delta = expectation(trial, h) - temp * entropy(trial) - free
            upper = -c * q + 4 * math.sqrt(contact * q) + temp * h2
            assert abs(np.trace(trial) - 1) < 1e-11
            assert contact >= lower - 1e-10
            assert entropy_loss <= h2 + 1e-10
            assert delta >= -1e-10
            assert delta <= upper + 1e-10
            max_entropy_violation = max(max_entropy_violation, entropy_loss - h2)
            max_upper_violation = max(max_upper_violation, delta - upper)
            min_variational = min(min_variational, delta)
            rows.append((h.shape[0], contact / lower))
    return len(rows), max(r[0] for r in rows), min(r[1] for r in rows), min_variational, max_entropy_violation, max_upper_violation


def pin_checks():
    rows = []
    b = 3
    mu = .8
    # One pinned b=3 vacuum cell, with a two-dimensional spectator bit.
    # The unpinned comparison retains that bit and has a three-state clock.
    # This tests the complete thermal trace and the pure-clock reset map.
    for vmax, tilt, coupling, temp in product((3, 5), (.55, .7), (.03, .08), (.01, .04)):
        gaps = gap_basis(b, vmax)
        a = [np.diag([float(g[j] == 0) for g in gaps]) for j in range(b)]
        common = bath(gaps, 1., tilt) + coupling * mu / 2 * (a[1] + a[2])
        hp = np.kron(common, np.eye(2))
        hu = np.kron(common, np.eye(2 * b))
        for j in range(b - 1):
            link = np.zeros((b, b))
            link[j, j + 1] = link[j + 1, j] = -1
            hu += coupling * np.kron(a[j], np.kron(np.eye(2), link))
        fp, rp = gibbs(hp, temp)
        fu, _ = gibbs(hu, temp)
        contact_edges = [expectation(rp, np.kron(a[j], np.eye(2))) for j in range(b - 1)]
        path = np.diag(contact_edges, 1) + np.diag(contact_edges, -1)
        top = np.linalg.eigvalsh(path)[-1]
        # Resetting only the clock preserves the spectator entropy in this model.
        margin = fp - fu - coupling * top
        assert margin >= -1e-10
        rows.append((hu.shape[0], margin))
    return len(rows), max(r[0] for r in rows), min(r[1] for r in rows)


def main():
    contact = contact_channel_checks()
    pin = pin_checks()
    text = (
        "# C92 finite checks\n\n"
        "[Numerical] These finite matrix checks test the entropy-channel algebra, not the infinite-volume mobile theorem.\n\n"
        f"Contact-channel cases: {contact[0]}. Maximum matrix dimension: {contact[1]}.\n\n"
        f"Minimum ratio of contact expectation to analytic lower bound: {contact[2]:.8g}.\n\n"
        f"Minimum trial free-energy excess: {contact[3]:.8g}.\n\n"
        f"Maximum entropy-bound violation: {contact[4]:.3g}. Maximum energy-bound violation: {contact[5]:.3g}.\n\n"
        f"Vacuum-clock pin trace cases: {pin[0]}. Maximum matrix dimension: {pin[1]}. Minimum excess above the contact-path bound: {pin[2]:.8g}.\n\n"
        "[Scope] The gap cutoff is on total vacancy number. Vacancy removal maps into that space, so the channel inequality remains exact at the cutoff. These tests use small clocks, zero port and infinite grammar penalty. No scheduled-band or full d=13 simulation is claimed.\n"
    )
    (HERE / "CHECKS.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
