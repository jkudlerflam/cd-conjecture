"""Item 5: is the seven-letter negative result ("at most n mod 2 bands cross the Fermi level")
a statement about zero field only?  The author's check (seven_check.py) fixes f = 0.
Here T and mu are taken from the author's seven-letter two-cell code at the same parameters,
and the uniform field f (an exact chemical potential for the logical fermions, available in the
bond through the charge term) is scanned.  Fermi points are counted by the referee's own method:
changes of the number of negative eigenvalues of h(k) = mu + f + T e^{ik} + T^dag e^{-ik} on a k grid.
"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SEVEN = os.path.join(os.path.dirname(HERE), "seven_letter_negative")
sys.path.insert(0, SEVEN)
os.chdir(SEVEN)
from twocell import analyse, one_particle

def count(mu, T, f, grid=4001):
    ks = np.linspace(-np.pi, np.pi, grid)
    neg = np.array([(np.linalg.eigvalsh(mu + f * np.eye(len(T)) + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)) < 0).sum() for k in ks])
    return int(np.abs(np.diff(neg)).sum())

res = []
for n in (2, 3):
    for label, gL in [("one rung", 0.0), ("two rungs with flux", 0.6 * np.exp(1.1j))]:
        p = dict(lam=1.0, gR=1.0, gL=gL, DA=0.7, DB=0.7, nu=0.0, f=0.0, chi=0.0)
        r = analyse(n, p); h, vacE = one_particle(r["h2"], n, r["bitlist"])
        T = h[:n, n:]; mu = h[:n, :n] + h[n:, n:]
        sc = np.linalg.norm(T, 2)
        best = (0, 0.0); table = []
        for x in np.linspace(-2.5, 2.5, 501):
            c = count(mu, T, x * sc, grid=1201)
            table.append((round(float(x), 3), c))
            if c > best[0]: best = (c, float(x))
        # reciprocity: spectrum(k) vs spectrum(-k)
        ks = np.linspace(0.1, 3.0, 30)
        rec_asym = max(abs(np.linalg.eigvalsh(mu + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k))
                           - np.linalg.eigvalsh(mu + T * np.exp(-1j * k) + T.conj().T * np.exp(1j * k))).max() for k in ks) / sc
        d = dict(n=n, design=label, fermi_points_at_zero_field=count(mu, T, 0.0), needed=2 * n,
                 max_fermi_points_over_field_scan=best[0], at_field_over_normT=best[1],
                 field_window_with_max_count=[min(x for x, c in table if c == best[0]), max(x for x, c in table if c == best[0])],
                 spectrum_k_vs_minus_k_asymmetry=float(rec_asym))
        res.append(d); print(json.dumps(d), flush=True)
json.dump(res, open(os.path.join(HERE, "item5_seven_field.json"), "w"), indent=1)
