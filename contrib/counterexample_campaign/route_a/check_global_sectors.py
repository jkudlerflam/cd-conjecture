"""Rational n=23 selector checks and finite diagnostics for README.md.

Run from any directory. Writes only this contribution's checks.json.
The all-length and nonzero-lambda statements are analytic proofs in README.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import math
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import scipy.linalg as la

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def interval_add(x, y):
    return x[0]+y[0], x[1]+y[1]


def interval_scale(c, x):
    z = c*x[0], c*x[1]
    return min(z), max(z)


def interval_sub(x, y):
    return interval_add(x, interval_scale(-1, y))


def atan_interval(x, terms):
    s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    nxt = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s, s+nxt), max(s, s+nxt)


def cos_interval_at(x, terms=35):
    s = sum(((-1)**j*x**(2*j)/F(math.factorial(2*j))
             for j in range(terms)), F(0))
    nxt = (-1)**terms*x**(2*terms)/F(math.factorial(2*terms))
    return min(s, s+nxt), max(s, s+nxt)


def cos_interval(x):
    assert 0 < x[0] <= x[1] < 2
    return cos_interval_at(x[1])[0], cos_interval_at(x[0])[1]


def rational_parameters():
    # Machin's formula; both arctan series have decreasing alternating terms.
    pi = interval_sub(interval_scale(16, atan_interval(F(1, 5), 65)),
                      interval_scale(4, atan_interval(F(1, 239), 22)))
    assert F(314159, 100000) < pi[0] < pi[1] < F(314160, 100000)
    aa = {}
    for b in (23, 24, 25):
        c = cos_interval(interval_scale(F(1, b+1), pi))
        aa[b] = interval_sub((F(2), F(2)), interval_scale(2, c))
    db = interval_sub(interval_scale(25, aa[24]), interval_scale(24, aa[25]))
    dprev = interval_sub(interval_scale(24, aa[23]), interval_scale(23, aa[24]))
    delta = interval_scale(F(1, 2), interval_add(db, dprev))
    gamma = interval_scale(F(1, 48), interval_sub(dprev, db))
    mu = interval_sub((F(2), F(2)), delta)
    eb = interval_scale(F(1, 24), interval_sub(aa[24], delta))
    # Delta-2 = 2-2cos(pi/25) = a_24; beta<1 suffices for t_d.
    kappa = interval_sub(aa[24], (F(1, 2**136)+F(1, 2**365),)*2)
    lmax = F(1, 2**20)
    B_upper = 2*lmax  # sqrt(3)<2, avoiding numerical radicals.
    assert mu[0] > 0 and eb[1] < 0 and kappa[0] > 0
    assert gamma[0] > 5*B_upper
    assert F(10) > F(4)+7*B_upper
    def below(x):
        return F((x.numerator*2**128)//x.denominator, 2**128)
    assert below(gamma[0]) > 5*B_upper and below(kappa[0]) > 0
    return dict(
        status="proved_rational_inequalities", b=24,
        gamma_interval=[float(x) for x in gamma],
        marker_potential_interval=[float(x) for x in mu],
        energy_density_interval=[float(x) for x in eb],
        excited_fragment_gap_lower=float(kappa[0]),
        lambda_max_power=-20, K=10,
        wrong_length_sector_gap_lower=float(gamma[0]-5*B_upper),
        illegal_sector_gap_lower=float(F(6)-7*B_upper),
        rational_witnesses={
            "gamma_lower": str(below(gamma[0])),
            "five_B_upper": str(5*B_upper),
            "kappa_lower": str(below(kappa[0])),
            "K_minus_four_minus_seven_B_upper": str(F(6)-7*B_upper),
        })


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def finite_diagnostics():
    cellmod = module(ROOT / "contrib/routeA_spatial/ycell.py", "route_a_ycell")
    p = cellmod.params(Delta=2.5, Dx=3.0, DM=4.0, delta=0.0,
        gp=.2, gm=-.15, td=.1, eps=.04, musel=.7,
        phs={"0": .13, "+": -.2, "-": .5},
        phx={"0": .7, "+": .6, "-": -.4},
        phsw={"0": .3, "+": -.8, "-": .9},
        phd={"0": .2, "+": .11, "-": -.43})
    records = []
    for n in range(5):
        c = cellmod.build_cell(n, p)
        H, exc = c["H"], c["exc"]
        mine = float(la.eigvalsh(H[np.ix_(exc, exc)], subset_by_index=[0,0])[0])
        minv = float(la.eigvalsh(H[np.ix_(~exc, ~exc)], subset_by_index=[0,0])[0])
        expected = p["musel"]-2*math.cos(math.pi/(n+2))
        assert abs(minv-expected) < 1e-11
        assert mine >= p["musel"]+.2-1e-11
        records.append(dict(n=n, dim=len(H), vacuum_min=minv,
                            excited_min=mine, proved_excited_lower=.9))
    fv = np.zeros((4,4), complex)
    fv[0,1] = fv[1,0] = 1
    fp = np.zeros((4,4), complex)
    fp[2,0] = fp[0,3] = 1
    fm = fp.conj().T
    V = np.kron(fv, fv)+1j*np.kron(fp, fm)-1j*np.kron(fm, fp)
    eig = la.eigvalsh(V)
    assert abs(max(abs(eig))-math.sqrt(3)) < 1e-12
    # Exact support test of all literal moves: inner pairs are matched and
    # external arrows are retained. This also checks the enlarged express rail.
    arrow = {"L": (0,0), "R": (1,1), "Xd": (0,1),
             "Xx": (0,1), "Me": (1,0), "M": (0,0)}
    entries = 0
    for (li,ri), outputs in cellmod.bond_moves(p).items():
        ai,bi = arrow[cellmod.LET[li][0]], arrow[cellmod.LET[ri][0]]
        assert ai[1] == bi[0]
        for _, (lo,ro) in outputs:
            ao,bo = arrow[cellmod.LET[lo][0]], arrow[cellmod.LET[ro][0]]
            assert ao[1] == bo[0] and ao[0] == ai[0] and bo[1] == bi[1]
            entries += 1
    return dict(status="numerical_diagnostics_and_exact_support_checks",
                cells=records, port_norm=float(max(abs(eig))),
                checked_directed_move_entries=entries)


if __name__ == "__main__":
    result = dict(selector=rational_parameters(), finite=finite_diagnostics())
    (HERE / "checks.json").write_text(json.dumps(result, indent=2)+"\n")
    compact = {"selector": {k:v for k,v in result["selector"].items()
                            if k != "rational_witnesses"}, "finite": result["finite"]}
    print(json.dumps(compact, indent=2))
