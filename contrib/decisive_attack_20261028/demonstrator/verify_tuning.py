#!/usr/bin/env python3
"""Exact second-order certificate for the tunable ten-letter demonstrator.

Standard library only. The 36-dimensional excited-pair resolvent is evaluated
with rational arithmetic, independently of routeA_demonstrator/effective.py.
"""
from fractions import Fraction as F
import json
from pathlib import Path


def solve(matrix, rhs):
    a = [row[:] + b[:] for row, b in zip(matrix, rhs)]
    n, nrhs = len(matrix), len(rhs[0])
    for col in range(n):
        pivot = next(row for row in range(col, n) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        divisor = a[col][col]
        a[col] = [v / divisor for v in a[col]]
        for row in range(n):
            if row != col and a[row][col]:
                factor = a[row][col]
                a[row] = [x - factor*y for x, y in zip(a[row], a[col])]
    return [row[n:n+nrhs] for row in a]


def certificate(d, kappa, jp, jz):
    # Excited states: p0,p1,a0,a1,b0,b1; p=Phi-, a=01, b=10.
    h = [[F(int(i == j))*d for j in range(6)] for i in range(6)]
    h[0][3] = h[3][0] = kappa
    h[1][4] = h[4][1] = -kappa
    pair = [[h[i//6][j//6]*int(i%6 == j%6)
             + h[i%6][j%6]*int(i//6 == j//6)
             for j in range(36)] for i in range(36)]
    src = [[F(0) for _ in range(4)] for _ in range(36)]
    for c1 in range(2):
        for c2 in range(2):
            c = 2*c1+c2
            src[6*c1+c2][c] += jz
            src[6*(2+c1)+4+c2][c] += jp/2
            src[6*(4+c1)+2+c2][c] += jp/2
    solution = solve(pair, src)
    actual = [[-sum(src[r][i]*solution[r][j] for r in range(36))
               for j in range(4)] for i in range(4)]
    a = 2*d
    r1 = 1/a
    r2 = a/(a*a-kappa*kappa)
    r4 = (a*a-2*kappa*kappa)/(a*(a*a-4*kappa*kappa))
    scalar = -jz*jz*r4-jp*jp*(2*r2+r1+r4)/8
    xy = jp*jz*kappa*kappa/(4*d*(d*d-kappa*kappa))
    zz = 3*jp*jp*kappa**4/(32*d*(d*d-kappa*kappa)*(4*d*d-kappa*kappa))
    expected = [[F(0) for _ in range(4)] for _ in range(4)]
    for i, sign in enumerate([1,-1,-1,1]):
        expected[i][i] = scalar+sign*zz
    expected[1][2] = expected[2][1] = xy
    assert actual == expected
    return {"D": str(d), "kappa": str(kappa), "Jp": str(jp), "Jz": str(jz),
            "Jxy": str(xy), "Jzz": str(zz),
            "Delta_abs": str(abs(2*zz/xy)) if xy else None,
            "exact_matrix_match": True}


def main():
    rows = []
    for d,k in [(F(1),F(1,5)),(F(1),F(3,5)),(F(2),F(1,2))]:
        for jp,jz in [(F(2),F(-1)),(F(1,10),F(-1)),(F(0),F(1)),(F(1),F(0))]:
            rows.append(certificate(d,k,jp,jz))
    result = {"status": "all exact rational assertions passed", "cases": rows}
    target = Path(__file__).with_name("tuning_checks.json")
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(f"Passed {len(rows)} exact 36-dimensional resolvent certificates; wrote {target.name}.")


if __name__ == "__main__":
    main()
