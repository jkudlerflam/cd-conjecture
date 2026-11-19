"""Rebuild the nine-orbital scalar Bloch matrix in d9_construction_notes.tex.

Run with python3 -B check_d9.py. The optional argument sets decimal precision.
Only the Python standard library is used. All energies use unit walker hopping.
"""

import math
import sys
from decimal import Decimal as D, getcontext


def band_factor(matrix, width):
    """In-place real Cholesky. Row i stores (i,i-offset) by offset."""
    for i, row in enumerate(matrix):
        for j in range(max(0, i - width), i + 1):
            other = matrix[j]
            value = row[i - j]
            for k in range(max(0, i - width), j):
                value -= row[i - k] * other[j - k]
            row[i - j] = value.sqrt() if i == j else value / other[0]


def band_solve(matrix, width, rhs):
    x = rhs.copy()
    for i, row in enumerate(matrix):
        value = x[i]
        for offset in range(1, min(i, width) + 1):
            value -= row[offset] * x[i - offset]
        x[i] = value / row[0]
    for i in reversed(range(len(x))):
        value = x[i]
        for j in range(i + 1, min(len(x), i + width + 1)):
            value -= matrix[j][j - i] * x[j]
        x[i] = value / matrix[i][0]
    return x


def scalar_blocks(n, beta, td):
    # Order each particle as bright_0, dark_0, bright_1, dark_1, ...
    m = 2 * (n + 1)
    vertices = m * m
    width = 4 * m
    zero = D(0)
    band = [[zero] * (width + 1) for _ in range(2 * vertices)]
    for p in range(vertices):
        l, r = divmod(p, m)
        diag = D(8 + 1000 * ((l % 2) + (r % 2)))
        band[2*p][0] = band[2*p+1][0] = diag
    # Both species gauge to real rails, with H[dark_n,bright_n] = -i td.
    edges = [(j + 2, j, D(-1), zero) for j in range(m - 2)]
    edges.extend([(1, 0, D(-1), zero), (m - 1, m - 2, zero, -td)])
    for high, low, re, im in edges:
        for other in range(m):
            for i, j in [(m*high + other, m*low + other),
                         (m*other + high, m*other + low)]:
                offset = 2*(i - j)
                band[2*i][offset] += re
                band[2*i+1][offset] += re
                if im:
                    band[2*i][offset-1] -= im
                    band[2*i+1][offset+1] += im
    band_factor(band, width)
    source = [zero] * (2*vertices)
    source[2*((m-1)*m)] = D(1)
    wave = band_solve(band, width, source)
    indices = [[m*2*j + k for k in range(m)] for j in range(n)]
    indices += [[m*k + 2*j for k in range(m)] for j in range(n)]
    gram = [[None] * (2*n) for _ in indices]
    for col, support in enumerate(indices):
        rhs = [zero] * (2*vertices)
        for p in support:
            rhs[2*p:2*p+2] = wave[2*p:2*p+2]
        solution = band_solve(band, width, rhs)
        for row, sample in enumerate(indices):
            re = im = zero
            for p in sample:
                a, b = wave[2*p:2*p+2]
                c, d = solution[2*p:2*p+2]
                re += a*c + b*d
                im += a*d - b*c
            gram[row][col] = (re, im)
    scale = -2 / beta**2
    hop = [[scale * gram[j][n+k][0] for k in range(n)] for j in range(n)]
    onsite = [[scale * (gram[j][k][1] + gram[n+j][n+k][1])
               for k in range(n)] for j in range(n)]
    # Use the antisymmetric part to remove roundoff on the diagonal of mu.
    onsite = [[(onsite[j][k] - onsite[k][j])/2 for k in range(n)] for j in range(n)]
    return hop, onsite


def determinant(matrix):
    """Real part of a complex determinant by row-pivoted elimination."""
    rows = [row.copy() for row in matrix]
    re, im = D(1), D(0)
    for j in range(len(rows)):
        idx = max(range(j, len(rows)), key=lambda i: max(map(abs, rows[i][j])))
        if idx != j:
            rows[idx], rows[j] = rows[j], rows[idx]
            re, im = -re, -im
        a, b = rows[j][j]
        re, im = re*a - im*b, re*b + im*a
        if not re and not im:
            return D(0)
        norm = a*a + b*b
        for i in range(j+1, len(rows)):
            c, d = rows[i][j]
            x, y = (c*a + d*b)/norm, (d*a - c*b)/norm
            for k in range(j+1, len(rows)):
                e, f = rows[j][k]
                u, v = rows[i][k]
                rows[i][k] = (u - x*e + y*f, v - x*f - y*e)
    return re


def null_slope(h, hp):
    """Solve the first n-1 pivot equations, then evaluate q* h' q / (q* q)."""
    n = len(h)
    rows = [row.copy() for row in h]
    for j in range(n-1):
        idx = max(range(j, n), key=lambda i: max(map(abs, rows[i][j])))
        rows[idx], rows[j] = rows[j], rows[idx]
        a, b = rows[j][j]
        norm = a*a + b*b
        for i in range(j+1, n):
            c, d = rows[i][j]
            x, y = (c*a + d*b)/norm, (d*a - c*b)/norm
            for k in range(j+1, n):
                e, f = rows[j][k]
                u, v = rows[i][k]
                rows[i][k] = (u - x*e + y*f, v - x*f - y*e)
    q = [(D(0), D(0)) for _ in h]
    q[-1] = (D(1), D(0))
    for j in reversed(range(n-1)):
        x = y = D(0)
        for (e, f), (u, v) in zip(rows[j][j+1:], q[j+1:]):
            x -= e*u - f*v
            y -= e*v + f*u
        a, b = rows[j][j]
        norm = a*a + b*b
        q[j] = ((x*a + y*b)/norm, (y*a - x*b)/norm)
    ray = D(0)
    for (a, b), row in zip(q, hp):
        for (c, d), (u, v) in zip(row, q):
            ray += a*(c*u-d*v) + b*(c*v+d*u)
    return ray / sum(a*a+b*b for a,b in q)


def crossings(hop, onsite, field):
    n = len(hop)
    even = [[hop[i][j] + hop[j][i] for j in range(n)] for i in range(n)]
    odd = [[hop[i][j] - hop[j][i] for j in range(n)] for i in range(n)]

    def matrices(t, sign):
        c = sign*(1-t*t)/(1+t*t)
        s = sign*2*t/(1+t*t)
        h = [[(c*even[i][j] + (field if i == j else 0), onsite[i][j] + s*odd[i][j])
              for j in range(n)] for i in range(n)]
        hp = [[(-s*even[i][j], c*odd[i][j]) for j in range(n)] for i in range(n)]
        return h, hp

    roots = []
    # Two rational charts cover the circle. The numerator has degree at most
    # 2n in z, so 2n distinct sign changes exhaust it.
    for sign in (1, -1):
        grid = [D(i)/128 for i in range(-128, 129)]
        signs = [determinant(matrices(t, sign)[0]) > 0 for t in grid]
        for i in range(1, len(grid)):
            if signs[i-1] == signs[i]:
                continue
            left, right = grid[i-1:i+1]
            for _ in range(120):
                mid = (left + right)/2
                if (determinant(matrices(mid, sign)[0]) > 0) == signs[i-1]:
                    left = mid
                else:
                    right = mid
            t = (left + right)/2
            h, hp = matrices(t, sign)
            # Differentiate det(h-EI) at a simple zero.
            dk = sum(determinant(h[:j] + [hp[j]] + h[j+1:]) for j in range(n))
            cof = sum(determinant([row[:j]+row[j+1:] for row in h[:j]+h[j+1:]])
                      for j in range(n))
            slope = dk/cof
            assert abs(null_slope(h, hp)/slope - 1) < D('1e-23')
            angle = (2*math.atan(float(t)) + (math.pi if sign == -1 else 0)) % (2*math.pi)
            roots.append((angle, slope))
    return sorted(roots)


def main():
    getcontext().prec = int(sys.argv[1]) if len(sys.argv) > 1 else 102
    n = 9
    a, b = 1, 1004
    for _ in range(n):
        a, b = b, 1004*b - a
    beta = 1/D(b)
    td = beta / 2**35
    alpha = math.sqrt(2/(n+2)) * math.sin(math.pi/(n+2))
    log_unit = -2*(359+163)*math.log10(2) + 4*math.log10(alpha) + 2*float(beta.log10())
    gap = 2*(math.cos(math.pi/(n+2)) - math.cos(2*math.pi/(n+2)))
    r = 10**(-2.74)
    print(f"clock gap {gap:.12g}, log10 E_u {log_unit:.8f}")
    print(f"supplied geometric model r {r:.9g}, sum/v1 {(1-r**n)/(1-r):.9f}")
    hop, onsite = scalar_blocks(n, beta, td)
    roots = crossings(hop, onsite, D(2)**(-83))
    assert len(roots) == 2*n and sum(v > 0 for _, v in roots) == n
    print("angle in radians, signed slope in E_u per cell wave number")
    for angle, slope in roots:
        print(f"{angle:.10f}  {slope:+.8E}")
    speeds = [abs(v) for _, v in roots]
    sep = min((roots[(i+1) % len(roots)][0] - roots[i][0]) % (2*math.pi)
              for i in range(len(roots)))
    print(f"least root spacing {sep:.9g}")
    sums = sorted((a[0]+b[0]) % (2*math.pi) for i, a in enumerate(roots) for b in roots[i:])
    pair_sep = min((sums[(i+1) % len(sums)] - sums[i]) % (2*math.pi) for i in range(len(sums)))
    print(f"least pair-sum spacing {pair_sep:.9g}")
    print(f"log10(min/max |slope|) {(min(speeds)/max(speeds)).log10():.8f}")
    print(f"A_free/(cell length times E_u) {sum(speeds)/2:.9E}")
    # The reported fit is a width hierarchy. For comparison compute the
    # extreme eigenvalues of the symmetric hopping part by inertia bracketing.
    s = [[-(hop[i][j]+hop[j][i])/2 for j in range(n)] for i in range(n)]
    def inertia(x):
        m = [row.copy() for row in s]
        for i in range(n):
            m[i][i] -= x
        count = 0
        for i in range(n):
            pivot = m[i][i]
            count += pivot < 0
            for j in range(i+1,n):
                factor = m[j][i]/pivot
                for k in range(i+1,j+1):
                    m[j][k] -= factor*m[k][i]
        return count
    bound = max(sum(map(abs, row)) for row in s) * 2
    widths = []
    for j in range(n):
        lo, hi = D(0), bound
        for _ in range(180):
            mid = (lo+hi)/2
            if inertia(mid) <= j:
                lo = mid
            else:
                hi = mid
        widths.append((lo+hi)/2)
    print(f"log10(min/max symmetric hopping eigenvalue) {(widths[0]/widths[-1]).log10():.8f}")
    mean = D(n-1)/2
    fit = sum((D(j)-mean)*w.log10() for j,w in enumerate(reversed(widths)))
    fit /= sum((D(j)-mean)**2 for j in range(n))
    print(f"symmetric hopping log10 slope per orbital {fit:.8f}")


if __name__ == "__main__":
    main()
