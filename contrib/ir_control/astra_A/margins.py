"""Closed-book numerical checks for RESULT.md. Run with python3 margins.py.

Only numpy and the Python standard library are used. Circle arithmetic requests
longdouble, including pi computed as arccos(-1). On ARM macOS longdouble is
float64; the output records the actual precision and is not interval-certified.
The search exhausts distinct four-index sets and all signs up to overall sign.
"""
import itertools
import json
import numpy as np

LD = np.longdouble
PI = np.arccos(LD(-1))
TAU = 2 * PI


def primes(count):
    out = []
    k = 2
    while len(out) < count:
        if all(k % p for p in out if p * p <= k):
            out.append(k)
        k += 1
    return out


def circ(x):
    return np.abs((x + PI) % TAU - PI)


def margins(points, labels):
    points = np.asarray(points, dtype=LD)
    size = len(points)
    best2 = (LD(10), None)
    for i in range(size):
        val = circ(2 * points[i])
        if val < best2[0]:
            best2 = (val, [(2, i)])
        for j in range(i + 1, size):
            for s in (-1, 1):
                val = circ(points[i] + s * points[j])
                if val < best2[0]:
                    best2 = (val, [(1, i), (s, j)])

    suffix = {}
    for start in range(2, size - 1):
        ii, jj = np.triu_indices(size - start, 1)
        suffix[start] = (ii + start, jj + start)
    best4 = (LD(10), None)
    signs = list(itertools.product((-1, 1), repeat=3))
    for i in range(size - 3):
        for j in range(i + 1, size - 2):
            kk, ll = suffix[j + 1]
            for sj, sk, sl in signs:
                values = circ(points[i] + sj * points[j]
                              + sk * points[kk] + sl * points[ll])
                pos = int(np.argmin(values))
                if values[pos] < best4[0]:
                    best4 = (values[pos], [(1, i), (sj, j),
                                          (sk, int(kk[pos])),
                                          (sl, int(ll[pos]))])

    def describe(best):
        distance, terms = best
        signed_sum = sum(LD(s) * points[i] for s, i in terms)
        winding = int(np.rint(signed_sum / TAU))
        residual = signed_sum - winding * TAU
        return dict(distance=str(distance),
                    terms=[dict(coefficient=s, index=i + 1, label=labels[i])
                           for s, i in terms],
                    winding=winding, signed_residual=str(residual))
    return dict(delta2=describe(best2), delta4=describe(best4))


def designed(n=40, irrational_jitter=True):
    size = 2 * n
    q = next(p for p in primes(size) if p >= size)
    indices = np.arange(size, dtype=np.int64)
    a = 2 * q * indices + indices ** 2 % q
    d = int(np.min(np.diff(a)))
    b = (d + 1) // 2
    period = 2 * int(a[-1]) + 4 * b + 2
    points = TAU * (a.astype(LD) + LD(b) + LD(1) / 6) / period
    jitter = LD('1e-10') * np.sqrt(np.array(primes(size), dtype=LD))
    if irrational_jitter:
        points += jitter
    # Pair point a with point n+a: the latter is the right Fermi point.
    k = (points[n:] - points[:n]) / 2
    theta = (points[n:] + points[:n]) / 2
    velocity = 2 * np.sin(k)  # t_a = 1
    labels = [f"a={i+1},-" for i in range(n)]
    labels += [f"a={i+1},+" for i in range(n)]
    data = dict(q=q, d=d, B=b, T=period, a_max=int(a[-1]),
                min_velocity=str(np.min(velocity)),
                scaffold_delta2=str(TAU * d / period),
                scaffold_delta4=str(TAU / (3 * period)),
                delta2_lower_bound=str(TAU * d / period
                                       - (2 * np.max(jitter) if irrational_jitter else 0)),
                delta4_lower_bound=str(TAU / (3 * period)
                                       - (4 * np.max(jitter) if irrational_jitter else 0)),
                irrational_jitter=irrational_jitter,
                max_jitter=str(np.max(jitter)) if irrational_jitter else '0',
                theta=[str(x) for x in theta], k=[str(x) for x in k])
    assert np.all(velocity >= 1)
    # Exact certificates for the rational scaffold's two attainers at n=40.
    if n == 40:
        numerator = 6 * (a + b) + 1
        assert abs(int(numerator[49] - numerator[50])) == 6 * d
        assert int(numerator[0] - numerator[30] - numerator[51]
                   - numerator[77]) == -6 * period - 2
    return points, labels, data


def sanity_checks():
    # A four-index exhaustive check against a separate direct implementation.
    sample = np.array([.12, .71, 1.37, 2.03, 2.61, 2.98], dtype=LD)
    result = margins(sample, [str(i) for i in range(len(sample))])
    reference = min(circ(sum(s * sample[i] for s, i in zip(sign, ids)))
                    for ids in itertools.combinations(range(len(sample)), 4)
                    for sign in itertools.product((-1, 1), repeat=4))
    assert abs(LD(result['delta4']['distance']) - reference) < LD('1e-17')
    # The opposite-chirality example loses its two roots at h = v delta / 2.
    v, delta = 2., .1
    for h in (.03, .15):
        discriminant = (delta / 2) ** 2 - (h / v) ** 2
        assert (discriminant > 0) == (h < v * delta / 2)


if __name__ == '__main__':
    sanity_checks()
    prime_list = np.asarray(primes(80), dtype=LD)
    theta = np.sqrt(prime_list[:40]) / 100
    k = PI / 2 + np.sqrt(prime_list[40:]) / 100
    points = np.empty(80, dtype=LD)
    points[0::2], points[1::2] = theta + k, theta - k
    labels = [f'a={a},{s}' for a in range(1, 41) for s in ('+', '-')]
    example = margins(points, labels)
    example['min_velocity_at_t_1'] = str(np.min(2 * np.sin(k)))
    pp, ll, meta = designed()
    optimized = margins(pp, ll)
    optimized.update(meta)
    output = dict(dtype=str(np.dtype(LD)),
                  precision_bits=int(np.finfo(LD).nmant),
                  prime_example=example, designed_example=optimized,
                  sanity_checks='passed')
    with open('numerical_results.json', 'w') as stream:
        json.dump(output, stream, indent=2)
    for title, result in [('prime example', example), ('designed example', optimized)]:
        print(title)
        print(json.dumps({k: v for k, v in result.items()
                          if k not in ('theta', 'k')}, indent=2))
