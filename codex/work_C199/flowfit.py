"""C199: two flow pictures for the qubit-chain size sequence, fixed before the converged L = 512 reading.
Inputs (converged, REPORT 9.178, 9.182, 9.186, 9.192): c_two(64->128) and c_two(128->256) at W = 1, 6, 12 for the two
keys with a chi = 1024 stage, and the K_s change per doubling from the pair estimate (REPORT 9.182).
(A) c = 2 with a correction from an irrelevant operator, c_two(L->2L) = 2 + a L^-p.
(B) a slow flow toward a pair liquid, c_two falling by a fixed amount per doubling (the 64->256 trend continued)."""
import numpy as np
data = {'p153c': {1: (2.196, 1.964), 6: (2.165, 1.998), 12: (2.234, 2.032)},
        'sym40_n316': {1: (2.312, 1.968), 6: (2.268, 1.999), 12: (2.356, 2.032)}}
print('(A) c = 2 plus a L^-p: deviation ratio per doubling and the predicted 256 -> 512 reading')
for k, d in data.items():
    for W, (a, b) in d.items():
        d1, d2 = a - 2, b - 2
        r = abs(d2 / d1); r_kt = 0.25
        # prediction: next deviation = d2 * r (fitted) or d2 * 1/4 (KT estimate below), bounded by |d2|
        print(f'  {k:11s} W={W:2d}: dev {d1:+.3f} -> {d2:+.3f}, ratio {r:.2f}; predicted 256->512 = {2 + d2 * max(r, r_kt):.3f} '
              f'(range {2 - abs(d2):.3f} to {2 + abs(d2):.3f})')
print('(B) linear decline per doubling continued')
for k, d in data.items():
    for W, (a, b) in d.items():
        print(f'  {k:11s} W={W:2d}: step {b - a:+.3f} per doubling -> 256->512 = {b + (b - a):.3f}, 512->1024 = {b + 2 * (b - a):.3f}')
print('K_s drift per doubling (pair estimate, REPORT 9.182): 64->128 -0.32 to -0.62, 128->256 -0.08 to -0.09')
for d1, d2 in ((-0.32, -0.09), (-0.62, -0.08)):
    r = d2 / d1
    print(f'  ratio {r:.2f} per doubling; remaining change if geometric: {d2 * r / (1 - r):+.3f}; '
          f'KT exponent: 2^(-4(K-1)) = r gives K = {1 - np.log2(r) / 4:.2f}')
