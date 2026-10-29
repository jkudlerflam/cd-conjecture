"""Parameter packing for the adversarial searches: real vector -> (A0, B)."""
import numpy as np


def unpack_nc(x, m):
    A0r = x[:m * m].reshape(m, m); A0i = x[m * m:2 * m * m].reshape(m, m)
    A0 = (A0r + A0r.T) / 2 + 1j * (A0i - A0i.T) / 2
    Br = x[2 * m * m:3 * m * m].reshape(m, m); Bi = x[3 * m * m:].reshape(m, m)
    return A0, Br + 1j * Bi


def unpack_mj(x, n):
    A0r = x[:n * n].reshape(n, n)
    A0 = (A0r - A0r.T) / 2
    B = x[n * n:].reshape(n, n)
    return A0, B
