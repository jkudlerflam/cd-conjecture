#!/usr/bin/env python3
"""Exact data for the XY + DM_z + h_z control (free fermions, no U(1), incommensurate, c = 1).

H = sum_j [Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + D (X_j Y_{j+1} - Y_j X_{j+1})] + h sum_j Z_j.
Jordan-Wigner gives xi_k = 2(Jx+Jy) cos k - 4 D sin k + 2h and pairing 2i(Jx-Jy) sin k on (k,-k).
Prints the ground-state energy per site and the Fermi points (zeros of the quasiparticle band).
"""
import sys
import numpy as np


def pair_energy(k, Jx, Jy, D, h):
    xi = lambda q: 2 * (Jx + Jy) * np.cos(q) - 4 * D * np.sin(q) + 2 * h
    a, b = xi(k), xi(-k)
    s = 0.5 * (a + b)
    Del = 2 * (Jx - Jy) * np.sin(k)
    even = s - np.sqrt(s ** 2 + Del ** 2)
    return np.minimum(np.minimum(even, a), b), xi


def e_inf(Jx, Jy, D, h, N=2_000_001):
    k = np.linspace(0, np.pi, N)
    E, _ = pair_energy(k, Jx, Jy, D, h)
    # sum over pairs (k,-k), k in (0, pi): (1/2pi) int_0^pi E dk per site, minus h from Z = 2n - 1
    return np.trapezoid(E, k) / (2 * np.pi) - h


def fermi_points(Jx, Jy, D, h, N=200001):
    k = np.linspace(-np.pi, np.pi, N)
    xi = lambda q: 2 * (Jx + Jy) * np.cos(q) - 4 * D * np.sin(q) + 2 * h
    s = 0.5 * (xi(k) + xi(-k)); a = 0.5 * (xi(k) - xi(-k))
    Ek = a + np.sqrt(s ** 2 + (2 * (Jx - Jy) * np.sin(k)) ** 2)
    z = np.where(np.sign(Ek[1:]) != np.sign(Ek[:-1]))[0]
    return k[z]


if __name__ == '__main__':
    Jx, Jy, D, h = map(float, sys.argv[1:5]) if len(sys.argv) > 4 else (1.0, 0.6, 0.5, 0.6)
    print('e_inf', e_inf(Jx, Jy, D, h), 'fermi points / pi', fermi_points(Jx, Jy, D, h) / np.pi)
