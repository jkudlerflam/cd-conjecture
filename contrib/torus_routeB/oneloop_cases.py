#!/usr/bin/env python3
"""Task 1, one-loop part: project the lattice interactions of models.py onto u(4) current-current
couplings at the free Fermi points, keep the resonant ones, and integrate the one-loop flow.

Checks first: Hubbard calibration of the sign (2 flavours, equal k: U > 0 bounded, U < 0 runaway),
Hermiticity and torus symmetry of the lattice terms, and the selection rule itself (the resonant
couplings of the torus cases contain no pair tunnelling).  Then, for every case and both signs,
the seed-0 flow at several g, and counts of runaway flows over 100 random draws.

Run: python3 oneloop_cases.py   (about a minute)  ->  oneloop_cases_stdout.txt
Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import numpy as np
import fermions as F
import tcases as M

T4, f4 = F.uN_basis(4)
T2, f2 = F.uN_basis(2)

print('margins (all 2- and 4-field momentum sums), K_DIST:',
      {k: round(v / np.pi, 4) for k, v in F.margins(M.K_DIST).items()}, '(units of pi)')

for U in (0.5, -0.5):
    k = np.array([0.3 * np.pi] * 2)
    G = F.project_uN([(U, ((1, 0, 0), (0, 0, 0), (1, 0, 1), (0, 0, 1)))], k)
    fate, l, Gf, _ = F.one_loop(G, F.resonance_mask(k), T=T2, f=f2, lmax=100)
    print(f'calibration: Hubbard U = {U:+.1f}, k = 0.3 pi: {fate} (l = {l:.2f}); expected '
          f'{"bounded (spin coupling marginally irrelevant)" if U > 0 else "runaway (spin gap)"}')

for case in ('i', 'ii', 'iii'):
    for s in (1, -1):
        k, monos, inter = M.case_monos(case, s, 0.3, 0)
        H = F.terms_matrix(monos, 2)
        Nt = F.total_number(2)
        print(f'lattice check case {case} s={s:+d}: |H-H^dag| = {np.linalg.norm(H - H.conj().T):.1e}, '
              f'|[H,N]| = {np.linalg.norm(H @ Nt - Nt @ H):.1e}, torus |[H,T]| = {F.check_torus(H, 2):.1e}')

print('\nseed 0: resonant couplings (Frobenius norms by type) and the one-loop fate; '
      'l = RG time at which max|G| exceeds max(1, 5 max|G0|), v = 1')
for g in (0.1, 0.2, 0.3):
    for case in ('i', 'ii', 'iii'):
        for s in (1, -1):
            k, monos, inter = M.case_monos(case, s, g, 0)
            G, other = F.project_uN(inter, k, report=True)
            cls = F.classify(G)
            fate, l, Gf, _ = F.one_loop(G, F.resonance_mask(k), T=T4, f=f4, lmax=60)
            print(f'  g={g} case {case:3s} s={s:+d}: ' + ' '.join(f'{kk}={v:.3f}' for kk, v in cls.items())
                  + f' | non-current resonant terms {len(other)} | {fate} at l = {l:.2f}')

print('\ncounts over seeds 0..99 at g = 0.1 (runaway / bounded within l = 40)')
for case in ('i', 'ii', 'iii'):
    for s in (1, -1):
        nr, ls = 0, []
        for seed in range(100):
            k, monos, inter = M.case_monos(case, s, 0.1, seed)
            G = F.project_uN(inter, k)
            fate, l, _, _ = F.one_loop(G, F.resonance_mask(k), T=T4, f=f4, lmax=40)
            if fate == 'runaway':
                nr += 1; ls.append(l)
        print(f'  case {case:3s} s={s:+d}: runaway in {nr} of 100' +
              (f', l* median {np.median(ls):.2f} (range {min(ls):.2f} to {max(ls):.2f})' if ls else ''))

print('\nscale-dependent detuning (all torus-allowed / U(1)-allowed couplings kept until their momentum '
      'exceeds the running cutoff e^{-l}), seed 0:')
for g in (0.1, 0.2, 0.3, 0.5):
    for case in ('i', 'ii'):
        for s in (1, -1):
            k, monos, inter = M.case_monos(case, s, g, 0)
            G = F.project_uN(inter, k, keep_all=True)
            fate, l, Gf, _ = F.one_loop(G, np.ones((4, 4, 4, 4)), T=T4, f=f4, lmax=60, Pmat=F.momentum_matrix(k))
            nc = np.abs(Gf - np.einsum('ab,cd->abcd', np.eye(4), np.eye(4)) * Gf).max()
            print(f'  g={g} case {case:3s} s={s:+d}: {fate} at l = {l:.2f}; final max non-Cartan |G| = {nc:.2e}')

# ---------------------------------------------------------------- so(8): umklapp-type couplings included
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'redteam_construction'))
import route_b_second_order as RB
f8 = RB.so_structure(8)
U8, U8i, q8 = F.so_charges(4)
for U in (0.5, -0.5):
    for kk in (0.3, 0.5):
        k2 = np.array([kk * np.pi] * 2)
        U4, U4i, q4 = F.so_charges(2)
        G = F.project_so([(U, ((1, 0, 0), (0, 0, 0), (1, 0, 1), (0, 0, 1)))], k2, N=2)
        fate, l, _ = F.so_flow(G.real, RB.so_structure(4), U4, U4i, F.so_momenta(q4, k2), lmax=100)
        print(f'so(4) calibration: Hubbard U = {U:+.1f}, k = {kk} pi: {fate} (l = {l:.2f})')
print('\nso(8) flows with every coupling (umklapp included) kept until its momentum exceeds e^{-l}, seed 0:')
for g in (0.1, 0.2, 0.3):
    for case in ('i', 'ii', 'iii'):
        for s in (1, -1):
            k, monos, inter = M.case_monos(case, s, g, 0)
            G = F.project_so(inter, k)
            assert np.abs(G.imag).max() < 1e-10
            P = F.so_momenta(q8, k)
            Gt = U8i @ G.real @ U8i.T
            umk = np.abs(Gt[(np.abs(q8.sum(1))[:, None] == 2) & (P < np.pi)]).max()
            fate, l, _ = F.so_flow(G.real, f8, U8, U8i, P, lmax=60)
            print(f'  g={g} case {case:3s} s={s:+d}: {fate} at l = {l:.2f} (largest umklapp-type amplitude {umk:.3f})')
print('\nso(8) counts over seeds 0..99 (runaway within l = 40)')
for g in (0.1, 0.2):
    for case in ('i', 'ii', 'iii'):
        for s in (1, -1):
            nr = 0
            for seed in range(100):
                k, monos, inter = M.case_monos(case, s, g, seed)
                G = F.project_so(inter, k)
                fate, l, _ = F.so_flow(G.real, f8, U8, U8i, F.so_momenta(q8, k), lmax=40)
                nr += fate == 'runaway'
            print(f'  g={g} case {case:3s} s={s:+d}: runaway in {nr} of 100')

# ---------------------------------------------------------------- locking criterion for case (ii)
print('\ncase (ii), seed 0: if flavours a, b had equal momenta, the so(8) flow would run away at l*_ab;'
      ' locking (commensurate-incommensurate transition of the conserved n_a - n_b) is expected when the'
      ' resulting gap ~ v e^{-l*_ab} exceeds the separation v |k_a - k_b|.  ratio = e^{-l*}/|k_a - k_b|:')
for g in (0.1, 0.2, 0.3):
    for s in (1, -1):
        k, monos, inter = M.case_monos('ii', s, g, 0)
        row = []
        for (a, b) in ((1, 2), (0, 1), (2, 3), (0, 2), (1, 3), (0, 3)):
            kk = k.copy(); kk[a] = kk[b] = (k[a] + k[b]) / 2
            fate, l, _ = F.so_flow(F.project_so(inter, kk).real, f8, U8, U8i, F.so_momenta(q8, kk), lmax=40)
            ratio = np.exp(-l) / abs(k[a] - k[b]) if fate == 'runaway' else 0.0
            row.append(f'({a},{b}) l*={l:.1f} ratio={ratio:.2f}')
        print(f'  g={g} s={s:+d}: ' + '; '.join(row))

# parity-only: the quadratic pairing f_{0a} f_{1a} contains psi_{aR} psi_{aL} at zero momentum
k = M.K_DIST
amp = [abs(np.exp(-1j * k[a]) - np.exp(1j * k[a])) for a in range(4)]
print('\nparity-only (case ip): f_{j,a} f_{j+1,a} -> (e^{-ik} - e^{ik}) psi_{aR} psi_{aL} + ..., amplitude',
      np.round(amp, 3), 'at zero momentum: a relevant (dimension 1) pairing mass for every flavour')
