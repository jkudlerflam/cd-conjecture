"""Figures for notes_qubit_c2.tex.
Every plotted number is copied from REPORT.md of the cd-conjecture repository: the two-size readings from
sections 9.178 to 9.205, the pair gap and infinite-chain readings from 9.188, 9.201 and 9.204.
Run: python3 make_figs.py -> notes_qubit_fig1.pdf, notes_qubit_fig2.pdf, notes_qubit_fig3.pdf
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#8a8a85', '#1f1f1d', '#6b6b66'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['CMU Serif', 'Computer Modern Roman', 'DejaVu Serif'],
                     'mathtext.fontset': 'cm', 'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.labelcolor': INK, 'legend.frameon': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})


def ref(ax, y):
    ax.axhline(y, color=GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)


def clean(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ Fig. 1: the band with two valleys
fig, a = plt.subplots(figsize=(3.3, 1.85))
k = np.linspace(-np.pi, np.pi, 801)
t1, t2 = 1.0, 1.0                                   # |t1| = |t2|, phi = pi/2, units of |t2|
a.plot(k / np.pi, 2 * t1 * np.cos(k + np.pi / 2), color=GRAY, lw=1.0, ls=(0, (3, 2)), label=r'$t_1$ only')
a.plot(k / np.pi, 2 * t1 * np.cos(k + np.pi / 2) - 2 * t2 * np.cos(2 * k), color=BLUE, lw=1.5,
       label=r'$t_1$ and $t_2$')
a.legend(loc='upper right', handlelength=2.0)
k0 = np.arcsin(1 / 4)
for kv in (k0, np.pi - k0):
    a.plot(kv / np.pi, -2.25, 'o', ms=4.5, color=BLUE, mec='white', mew=0.8, zorder=3)
a.text(k0 / np.pi, -2.85, r'$k_0$', fontsize=8.5, color=INK, ha='center')
a.text((np.pi - k0) / np.pi, -2.85, r'$\pi-k_0$', fontsize=8.5, color=INK, ha='center')
a.set_xlim(-1, 1); a.set_ylim(-3.0, 4.4); a.set_xlabel(r'$k/\pi$'); a.set_ylabel(r'$\epsilon(k)/|t_2|$')
a.set_xticks([-1, -0.5, 0, 0.5, 1]); clean(a)
fig.savefig('notes_qubit_fig1.pdf')

# ------------------------------------------------------------------ Fig. 2: where the pair gap is (effective model)
n_g = np.array([1 / 48, 1 / 32, 1 / 24, 1 / 20, 1 / 16, 1 / 8, 5 / 32, 3 / 16])
d_g = np.array([0.0129, 0.0085, 0.0049, 0.0027, 0.0012, 0.0020, -0.0005, 0.0031])
e_g = np.array([0.0028, 0.0029, 0.0034, 0.0010, 0.0013, 0.0056, 0.0007, 0.0069])
n_c = np.array([1 / 32, 1 / 24, 1 / 20, 1 / 16, 3 / 40, 1 / 12, 3 / 32, 1 / 10, 1 / 8, 1 / 6, 3 / 16])
c_c = np.array([1.01, 1.04, 1.11, 1.49, 1.97, 1.98, 2.00, 2.00, 1.995, 1.993, 1.996])
fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 1.6), gridspec_kw=dict(wspace=0.34))
for ax in (a, b):
    ax.axvspan(1 / 16, 3 / 40, color='#e8eef8', zorder=0, lw=0)
    ax.axvline(3 / 16, color=ORANGE, lw=0.8, ls=(0, (1, 1.5)), zorder=0)
    ax.set_xlim(-0.004, 0.2); ax.set_xlabel(r'density $n$'); clean(ax)
ref(a, 0)
a.errorbar(n_g, d_g, yerr=e_g, fmt='o', ms=4, color=BLUE, mec='white', mew=0.8, elinewidth=0.8, capsize=0)
a.plot([0], [0.0195], marker='*', ms=8, color=INK, ls='none', clip_on=False)
a.annotate(r'$E_b$ (two-body)', (0.002, 0.0195), xytext=(0.014, 0.0192), fontsize=7.5, color=INK, va='center')
a.set_ylabel(r'$\Delta_{\rm pb}(L\to\infty)/|t_2|$'); a.set_ylim(-0.009, 0.024)
a.text(3 / 16, 0.0222, r'$n=3/16$', fontsize=7.5, color=ORANGE, ha='center')
a.text(0.0688, 0.0222, r'$n_c$', fontsize=8, color=BLUE, ha='center')
a.text(-0.24, 1.0, '(a)', transform=a.transAxes, fontsize=9)
ref(b, 2); ref(b, 1)
b.plot(n_c, c_c, color=BLUE, lw=1.5, marker='o', ms=4.5, mec='white', mew=0.8)
b.set_ylabel(r'$c$, infinite chain'); b.set_ylim(0.8, 2.25); b.set_yticks([1.0, 1.5, 2.0])
b.text(-0.24, 1.0, '(b)', transform=b.transAxes, fontsize=9)
fig.savefig('notes_qubit_fig2.pdf')

# ------------------------------------------------------------------ Fig. 3: the qubit chain
fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 1.9), gridspec_kw=dict(wspace=0.32, width_ratios=[1.15, 1]))
L2 = np.array([128, 256, 512])
qp = [2.165, 2.032, 1.954]          # phi = 1.53 point; L = 256 at chi = 1536, L = 512 at chi = 2048
qs = [2.268, 2.025, 1.983]          # symmetric point
qp_ext, qs_ext = 2.021, 2.041       # 256 -> 512, linear extrapolation in the discarded weight
eff = [1.817, 1.877, 1.93]          # effective model, open chains, n = 3/16
c2 = [1.917, 1.898, 1.952]          # exact c = 2 control
c1 = [0.940, 0.973, 0.988]          # exact c = 1 control
ref(a, 2); ref(a, 1)
a.plot(L2, c2, color=GRAY, lw=1.0, marker='s', ms=4, mfc='white', mec=GRAY, label=r'exact $c=2$')
a.plot(L2, c1, color=GRAY, lw=1.0, marker='D', ms=3.5, mfc='white', mec=GRAY, label=r'exact $c=1$')
a.plot(L2, eff, color=AQUA, lw=1.5, marker='^', ms=5, mec='white', mew=0.8, label='effective model')
a.plot(L2, qs, color=BLUE, lw=1.5, marker='o', ms=5, mec='white', mew=0.8, label=r'qubit chain, $\varphi=\pi/2$')
a.plot(L2, qp, color=ORANGE, lw=1.5, marker='s', ms=4.5, mec='white', mew=0.8, label=r'qubit chain, $\varphi=1.53$')
a.plot([512 * 1.05], [qs_ext], marker='o', ms=4.5, mfc='white', mec=BLUE, mew=1.0, ls='none')
a.plot([512 * 0.95], [qp_ext], marker='s', ms=4, mfc='white', mec=ORANGE, mew=1.0, ls='none')
a.text(512, 2.13, 'extrap.', fontsize=7, color=MUTED, ha='center')
a.set_xscale('log', base=2); a.set_xticks(L2)
a.set_xticklabels([r'$64\to128$', r'$128\to256$', r'$256\to512$'])
a.set_xlim(100, 680); a.set_ylim(0.6, 2.5); a.set_ylabel(r'$c(L\to2L)$, $W=6$')
a.set_yticks([1.0, 1.5, 2.0, 2.5]); a.minorticks_off()
a.legend(loc='center', bbox_to_anchor=(0.5, 0.37), handlelength=2.0, ncol=2, columnspacing=0.8, fontsize=6.8)
a.text(-0.2, 1.0, '(a)', transform=a.transAxes, fontsize=9); clean(a)

data = [(r'$\varphi=1.53$', ORANGE, 's', [8.9e-10, 4.4e-10, 1.5e-10, 6.0e-11], [1.202, 1.550, 1.851, 1.954], 2.021),
        (r'$\varphi=\pi/2$', BLUE, 'o', [7.9e-10, 3.8e-10, 1.3e-10, 5.0e-11], [1.290, 1.622, 1.893, 1.983], 2.041)]
ref(b, 2)
for lab, col, mk, e, c, ext in data:
    e = np.array(e) * 1e10
    b.plot(e, c, color=col, lw=1.5, marker=mk, ms=4.5, mec='white', mew=0.8, label=lab)
    b.plot([0], [ext], marker=mk, ms=4.5, mfc='white', mec=col, mew=1.0, ls='none', clip_on=False)
    b.plot([0, e[-1]], [ext, c[-1]], color=col, lw=0.8, ls=(0, (2, 2)))
lab_pts = [(r'$\chi=768$', (8.9, 1.202), (8.4, 1.52)), ('1024', (4.4, 1.550), (5.0, 1.30)),
           ('1536', (1.5, 1.851), (3.0, 1.93)), ('2048', (0.6, 1.954), (1.9, 2.13))]
for txt, xy, xyt in lab_pts:
    b.annotate(txt, xy, xytext=xyt, fontsize=7, color=MUTED, ha='center', va='center',
               arrowprops=dict(arrowstyle='-', lw=0.4, color=GRAY, shrinkA=1, shrinkB=3))
b.set_xlabel(r'discarded weight at $L=512$ ($\times10^{-10}$)'); b.set_ylabel(r'$c(256\to512)$, $W=6$')
b.set_xlim(-0.3, 9.8); b.set_ylim(1.0, 2.2); b.legend(loc='lower left')
b.text(-0.2, 1.0, '(b)', transform=b.transAxes, fontsize=9); clean(b)
fig.savefig('notes_qubit_fig3.pdf')
print('ok')
