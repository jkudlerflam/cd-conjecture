"""Figures for notes_unbounded_c.tex.  Every number is read from a file in
contrib/markerless_d9_d10/ of this repository; nothing is typed in by hand except the schematic.
Run:  python3 make_figs.py   ->  figU1.pdf, figU2.pdf, figU3.pdf   (needs numpy and matplotlib)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'markerless_d9_d10')
BLUE, ORANGE, AQUA, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#8a8a85', '#1f1f1d', '#6b6b66'
PALE_BLUE, PALE_ORANGE, PALE_GRAY = '#dbe9f9', '#fbe3d8', '#ececea'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['CMU Serif', 'Computer Modern Roman', 'DejaVu Serif'],
                     'mathtext.fontset': 'cm', 'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 6.5,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.labelcolor': INK, 'legend.frameon': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
W = 3.375  # PRL column width, inches


def clean(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ Fig. 1: soliton lattice, pair creation, kink states
fig = plt.figure(figsize=(W, 3.75))
a = fig.add_axes([0.0, 0.46, 1.0, 0.54]); b = fig.add_axes([0.0, 0.0, 1.0, 0.41])
for ax in (a, b):
    ax.set_axis_off()
FILL = {'L': PALE_BLUE, 'R': PALE_ORANGE, 'K': '#ffffff', 'C': '#ffffff'}
EDGE = {'L': BLUE, 'R': ORANGE, 'K': INK, 'C': AQUA}


def site(ax, x, y, kind, occ=None, label=None, w=0.84, h=0.62):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.0,rounding_size=0.10',
                                fc=FILL[kind], ec=EDGE[kind], lw=1.1 if kind in 'KC' else 0.7))
    if label is not None:
        ax.text(x, y - 0.01, label, ha='center', va='center', fontsize=6.6, color=INK)
    elif occ is not None:
        ax.plot(x, y, 'o', ms=4.2, mfc=INK if occ else 'white', mec=INK, mew=0.7, zorder=3)


def row(ax, y, cells, tag):
    """cells: list of two lists of (kind, occupation or label)"""
    x = 0.0
    for c, cell in enumerate(cells):
        for kind, v in cell:
            site(ax, x, y, kind, occ=v if kind in 'LR' else None, label=v if kind in 'KC' else None); x += 1.0
        if c == 0:
            ax.plot([x - 0.5, x - 0.5], [y - 0.47, y + 0.47], color=MUTED, lw=0.7, ls=(0, (2, 2)))
    ax.text(-0.75, y, tag, ha='right', va='center', fontsize=6.6, color=MUTED)


K = ('K', r'$K$')
row(a, 2.3, [[('L', 1), ('L', 0), K, ('R', 1), ('R', 0)], [('L', 0), K, ('R', 1), ('R', 1), ('R', 0)]], '(i)')
row(a, 1.15, [[('L', 1), ('L', 0), ('L', 1), ('L', 0), K], [K, ('R', 0), ('R', 1), ('R', 1), ('R', 0)]], '(ii)')
row(a, 0.0, [[('L', 1), ('L', 0), ('L', 1), ('L', 0), ('C', r'$\tilde K_{\!-}$')],
             [('C', r'$K_{\!+}$'), ('R', 0), ('R', 1), ('R', 1), ('R', 0)]], '(iii)')
for x0 in (2.0, 6.0):
    a.add_patch(FancyArrowPatch((x0 - 0.42, 2.3 + 0.52), (x0 + 0.42, 2.3 + 0.52), arrowstyle='<->', mutation_scale=5,
                                lw=0.6, color=INK))
a.text(2.0, 3.02, 'kink', ha='center', va='bottom', fontsize=6.2, color=INK)
a.text(0.5, 3.02, r'$L$ domain', ha='center', va='bottom', fontsize=6.2, color=BLUE)
a.text(3.5, 3.02, r'$R$ domain', ha='center', va='bottom', fontsize=6.2, color=ORANGE)
a.text(4.5, -0.52, r'antikink (immobile)', ha='center', va='top', fontsize=6.2, color=MUTED)
a.plot(7.55, 3.17, 'o', ms=4.2, mfc=INK, mec=INK, mew=0.7); a.text(7.75, 3.17, 'occupied', fontsize=6.2, va='center', color=INK)
a.plot(7.55, 2.86, 'o', ms=4.2, mfc='white', mec=INK, mew=0.7); a.text(7.75, 2.86, 'empty', fontsize=6.2, va='center', color=INK)
a.add_patch(FancyArrowPatch((4.5, 1.15 - 0.38), (4.5, 0.0 + 0.38), arrowstyle='<->', mutation_scale=6, lw=0.7, color=AQUA))
a.text(4.72, 0.575, r'pair creation $\lambda$', ha='left', va='center', fontsize=6.4, color=INK)
a.set_xlim(-1.6, 9.7); a.set_ylim(-0.95, 3.45)
a.text(-1.55, 3.3, '(a)', fontsize=8, color=INK, va='top')

# (b) states of a charged kink inside its cell: dark branch on top, bright branch below
n = 5
xs = np.arange(n + 1) * 1.25
yX, yD = 1.25, 0.0
b.plot(xs, np.full(n + 1, yX), color=GRAY, lw=1.0, zorder=1)
b.plot(xs, np.full(n + 1, yD), color=AQUA, lw=1.4, zorder=1)
b.plot([xs[0], xs[0]], [yD, yX], color=INK, lw=1.0, zorder=1)
b.plot([xs[-1], xs[-1]], [yD, yX], color=ORANGE, lw=0.9, ls=(0, (2, 2)), zorder=1)
for j, x in enumerate(xs):
    b.plot(x, yX, 'o', ms=5.5, mfc='white', mec=GRAY, mew=1.0, zorder=3)
    b.plot(x, yD, 'o', ms=5.5, mfc='white' if j == n else AQUA, mec=AQUA, mew=1.0, zorder=3)
    if j < n:
        b.text(x, yD - 0.36, r'$g$', ha='center', va='top', fontsize=6.4, color=INK)
    if j in (0, n):
        b.text(x, yX + 0.26, r'$0$' if j == 0 else r'$n$', ha='center', va='bottom', fontsize=6.2, color=MUTED)
b.text(xs[0] - 0.32, yX, 'dark', ha='right', va='center', fontsize=6.6, color=MUTED)
b.text(xs[0] - 0.32, yD, 'bright', ha='right', va='center', fontsize=6.6, color=INK)
b.text(xs[2] + 0.62, yX + 0.26, r'position of the kink in its cell', ha='center', va='bottom', fontsize=6.2, color=MUTED)
b.text(xs[2] + 0.62, yX - 0.2, r'energy $D$, no coupling to the bosons', ha='center', va='top', fontsize=6.2, color=MUTED)
b.text(xs[2] + 0.62, yD + 0.2, r'exchanges charge with the sites', ha='center', va='bottom', fontsize=6.2, color=INK)
b.text(xs[-1] + 0.12, (yX + yD) / 2, r'$t_d\,e^{i\psi}$', ha='left', va='center', fontsize=6.6, color=INK)
b.add_patch(FancyArrowPatch((xs[-1] + 1.55, yX + 0.55), (xs[-1] + 0.14, yX + 0.12), arrowstyle='->', mutation_scale=6,
                            lw=0.7, color=INK))
b.text(xs[-1] + 1.6, yX + 0.58, 'kink of the\nleft cell', ha='center', va='bottom', fontsize=6.2, color=INK, linespacing=0.95)
b.add_patch(FancyArrowPatch((xs[0] - 1.45, yD - 0.6), (xs[0] - 0.14, yD - 0.12), arrowstyle='->', mutation_scale=6,
                            lw=0.7, color=INK))
b.text(xs[0] - 1.5, yD - 0.62, 'kink of the\nright cell', ha='center', va='top', fontsize=6.2, color=INK, linespacing=0.95)
b.set_xlim(-2.9, xs[-1] + 2.6); b.set_ylim(-1.35, 2.5)
b.text(-2.85, 2.45, '(b)', fontsize=8, color=INK, va='top')
fig.savefig(os.path.join(HERE, 'figU1.pdf')); plt.close(fig)

# ------------------------------------------------------------------ Fig. 2: bands of the certified n = 9 member
z = np.load(os.path.join(SRC, 'certificate', 'full_loop_whitened_n9.npz'), allow_pickle=True)
T, mu, field = z['T'], z['mu'], z['field']                      # whitened: Herm part of T is -1/2
bands = json.load(open(os.path.join(SRC, 'certificate', 'full_loop_bands_n9.json')))
roots = np.array(sorted(bands['roots']))


def spectrum(k):
    return np.array([np.linalg.eigvalsh(mu + field + T * np.exp(1j * x) + T.T * np.exp(-1j * x)) for x in k])


fig = plt.figure(figsize=(W, 3.5))
top = fig.add_axes([0.14, 0.60, 0.84, 0.38])
bl = fig.add_axes([0.14, 0.11, 0.37, 0.36]); br = fig.add_axes([0.61, 0.11, 0.37, 0.36])
k = np.linspace(-np.pi, np.pi, 1601)
ev = spectrum(k)
for j in range(ev.shape[1]):
    top.plot(k / np.pi, ev[:, j], color=BLUE, lw=0.7)
top.axhline(0, color=GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
top.set_xlim(-1, 1); top.set_ylim(-1.25, 1.35); top.set_xticks([-1, -0.5, 0, 0.5, 1])
top.set_xlabel(r'$k/\pi$', labelpad=1); top.set_ylabel(r'$\tilde\epsilon_a(k)$')
top.text(0.0, 1.05, r'$n=9$', fontsize=6.6, color=INK, ha='center')
top.text(-0.19, 1.02, '(a)', transform=top.transAxes, fontsize=8)
clean(top)
for ax, sel, lab, ttl in ((bl, roots < 0, '(b)', '9 left-movers'), (br, roots > 0, '(c)', '9 right-movers')):
    r = roots[sel]
    lo, hi = r.min() - 0.035, r.max() + 0.035
    kk = np.linspace(lo, hi, 801); e2 = spectrum(kk)
    for j in range(e2.shape[1]):
        ax.plot(kk / np.pi, e2[:, j], color=BLUE, lw=0.7)
    ax.axhline(0, color=GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.plot(r / np.pi, np.zeros_like(r), 'o', ms=3.4, mfc=ORANGE, mec='white', mew=0.5, zorder=4)
    ax.set_xlim(lo / np.pi, hi / np.pi); ax.set_ylim(-0.085, 0.085)
    ax.set_xlabel(r'$k/\pi$', labelpad=1); clean(ax)
    ax.locator_params(axis='x', nbins=4)
    ax.text(-0.02, 1.04, lab, transform=ax.transAxes, fontsize=8, ha='right')
    ax.text(0.5, 1.04, ttl, transform=ax.transAxes, fontsize=6.6, color=INK, ha='center')
    top.axvspan(lo / np.pi, hi / np.pi, color=PALE_GRAY, zorder=-1, lw=0)
bl.set_ylabel(r'$\tilde\epsilon_a(k)$')
fig.savefig(os.path.join(HERE, 'figU2.pdf')); plt.close(fig)

# ------------------------------------------------------------------ Fig. 3: c against d, and the velocity hierarchy
ns = [9, 10, 11, 12, 14, 16, 20]
cert = {m: json.load(open(os.path.join(SRC, 'certificate', f'root_certificate_n{m}.json'))) for m in ns}
assert all(cert[m]['status'] == 'proved_integer_interval_certificate' and cert[m]['root_count'] == 2 * m for m in ns)
vel = json.load(open(os.path.join(SRC, 'data', 'hopping_eigenvalues_n9.json')))['log10_eigenvalues']
vel = np.array(vel) - vel[0]
fig, (a, b) = plt.subplots(1, 2, figsize=(W, 1.95), gridspec_kw=dict(wspace=0.5, left=0.12, right=0.985, bottom=0.2, top=0.9))
dd = np.arange(2, 15)
a.plot(dd, dd - 1, color=GRAY, lw=1.0)
a.fill_between(dd, dd - 1, 23, color=PALE_GRAY, lw=0, zorder=0)
a.plot(dd, dd - 1, 'o', ms=2.6, color=GRAY, zorder=2)
a.plot(np.full(len(ns), 10), [m + 1 for m in ns], 'o', ms=3.6, mfc=BLUE, mec='white', mew=0.5, zorder=4)
a.plot(np.full(len(ns), 9), ns, 'o', ms=3.6, mfc='white', mec=BLUE, mew=0.9, zorder=3)
a.text(11.0, 19.6, r'$c=n{+}1$', fontsize=6.4, color=BLUE, ha='left')
a.text(11.0, 17.6, r'$d=10$', fontsize=6.4, color=BLUE, ha='left')
a.text(8.4, 13.3, r'$c=n$,', fontsize=6.4, color=INK, ha='right')
a.text(8.4, 11.5, r'$d=9$', fontsize=6.4, color=INK, ha='right')
a.text(12.6, 6.6, r'$c=d-1$', fontsize=6.4, color=MUTED, ha='center', rotation=38)
a.text(3.0, 20.0, r'$c>d-1$', fontsize=6.4, color=MUTED, ha='left', va='top')
a.set_xlim(2, 14); a.set_ylim(0, 23); a.set_xticks([2, 6, 10, 14]); a.set_yticks([0, 5, 10, 15, 20])
a.set_xlabel(r'local dimension $d$', labelpad=1); a.set_ylabel(r'central charge $c$', labelpad=1)
a.text(-0.3, 1.04, '(a)', transform=a.transAxes, fontsize=8); clean(a)
aa = np.arange(1, len(vel) + 1)
slope = np.polyfit(aa, vel, 1)
b.plot(aa, np.polyval(slope, aa), color=GRAY, lw=0.8, ls=(0, (3, 2)), zorder=1)
b.plot(aa, vel, 'o', ms=3.6, mfc=AQUA, mec='white', mew=0.5, zorder=3)
b.text(8.9, -3.0, r'$v_a\propto r^{a}$', fontsize=6.6, color=INK, ha='right')
b.text(8.9, -6.2, r'$r\simeq10^{%.2f}$' % slope[0], fontsize=6.6, color=INK, ha='right')
b.set_xlim(0.4, 9.6); b.set_ylim(-24, 1.5); b.set_xticks([1, 3, 5, 7, 9]); b.set_yticks([-20, -15, -10, -5, 0])
b.set_xlabel(r'band $a$  ($n=9$)', labelpad=1); b.set_ylabel(r'$\log_{10}(v_a/v_1)$', labelpad=1)
b.text(-0.36, 1.04, '(b)', transform=b.transAxes, fontsize=8); clean(b)
fig.savefig(os.path.join(HERE, 'figU3.pdf')); plt.close(fig)
print('wrote figU1.pdf figU2.pdf figU3.pdf;  velocity ratio per band: 10^%.3f' % slope[0])
