"""C181: summary of run.py outputs. Per tag (last sweep): L, N, E, truncation, S_mid, c_prof, valley polarization p,
density homogeneity, dominant momentum of the pair correlation D(r) (separation 2, c frame), |G_f| decay.
Pairs of sizes with the same label give c_two. Parity gaps from tags with N-1, N, N+1.
Usage: python3 summarize.py <dir> [<dir> ...]"""
import json, glob, sys, os, re, numpy as np
import analyze as A
recs = {}
for d in sys.argv[1:]:
    for f in sorted(glob.glob(os.path.join(d, '*.jsonl'))):
        lines = open(f).read().strip().split('\n')
        if not lines or not lines[-1]: continue
        r = json.loads(lines[-1]); recs[r['tag']] = r
def lab(tag):
    m = re.match(r'(.*)_L(\d+)_N(\d+)$', tag); return m.group(1), int(m.group(2)), int(m.group(3))
rows = {}
print("# tag sweeps chi E trunc Smid c_prof p homog(maxdev) Qpair Kp_dist |Gf|(L/8)/|Gf|(L/4)")
for tag, r in sorted(recs.items()):
    L, N = r['cfg']['L'], r['cfg']['N']
    S = np.array(r['S'])
    out = dict(E=r['E'], tr=r['trunc'], S=S, sweep=r['sweep'], chi=r['chi'], L=L, N=N)
    s = f"{tag:32s} {r['sweep']:2d} {r['chi']:4d} {r['E']:.8f} {r['trunc']:.1e} {S[L//2-1]:.4f} {A.c_prof(S, L):.3f}"
    if 'meas' in r:
        m = r['meas']
        Wp, Wm = m['Wa_f'], m['Wb_f']; p = abs(Wp - Wm) / (Wp + Wm)
        pb = abs(m['Wa_b'] - m['Wb_b']) / (m['Wa_b'] + m['Wb_b'])
        d = np.array(m['dens']); k = np.ones(16) / 16; ds = np.convolve(d, k, 'valid')
        c = ds[len(ds)//4: 3*len(ds)//4]; hom = (c.max() - c.min()) / c.mean()
        # pair correlation from the first reference, separation 2
        key = sorted([kk for kk in m['Dpair'] if kk.startswith('2_')])[0]
        D = np.array(m['Dpair'][key]); Dc = D[:, 0] + 1j * D[:, 1]
        n = len(Dc); qs = np.linspace(-np.pi, np.pi, 721)
        rr = np.arange(n); win = np.hanning(n)
        F = np.abs([(Dc * win * np.exp(-1j * q * rr)).sum() for q in qs])
        Q = qs[np.argmax(F)]
        G = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
        W = G.shape[0]; g = np.abs(np.array([np.mean(np.abs(np.diag(G, k))) for k in range(W)]))
        rat = g[max(1, L // 8)] / g[max(1, L // 4)] if L // 4 < W else np.nan
        out.update(p=p, pb=pb, hom=hom, Q=Q, rat=rat)
        s += f" p={p:.3f} pb={pb:.3f} hom={hom:.3f} Q={Q:+.3f} Gratio={rat:.2f}"
    print(s)
    rows[tag] = out
print("# two-size and parity")
groups = {}
for tag, o in rows.items():
    b, L, N = lab(tag); groups.setdefault(b, {})[(L, N)] = o
for b, g in sorted(groups.items()):
    Ls = sorted(set(L for (L, N) in g))
    main = {}
    for L in Ls:
        Ns = sorted(N for (LL, N) in g if LL == L)
        if len(Ns) == 3: main[L] = Ns[1]
        elif len(Ns) == 1: main[L] = Ns[0]
    for L in Ls:
        if L in main and all((L, main[L] + d) in g for d in (-1, 1)):
            N = main[L]; dP = g[(L, N+1)]['E'] + g[(L, N-1)]['E'] - 2 * g[(L, N)]['E']
            print(f"{b:28s} L={L} N={N} parity gap {dP:.5f} L*gap {L*dP:.3f}")
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            L1, L2 = Ls[i], Ls[j]
            if L1 in main and L2 in main and abs(main[L1] / L1 - main[L2] / L2) < 1e-9:
                S1, S2 = g[(L1, main[L1])]['S'], g[(L2, main[L2])]['S']
                cs = [6 * (A.Sbar(S2, L2, W) - A.Sbar(S1, L1, W)) / np.log(L2 / L1) for W in (1, 6, 12)]
                print(f"{b:28s} c_two({L1}->{L2}) W=1,6,12: " + " ".join(f"{c:.3f}" for c in cs))
