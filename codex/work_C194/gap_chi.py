"""C194: Delta_pb at the end of each chi stage (convergence of the gap in chi), for every complete triple.
Usage: python3 -E -s gap_chi.py [outdir ...]"""
import json, glob, os, sys
st = {}
for d in sys.argv[1:] or ['out_main']:
    for f in glob.glob(os.path.join(d, '*.jsonl')):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        fin = [r for r in recs if r.get('kind') == 'final']
        if not fin: continue
        c = fin[0]['cfg']; sw = [r for r in recs if r.get('kind') == 'sweep']
        stage = {}
        for r in sw: stage[r['chi']] = r['E']
        st[(c['nom'], c['L'], c['N'], d)] = stage
for (nom, L, N, d), s in sorted(st.items(), key=lambda kv: (eval(kv[0][0]), kv[0][1])):
    if N % 2: continue
    m, p = st.get((nom, L, N - 1, d)), st.get((nom, L, N + 1, d))
    if not (m and p): continue
    chis = sorted(set(s) & set(m) & set(p))
    print(f"n={nom:5s} L={L:4d} " + '  '.join(f"chi{c}: {L*(p[c]+m[c]-2*s[c]):.4f}" for c in chis if c >= 64) + f"   [{d}]")
