import sys, json
for l in sys.stdin:
    try: d = json.loads(l)
    except Exception: print(l, end=''); continue
    if 'tuned' in d:
        t = d['tuned']; print(f"TUNED u={t['u']} v={t['v']} {t['rmode']} phi={t['phi']:.5f} mu={t['mu']:.5f} ka={t['ka']:.4f} kb={t['kb']:.4f} ea={t['ea']:.1e} eb={t['eb']:.1e}", flush=True); continue
    print(f"u={d['u']} v={d['v']} {d['rmode']} L={d['L']} Nmax={d['Nmax']} K={d['K']:.3f} T={d['T']:.5f} Elow={d['E_low']:.5f} N={d['N_low']:.2f} B={d['B']:.5f} E1min={d['E1min']:.4f}@{d['q1min']:.3f} levels={[round(x,4) for x in d['levels']]} Ns={[round(x,2) for x in d['Ns']]} {d['sec']}s", flush=True)
