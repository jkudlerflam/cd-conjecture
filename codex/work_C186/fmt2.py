import sys, json
for l in sys.stdin:
    try: d = json.loads(l)
    except Exception: print(l, end=''); continue
    print(f"u={d['u']} v={d['v']:.4f} r={d['r']:.3f} {d['rmode']:6s} L={d['L']} Nmax={d['Nmax']} phi={d['phi']:.4f} mu={d['mu']:.4f} qa={d['qa']:.3f} qb={d['qb']:.3f} ea={d['ea']:.4f} eb={d['eb']:.4f} K={d['K']:.3f} T={d['T']:.4f} Elow={d['E_low']:.5f} N={d['N_low']:.2f} B={d['B']:+.5f} lv={[round(x,4) for x in d['levels']]} Ns={[round(x,2) for x in d['Ns']]} {d['sec']}s", flush=True)
