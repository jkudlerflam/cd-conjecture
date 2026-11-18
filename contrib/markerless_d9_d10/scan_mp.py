import mpmath as mp, sys, json
from analyze_mp import analyse
mp.mp.dps = int(sys.argv[2]) if len(sys.argv)>2 else 60
n = int(sys.argv[1])
for a in sys.argv[3].split(','):
    for D in sys.argv[4].split(','):
        for eta in sys.argv[5].split(','):
            for fg in sys.argv[6].split(','):
                r = analyse(n, a, D, mp.mpf(eta), mp.mpf(fg))
                e = r['eigS_over_T']
                print(f"n={n} a={a} D={D} eta={eta} f={fg}: eigS min/max={e[0]/e[-1]:.2e} antiherm={r['antiherm']:.2e} mu/T={r['mu_over_T']:.2e} "
                      f"on_circle={r.get('on_circle')}/{2*n} d2={r.get('d2',0):.3e} d4={r.get('d4',0):.3e} log10T={r['log10T']:.1f}", flush=True)
