import json, sys
for l in open(sys.argv[1]):
    d = json.loads(l)
    if 'init' in d or 'start' in d or 'resume' in d: continue
    print(d['chi'], d['t'], '%.11f' % d['e_site'], 'S=%.5f' % d['Sbar'], 'N=%.4f var=%.4f' % (d['Ncell'], d['varN']),
          'rp=%.1e fid=%.7f res=%.1e mv=%.0f D=%d disc=%.1e %.0fs' % (d['r_pred'], d['fid0'], d['res_max'], d['mv_mean'], d['D'], d['disc'], d['secs']))
