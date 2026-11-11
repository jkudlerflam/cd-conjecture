"""Build the chi = 1536 continuation configs at L = 256 (REPORT/RESULT R13). Reads outX1024 records, writes cfg/<tag>.json."""
import json, os
SRC = '/scratch/lustre/jkudlerflam/c18x/c190g/outX1024'
os.makedirs('cfg', exist_ok=True)
for tag in ('p153c_L256_P1_rand1', 'sym40_n316_L256_P1_rand1'):
    old = json.load(open(f'{SRC}/{tag}.json'))
    assert old['status'] == 'done' and old['hist'][-1]['chi'] == 1024 and old['cfg']['L'] == 256
    ck = f'{SRC}/ckpt/{tag}.npz'; assert os.path.exists(ck)
    cfg = dict(old['cfg']); cfg.pop('conv', None)
    cfg['init'] = 'file:' + ck
    cfg['sched'] = []
    cfg['final'] = dict(chi=1536, res_tol=1e-9, min_sweeps=3, max_sweeps=8, tolE=2e-9, tolS=5e-4)
    last = old['hist'][-1]
    cfg['continues'] = dict(src=SRC, chi=1024, sweeps=len(old['hist']), Sbar6=last['Sbar6'], E=last['E'], trunc=last['trunc'])
    json.dump(cfg, open(f'cfg/{tag}.json', 'w'), indent=1)
    print(tag, {k: cfg[k] for k in ('model', 'hz', 'L', 'P', 'gauge') if k in cfg}, cfg['final'])
