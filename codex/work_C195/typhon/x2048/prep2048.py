"""Build the chi = 2048 continuation configs at L = 512 (C195x of codex/work_C195/RESULT.md)."""
import json, os
os.makedirs('cfg', exist_ok=True)
for d, tag in (('c195', 'p153c_L512_P0_rand1'), ('c195b', 'sym40_n316_L512_P1_rand1')):
    src = f'/scratch/lustre/jkudlerflam/c19x/{d}/out512'
    old = json.load(open(f'{src}/{tag}.json'))
    assert old['status'] == 'done' and old['cfg']['L'] == 512
    ck = f'{src}/ckpt/{tag}.npz'; assert os.path.exists(ck)
    last = [h for h in old['hist'] if 'Sbar6' in h][-1]; assert last['chi'] == 1536
    cfg = dict(old['cfg']); cfg.pop('conv', None)
    cfg['init'] = 'file:' + ck; cfg['sched'] = []
    cfg['final'] = dict(chi=2048, res_tol=1e-9, min_sweeps=3, max_sweeps=6, tolE=4e-9, tolS=5e-4)
    cfg['continues'] = dict(src=src, chi=1536, sweeps=len(old['hist']), E_phys=last.get('E_phys'), Sbar1=last['Sbar1'],
                            Sbar6=last['Sbar6'], Sbar12=last['Sbar12'], trunc=last['trunc'])
    json.dump(cfg, open(f'cfg/{tag}.json', 'w'), indent=1)
    print(tag, cfg['model'], cfg.get('gauge'), cfg['hz'], cfg['P'], cfg['final'])
