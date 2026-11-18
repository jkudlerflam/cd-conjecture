from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIOR=ROOT/'contrib/v4_bulk_diagnosis_20261117'
OLD=ROOT/'contrib/simple_qubit_candidates_20261117'

def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sources():
    files=list(PRIOR.rglob('*.json'))+[OLD/'outputs'/f'legacy_V4_a1_L96_P{P}.json' for P in (0,1)]
    files=[p for p in files if 'runtime' not in p.parts]
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(files)}

def candidate_records():
    result=[]
    for r in read(PRIOR/'analysis.json')['best_original_states']:
        d=read(PRIOR/r['file']);S=d.get('S')
        if S is None:S=read(OLD/'outputs'/f"legacy_V4_a1_L96_P{d['P']}.json")['S']
        result.append(dict(tag=Path(r['file']).stem,L=d['L'],P=d['P'],chi=d['chi'],E=d['E'],
            S=S,operators=d['operators'],state_status='completed',source=r['file']))
    raw=read(PRIOR/'outputs/v4_L144_P1_s1_c320.json');last=raw['hist'][-1]
    # No completed chi320 bond matrix exists; do not pair it with a chi256
    # profile and represent that as a same-state residual test.
    result.append(dict(tag='v4_L144_P1_s1_c320_interrupted',L=144,P=1,chi=320,E=last['E_phys'],
        S=last['Sprof'],density=[(1-z)/2 for z in last['Zprof']],state_status='interrupted',
        source='outputs/v4_L144_P1_s1_c320.json'))
    return result
