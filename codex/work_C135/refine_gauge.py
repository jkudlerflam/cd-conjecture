"""Direct dual scale refinement, independent of primal cone multipliers."""
from core import *
for path in Path('results').glob('*.json'):
    rec=json.loads(path.read_text());h=unpack(rec['h']);m=rec['M']
    if m['gap']<1e-7:continue
    d=cp.Variable((4,4),hermitian=True)
    cons=[cp.real(cp.trace(d))==0,cp.normNuc(d)<=1]+[cp.real(cp.trace(d@g))==0 for g in G]
    p=cp.Problem(cp.Maximize(cp.real(cp.trace(h@d))),cons);p.solve(solver='SCS',eps=2e-10,max_iters=30000)
    D=herm(d.value);D-=np.trace(D)*np.eye(4)/4
    for g in G:D-=np.trace(D@g).real/np.trace(g@g).real*g
    D/=max(1.,abs(la.eigvalsh(D)).sum());lo=float(np.trace(D@h).real)
    m['direct_dual_status']=p.status;m['lower']=max(m['lower'],lo);m['gap']=m['upper']-m['lower'];write_json(path,rec)
    print(rec['name'],m['gap'],p.status)
