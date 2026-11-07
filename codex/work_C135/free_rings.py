"""Analytic finite free-chain spectra, no exponential Hilbert space."""
from core import *
from scipy.special import ellipe
def xx_energy(L,b):
    vals=[]
    for m in range(L+1):
        phi=0 if m%2 else .5
        eps=np.sort(np.cos(2*np.pi*(np.arange(L)+phi)/L))
        vals.append(float(eps[:m].sum()+b*(m-L/2)))
    return min(vals)
def tfi_energy(L,g):
    k=2*np.pi*(np.arange(L)+.5)/L
    return float(-np.sqrt(1+g*g-2*g*np.cos(k)).sum())
checks=[]
for path in Path('results').glob('*.json'):
    rec=json.loads(path.read_text());name=rec['name']
    if name in ['XX','XX_field']:
        b=0 if name=='XX' else .5
        fn=lambda L:xx_energy(L,b)
        exact=-(np.sqrt(1-b*b)+b*np.arcsin(b))/np.pi
    elif name.startswith('TFIM'):
        h=unpack(rec['h']);g=float(-np.trace(h@B[8]).real/(2*np.trace(B[8]@B[8]).real))
        fn=lambda L:tfi_energy(L,g)
        exact=-2*(1+g)*ellipe(4*g/(1+g)**2)/np.pi
    else:continue
    discrepancy=max(abs(fn(int(L))-r['E']) for L,r in rec['rings'].items())
    assert discrepancy<1e-8
    for L in range(3,21):
        if str(L) not in rec['rings']:
            e=fn(L);rec['rings'][str(L)]=dict(L=L,E=e,e=e/L,residual=None,method='analytic free spectrum, evaluated in floating point',petz=[])
    rec['exact_free_density']=float(exact);rec['free_spectrum_ED_discrepancy']=discrepancy
    write_json(path,rec);checks.append(dict(name=name,maximum_ED_difference=discrepancy,e0=float(exact)))
write_json('free_ring_checks.json',checks);print('Free ring formulas match all available ED results; extended energies through L=20.')
