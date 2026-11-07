import json
from pathlib import Path
import numpy as np
rows=[json.loads(s) for s in Path('ed_results.jsonl').read_text().splitlines()]
groups={}
for a in rows:groups.setdefault((a['family'],a['parameter']),[]).append(a)
lines=['# C133 full-ring numerical results','', '[numerical] Ground energies come from complete d^L Hilbert spaces, without charge-sector restriction. Reflection symmetry proves the M used here is gauge optimal. Each Q uses the same bond on L and 2L. A_eff=4Q/pi is a finite-size estimator, not an established infrared A. The independent fit uses E_L/L=e_inf-pi A/(6 L^2)+b/L^4 where at least four compatible sizes are available. Three sizes give no fit diagnostic, so no extrapolation is reported for them.','', '| Family | Parameter | d | Doubling L | Q_L | A_eff/M | Fit A/M | pi H_(d-1) |','| :-- | --: | --: | :-- | :-- | :-- | :-- | --: |']
summary=[]
for (fam,p),rs in groups.items():
    rs.sort(key=lambda a:a['L']);m=rs[0]['M'];d=rs[0]['d'];by={a['L']:a for a in rs}
    pairs=[(l,l*(by[2*l]['E']-2*a['E'])/m) for l,a in by.items() if 2*l in by]
    # SU(4) vacuum requires multiple four, ULS multiple three, TB/XXZ even.
    if fam=='spin_orbital': pairs=[(l,q) for l,q in pairs if l%4==0]
    if fam=='spin1_blbq': pairs=[(l,q) for l,q in pairs if l%(2 if p<0 else 3)==0]
    Afit=None
    if len(rs)>=4:
        x=1/np.array([a['L'] for a in rs],float)**2;y=np.array([a['E']/a['L'] for a in rs]);mat=np.column_stack([np.ones(len(x)),x,x*x]);fit=np.linalg.lstsq(mat,y,rcond=None)[0];Afit=float(-6*fit[1]/(np.pi*m))
    env=float(np.pi*sum(1/j for j in range(1,d)))
    lines.append(f'| {fam} | {p:.6g} | {d} | '+', '.join(str(l) for l,q in pairs)+' | '+', '.join(f'{q:.6f}' for l,q in pairs)+' | '+', '.join(f'{4*q/np.pi:.6f}' for l,q in pairs)+(f' | {Afit:.6f}' if Afit is not None else ' | unresolved')+f' | {env:.6f} |')
    summary.append(dict(family=fam,parameter=p,d=d,M=m,pairs=[dict(L=l,Q=q,Aeff_over_M=4*q/np.pi) for l,q in pairs],fit_A_over_M=Afit,envelope=env))
lines+=['','[numerical] Largest eigenvector residual: '+f'{max(a["residual"] for a in rows):.3e}. '+f'{len(rows)} ring ground states. Maximum Hilbert dimension {max(a["dimension"] for a in rows)}.','', '[assumed, infrared interpretation] The self-dual ferromagnetic FZ and Ashkin-Teller rings use the untwisted vacuum, XXZ and the spin-1 TB point use even rings, ULS uses multiples of three, and spin-orbital comparisons use multiples of four. These are the known vacuum-compatible sequences at the named critical points. Away from known critical points Q is still a finite-ring observable, and is not labeled an infrared central-charge measurement.','', '[numerical, limitation] The fit at XXZ Delta=1.05 remains close to critical on these rings although its asymptotic phase is gapped. The Heisenberg finite-size estimator exceeds pi slightly. These are controls against identifying small-ring drift with a counterexample. The q=5,6 results only reach L=6 and do not establish asymptotic A. No DMRG was run.']
Path('NUMERICS.md').write_text('\n'.join(lines)+'\n');Path('summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('\n'.join(lines))
