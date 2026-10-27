"""Create the final report from named, locally generated outputs."""
import json, math
from pathlib import Path
from d2 import entropy_fit, dump
import numpy as np

sources=['best_opt','candidate_refined','candidate854','candidate628','candidate871',
         'extra_candidates','extra_candidates_b','extra_candidates_c','structured_refined',
         'critical_refined','optimized16_refined','max_raw_refined','ferromagnet_half','structured_dm_refined']
allc={}
for name in sources:
    path=Path('results/'+name+'.json')
    if path.exists(): allc.update(json.loads(path.read_text()))
scans={name:json.loads(Path('results/'+name+'.json').read_text()) for name in
       ['global_scan','boundary_scan','optimization','optimization16','shortlist','global_shortlist']}
raw=[]
for name,data in scans.items():
    for p in data:
        for r in p.get('rows',[]): raw.append((r['c'],name,p,r))
maximum=max(raw,key=lambda x:x[0]); maxc,name,par,row=maximum

def metrics(item):
    rows=item['rows']; r=rows[-1]; tail=[x for x in rows if x['L']>=16]
    L=np.array([x['L'] for x in tail],float)
    gaps=np.array([np.pi*(x['v_q1_plus']+x['v_q1_minus'])/x['L'] for x in tail])
    if len(L)>=3 and min(gaps)>1e-10:
        coeff,cov=np.polyfit(np.log(L),np.log(gaps),1,cov=True)
        z=-float(coeff[0]); zse=float(np.sqrt(cov[0,0]))
    else: z=zse=float('nan')
    windows=[r['c'],entropy_fit(r['entropy'],r['L'],3)['c'],entropy_fit(r['entropy'],r['L'],max(2,r['L']//4))['c']]
    drift=max(abs(r['c']-x['c']) for x in rows[-3:])
    err=max(r['fit_se'],max(windows)-min(windows),drift)
    if len(L)>=3:
        halfcoef,halfcov=np.polyfit(np.log(L),[x['entropy'][-1] for x in tail],1,cov=True)
        half=float(3*halfcoef[0]); halferr=float(3*np.sqrt(halfcov[0,0]))
        for n in (3,4,5):
            if len(rows)<n: continue
            rr=rows[-n:]
            other=3*np.polyfit(np.log([x['L'] for x in rr]),[x['entropy'][-1] for x in rr],1)[0]
            halferr=max(halferr,abs(other-half))
    else: half=halferr=float('nan')
    return dict(c=r['c'],sensitivity=err,z_q1=z,z_fit_se=zse,c_half=half,c_half_sensitivity=halferr,L=r['L'],
                raw_c_energy=item.get('energy_fit',{}).get('c_energy'),
                c_MV_plus=r['stress'][0]['selected']['c_MV'],c_MV_minus=r['stress'][1]['selected']['c_MV'])

metrics_all={key:metrics(item) for key,item in allc.items()}
dump('results/candidate_metrics.json',metrics_all)
generic22=[key for key,item in allc.items() if item['rows'][-1]['L']==22 and not key.startswith('XXZ_')]
generic22.sort(key=lambda key:metrics_all[key]['c'],reverse=True)
dump('results/best_candidates.json',[dict(**allc[key]['parameters'],summary=metrics_all[key],rows=allc[key]['rows']) for key in generic22])
dump('results/best_candidate_points.json',[allc[key]['parameters'] for key in generic22])
largest22=max(generic22,key=lambda key:metrics_all[key]['c']) if generic22 else None
above=sum(p.get('c_eff',0)>1 for p in scans['global_scan'])
lines=['# D2 verdict and numerical results','',
       'No scale-consistent counterexample to c_T ≤ 1 has been established. This is a finite numerical search, not a proof of the conjecture. Large entropy slopes are abundant, but the examined outliers do not yet support a common conformal description with stable velocities and justified stress-tensor identification. Some remain unresolved beyond the available sizes.', '',
       f'The largest raw entropy slope found was c_eff={maxc:.6f} at L={row["L"]}, point {par["id"]}. Its within-chain regression standard error was {row["fit_se"]:.6f}; this is not a finite-size error bar. The largest L=12 slope was 2.268710 at point 854. Both points were selected for momentum-resolved refinement through L=22.', '',
       'The largest controlled central charge is still c=1 in the benchmark chains. The four reference models extrapolate to 1/2 or 1 with the conservative sensitivity estimates in calibration.md. No generic point is assigned a central charge above one here.', '',
       '![Calibration, candidate size dependence, and ferromagnetic crossover](results/diagnostics.png)', '',
       '## Scope', '',
       f'The global scan contains {len(scans["global_scan"])} parameter samples, 138 structured samples and 1,200 random canonical samples, at L=10,12. An additional {len(scans["boundary_scan"])} boundary samples use L=12,14. The pilot is a subset of the global scan and is not counted again. The first local optimization made {len(scans["optimization"])} objective evaluations; the second, at L=14,16 with a finite-momentum velocity penalty, made {len(scans["optimization16"])} evaluations.', '',
       f'A total of {len(scans["shortlist"])+len(scans["global_shortlist"])} global points were followed at L=14,16. The velocity screen covers 254 global and boundary samples. There are {len(generic22)} generic/structured parameter points with completed L=22 refinements, in addition to four calibration chains and the exact and near-ferromagnetic controls. All ground momentum sectors were searched for generic refinements; benchmark momentum and charge assumptions are stated in methods.md.', '',
       'The pilot follow-up selected c_eff(L=12)>1.12 or initial score>1.02. The subsequent global follow-up selected c_eff(L=12)>1.15 or score>1.05, excluding those already followed. These are screening choices, not physical exclusion criteria. The second optimizer used seeds 883 and 1144 and added 0.4 min(|z_q-1|,5) to the penalty. Every size sequence holds its couplings fixed; optimized couplings were frozen before the larger-L checks.', '',
       'The structured families are XXZ in a z field, XYZ in a general field, XYZ with DM exchange, uniform Kitaev-Heisenberg-Gamma exchange with a field, and Ising with transverse and longitudinal fields. A uniform symmetric exchange matrix is rotation-equivalent to XYZ, so uniform KHG is a chart of the same family. Alternating-bond KHG would not satisfy the one-site translation hypothesis and was not included.', '',
       'The boundary grids cover the XY/Ising/saturation meeting near gamma=0 and h_z=2, the XXZ boundaries near Delta=-1 and Delta=1, and DM versus pairing competition near D=gamma, including interacting variants. They locate finite-size crossovers and candidate critical regions, not thermodynamic phase boundaries with certified errors.', '',
       '## Largest-size diagnostics', '',
       'The uncertainty after c_eff is the largest of its fit standard error, cut-window spread, and drift across the last three sizes. The half-chain uncertainty includes regression error and changes among the last three, four, or five sizes. These are sensitivity scales, not confidence intervals. z_q is fitted from the average lowest gaps at relative momenta ±2 pi/L over L≥16; its displayed error is only the regression standard error. A conformal tower would require stable gap ratios and z_q tending to one. State changes and finite-size crossings can invalidate this simple exponent fit.', '',
       '| point | L | entropy c_eff, sensitivity | half-chain scaling slope | z_q | raw c_MV,+ | raw c_MV,- | naive energy c |',
       '|---|---:|---:|---:|---:|---:|---:|---:|']
for key in generic22:
    m=metrics_all[key]
    lines.append(f"| {key} | {m['L']} | {m['c']:.4f} ± {m['sensitivity']:.4f} | {m['c_half']:.3f} ± {m['c_half_sensitivity']:.3f} | {m['z_q1']:.3f} ± {m['z_fit_se']:.3f} | {m['c_MV_plus']:.4g} | {m['c_MV_minus']:.4g} | {m['raw_c_energy']:.3f} |")
if largest22:
    m=metrics_all[largest22]
    lines += ['',f'The largest raw L=22 entropy slope among these points is {m["c"]:.6f}, at {largest22}. This is an effective slope, not an identified c_T. The complete size sequences, cut-window comparisons, velocities, and parameters are in candidates.md and the JSON files.']
lines += ['', 'The negative or very large energy and matrix-element entries are diagnostic failures. They are not central charges and are not used as upper or lower bounds. A small gap in a non-stress state can make the proposed matrix-element ratio arbitrarily large. For multiple CFTs, a single selected stress branch would also miss part of c_T.', '',
          '## Parameters of the raw maximum', '', '```json',json.dumps({k:par[k] for k in ('id','J','field')},indent=2),'```','',
          '## Three checks that prevent an overclaim','',
          '1. The exact ferromagnetic boundary has a translation eigenstate with c_eff=1.360032 at L=22. Its energy is E0=-L, every nonzero Fourier energy-density mode annihilates it, and its magnon dispersion is 4(1-cos q), with z=2. The logarithmic entropy does not imply a CFT. At Delta=-0.999 the finite-size entropy slope remains near this value while the Fourier stress estimate is close to zero.', '',
          '2. The literal left-site-field, largest-overlap prescription selects a non-stress level in critical Ising at the tested sizes. At L=22 it returns 0.873678. Splitting the field equally between the two sites changes no Hamiltonian but selects the stress branch and returns 0.484874, extrapolating toward 1/2. The density improvement and its scaling are derived in methods.md. The main tables use the symmetric density, with the literal-density check retained separately.', '',
          '3. The naive Casimir fit fails even for a c=1 XX chain at finite density. Exact ring energies at L=8,...,22 give c_E=-0.38670 for h_z=0.3 and c_E=1.71226 for h_z=1.1, using the correct velocity. Particle-number changes add a zero-mode term of the same order as the Casimir energy. The largest single Fourier overlap also oscillates when finite-size curvature splits the conformal multiplet. Consequently, disagreement among the three raw formulas is not by itself proof that a generic candidate is nonconformal.', '',
          '## What was and was not completed','',
          'The ED toolkit, all-momentum resolution through L=22, calibration, structured and random scans, two rounds of local optimization, candidate refinements, and cluster array script were completed. Full-space and momentum-space spectra and Fourier actions agree to about 1e-14 on the independent small-system checks. On-site rotation, density-gauge invariance of H, charge-block entropy, and exact XX-ring energies were also checked. Numerical residuals are saved separately from finite-size uncertainty.', '',
          f'The literal request to take every value above one to L=22 with all three estimators was not completed: {above} of the {len(scans["global_scan"])} global samples already had c_eff(L=12)>1. The local budget was used for staged screening and the strongest or most scale-consistent candidates. The array script extends all global points to L=22. It was written but not submitted. No exhaustive statement about the continuous eight-parameter family is possible from these samples.', '',
          'The nearly degenerate boundary control b43 was completed through L=18. Its slow L=20 calculation was interrupted, with no L=20 result reported, to reserve time for the larger outliers. The boundary-file cluster array can finish it. Two earlier implementation/calibration failures and their superseded files are identified in README.md.', '',
          'The remaining scientific uncertainty is located at the extrapolation and conformal-family identification, not at the small-system matrix construction. Longer chains, zero-mode-aware spectral fits, and independent tensor-network checks are needed for the unresolved candidates. No theorem or numerical exclusion of c_T>1 is claimed.', '',
          'All files were created in the current directory. No pre-existing project files were read or listed. The local elapsed time and process completion record are in STATE.md.']
Path('report.md').write_text('\n'.join(lines)+'\n')
print('Wrote report.md; largest raw c_eff',maxc,'at',par['id'],'L',row['L'])
