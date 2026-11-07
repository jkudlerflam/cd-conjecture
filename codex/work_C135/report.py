"""Build complete coverage tables, explicitly retaining uncomputed entries."""
from core import *
def fmt(v):return 'NA' if v is None else f'{v:.7g}'
def report():
    records=[json.loads(p.read_text()) for p in sorted(Path('results').glob('*.json'))]
    dictionary=[]
    for rec in records:
        for r in rec['rings'].values():
            for p in r['petz']:
                if p['diagnostics']['trace_residual']<1e-8:dictionary.append(unpack(p['tau']))
    for s in P[1:]:
        for sign in [-1,1]:q=(I+sign*s)/2;dictionary.append(np.kron(q,q))
    lines=['# C135 numerical tables','', '[Numerical] All entries are floating-point computations, not interval certificates. NA means not computed. The trace norm has no factor 1/2. Energy bounds use corrected dual eigenvalues and repaired primal states; solver statuses and residuals are retained in results/*.json. A status other than optimal is exploratory even when a corrected bound has a small gap.','', '[Numerical] e0_up is the least energy over the shared dictionary of computed Petz-periodized marginals and six pure product states. The optional exact-free column uses the analytic dispersion integral. It is not a ring extrapolation. Bounds involving a numerical ring ground energy assume the eigensolver found the lowest eigenvalue; its vector residual is recorded but does not alone prove that assumption.','', '[Numerical] Q_low uses the strongest computed dual bound at n>L. K_low(n) uses lower(n+2)-upper(n), so subtraction has the correct direction. K_up uses e0_up-lower(n). A fixed finite dictionary has a variational offset, so growth of its K_up with n is not evidence of failure.','']
    summary=dict(models=len(records),random_models=sum(r['name'].startswith('random_') for r in records),dictionary_size=len(dictionary),energy_rows=0,lti_rows=0,distance_rows=0,inaccurate_lti=[],worst_Q_low=None,worst_Q_up=None,worst_K_low=None,worst_K_up=None,worst_distance_up=None,max_ring_residual=0,max_lti_residual=0,max_lti_dual_gap=0,max_petz_trace_residual=0,validation_failures=[],coverage=[])
    def update(key,value,where):
        if value is not None and (summary[key] is None or value>summary[key]['value']):summary[key]=dict(value=float(value),**where)
    tables={}
    for rec in records:
        name=rec['name'];h=unpack(rec['h']);M=rec['M']['upper'];Mlo=rec['M'].get('lower',M);up=min(float(np.trace(h@tau).real) for tau in dictionary)
        def normlow(v):return v/(M if v>=0 else Mlo)
        def normup(v):return v/(Mlo if v>=0 else M)
        exact=rec.get('exact_free_density');rowsA=[];rowsC=[]
        lines += [f'## {name}', '',f'[Numerical] M={M:.10g}, e0_up={up:.12g}. Computed ring lengths {", ".join(sorted(rec["rings"],key=int))}; LTI levels {", ".join(sorted(rec["lti"],key=int))}.','', '| L | Q lower | Q upper | Q exact free | n for lower |','| :--- | ---: | ---: | ---: | ---: |']
        for L in range(3,21):
            row=dict(L=L,lower=None,upper=None,exact_free=None,n=None)
            if str(L) in rec['rings']:
                r=rec['rings'][str(L)];e=r['e'];row['upper']=normup(L*L*(up-e));row['exact_free']=L*L*(exact-e)/M if exact is not None else None
                candidates=[(v['lower'],int(n)) for n,v in rec['lti'].items() if int(n)>L]
                if candidates:
                    lo,n=max(candidates);row['lower']=normlow(L*L*(lo-e));row['n']=n
                summary['energy_rows']+=1;summary['max_ring_residual']=max(summary['max_ring_residual'],r.get('residual') or 0)
                for p in r['petz']:summary['max_petz_trace_residual']=max(summary['max_petz_trace_residual'],p['diagnostics']['trace_residual'])
                for n,v in rec['lti'].items():
                    if int(n)<=L and v['lower']>e+1e-6:summary['validation_failures'].append([name,n,L,'LTI lower above compatible ring'])
                update('worst_Q_low',row['lower'],dict(name=name,L=L));update('worst_Q_up',row['upper'],dict(name=name,L=L))
                if row['lower'] is not None and row['lower']>row['upper']+1e-6:summary['validation_failures'].append([name,L,'reversed Q bracket'])
            rowsA.append(row);lines.append(f'| {L} | {fmt(row["lower"])} | {fmt(row["upper"])} | {fmt(row["exact_free"])} | {row["n"] or "NA"} |')
        lines+=['','| n | E_LTI numerical | K lower via n+2 | K upper via dictionary | K exact free | status |','| :--- | ---: | ---: | ---: | ---: | :--- |']
        for n in range(2,9):
            row=dict(n=n,value=None,lower=None,upper=None,exact_free=None,status='not computed')
            if str(n) in rec['lti']:
                v=rec['lti'][str(n)];row.update(value=v['value'],upper=normup(n*n*(up-v['lower'])),status=v['status']);summary['lti_rows']+=1
                if str(n+2) in rec['lti']:row['lower']=normlow(n*n*(rec['lti'][str(n+2)]['lower']-v['upper']))
                if exact is not None:row['exact_free']=n*n*(exact-v['value'])/M
                if v['status']!='optimal':summary['inaccurate_lti'].append([name,n,v['status']])
                summary['max_lti_residual']=max(summary['max_lti_residual'],v['lti_residual']);summary['max_lti_dual_gap']=max(summary['max_lti_dual_gap'],v['dual_gap'])
                update('worst_K_low',row['lower'],dict(name=name,n=n));update('worst_K_up',row['upper'],dict(name=name,n=n))
                if v['lower']>v['upper']+1e-6:summary['validation_failures'].append([name,n,'reversed LTI bracket'])
            rowsC.append(row);lines.append(f'| {n} | {fmt(row["value"])} | {fmt(row["lower"])} | {fmt(row["upper"])} | {fmt(row["exact_free"])} | {row["status"]} |')
        dp=Path(f'results_distance/{name}.json');dist=json.loads(dp.read_text()) if dp.exists() else dict(rows={})
        lines+=['','| L | L² d_low | L² d_up | lower status | hull status |','| :--- | ---: | ---: | :--- | :--- |']
        for L in range(3,9):
            r=dist['rows'].get(str(L),{});lo=r.get('L2_lower');hi=r.get('L2_upper')
            if hi is not None:
                summary['distance_rows']+=1;update('worst_distance_up',hi,dict(name=name,L=L))
                if lo is not None and lo>hi+1e-5:summary['validation_failures'].append([name,L,'reversed distance bracket'])
            lines.append(f'| {L} | {fmt(lo)} | {fmt(hi)} | {r.get("lti",{}).get("status","not computed")} | {r.get("hull",{}).get("status","not computed")} |')
        lines+=[''];tables[name]=dict(M=M,e0_up=up,directional=rowsA,lti=rowsC,distance=dist['rows'])
        summary['coverage'].append(dict(name=name,ring_max=max(map(int,rec['rings'])),lti_max=max(map(int,rec['lti'])),distance_lower_L=[int(L) for L,r in dist['rows'].items() if 'lti' in r]))
    Path('TABLES.md').write_text('\n'.join(lines)+'\n');write_json('tables.json',tables);write_json('summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['coverage','inaccurate_lti']},indent=2));print('Inaccurate LTI entries:',len(summary['inaccurate_lti']))
if __name__=='__main__':report()
