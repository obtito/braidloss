"""Predeclared finite parameter search, gated by baseline and independent FEM.

It minimizes the physical impedance result, with explicit geometry bounds.
There is no fitted pitch penalty, equal-current input, or unconstrained change
of copper area. A result here is a screening candidate, not a verified optimum.
"""
from pathlib import Path
from math import lcm
from itertools import product
import csv
import json
import time
import numpy as np
from recursive_geometry import RecursiveCable
from clearance_geometry import periodic_curve_clearance
from multipole_rom import evaluate
from validate_rom import source_digest

ROOT=Path(__file__).resolve().parents[1]

def designs():
    for outer in [64.,96.,128.]:yield (64,),(outer,)
    for inner,outer,sign in product([16.,32.,64.],[64.,96.,128.],[-1,1]):
        yield (4,16),(inner*sign,outer)
    for inner,middle,outer,s1,s2 in product([16.,32.,64.],[16.,32.,64.],[64.,96.,128.],[-1,1],[-1,1]):
        yield (4,4,4),(inner*s1,middle*s2,outer)

def design_id(f,p):
    return 'g'+'x'.join(str(x) for x in f)+'_p'+'_'.join(('m' if x<0 else '')+str(int(abs(x))) for x in p)

def cell_length(f,p):
    if len(f)==1:return .0005
    # All inner group factors in this predeclared study have fourfold symmetry.
    # This formula only proposes a cell; geometry and tangents are verified next.
    lengths=[round(abs(x)*1000/4) for x in p[:-1]]
    return lcm(*lengths)*1e-6

def run_search():
    validation=json.loads((ROOT/'rom_validation.json').read_text(encoding='utf-8'))
    if validation['status']!='screening_accepted' or validation['source_sha256']!=source_digest():
        raise RuntimeError('Independent 3D validation has not accepted this ROM version')
    baseline=json.loads((ROOT/'baseline_frozen.json').read_text(encoding='utf-8'))
    post_file=ROOT/'data/rom_post_confirmation.json'
    post=json.loads(post_file.read_text(encoding='utf-8')) if post_file.exists() else {}
    out=ROOT/'data/search_64';out.mkdir(exist_ok=True)
    rows=[]
    for f,p in designs():
        name=design_id(f,p);path=out/(name+'.json')
        if path.exists():
            record=json.loads(path.read_text(encoding='utf-8'))
            if record['rom_sha256']!=source_digest():raise RuntimeError('Old search record has a different ROM version')
            rows.append(record);continue
        if post.get('status')=='restricted_after_new_FEM_evidence':
            raise RuntimeError('New FEM evidence exceeds the 5% ROM target. Existing records may be replayed, but new predictions require further validation.')
        started=time.perf_counter()
        c=RecursiveCable(f,p,gap_um=40.)
        d=c.diagnostics();ell=cell_length(f,p)
        record=dict(id=name,factors=list(f),signed_pitches_mm=list(p),cell_length_mm=ell*1000,
                    rom_sha256=source_digest(),geometry=d)
        exclude=[]
        if d['outer_diameter_bound_mm']>6:exclude.append('outer_diameter_above_6_mm')
        if d['growth_mean']>1.10:exclude.append('mean_length_growth_above_10_pct')
        if d['max_tangent_angle_deg']>30:exclude.append('tangent_angle_above_30_deg')
        if d['min_bend_radius_mm']<10*2*c.a*1000:exclude.append('bend_radius_below_10_wire_diameters')
        if not exclude:
            clearance=periodic_curve_clearance(c)
            record['clearance']=clearance
            if clearance['clearance_lower_bound_um']<10:exclude.append('spline_clearance_below_10_um')
        mapping=c.screw_mapping(ell,2*np.pi*ell/(p[-1]*1e-3))
        record['mapping']=mapping
        if not mapping['geometry_mapping_pass']:exclude.append('periodic_geometry_mapping_failed')
        if exclude:
            record.update(status='excluded_geometry',reasons=exclude)
        else:
            r=evaluate(c,ell,slices=16,modes=3)
            record['ROM']={k:r[k] for k in ['Rac','Rdc','K','imbalance_complex','passivity_min_eigenvalue','slices','modes']}
            record['ROM'].update(Pcu_20A_W_per_m=400*r['Rac'],Z_real=r['Z'].real,Z_imag=r['Z'].imag,
                currents_fraction_real=r['currents_fraction'].real.tolist(),
                currents_fraction_imag=r['currents_fraction'].imag.tolist())
            if r['passivity_min_eigenvalue']<=0:record.update(status='rejected_nonpassive',reasons=['real_impedance_not_positive'])
            elif r['Rac']>baseline['results']['Rac']:record.update(status='excluded_loss',reasons=['Rac_above_frozen_baseline'])
            else:record.update(status='feasible_screening_candidate',reasons=[])
        record['elapsed_seconds']=time.perf_counter()-started
        path.write_text(json.dumps(record,indent=2),encoding='utf-8')
        rows.append(record)
        print(name,record['status'],record.get('ROM',{}).get('K',''),flush=True)
    acceptable=sorted([r for r in rows if r['status']=='feasible_screening_candidate'],key=lambda r:r['ROM']['K'])
    if not acceptable:raise RuntimeError('No feasible candidate')
    fields=['id','structure','pitches_mm','status','reason','Rdc','Rac','K','Pcu','OD_mm','length_growth_pct','angle_deg','clearance_um','time_s']
    with (ROOT/'parameter_search.csv').open('w',newline='',encoding='utf-8-sig') as fp:
        writer=csv.DictWriter(fp,fieldnames=fields);writer.writeheader()
        for r in rows:
            g=r['geometry'];v=r.get('ROM',{})
            writer.writerow(dict(id=r['id'],structure='x'.join(map(str,r['factors'])),pitches_mm='/'.join(map(str,r['signed_pitches_mm'])),
                status=r['status'],reason=';'.join(r['reasons']),Rdc=v.get('Rdc'),Rac=v.get('Rac'),K=v.get('K'),Pcu=v.get('Pcu_20A_W_per_m'),
                OD_mm=g['outer_diameter_bound_mm'],length_growth_pct=100*(g['growth_mean']-1),angle_deg=g['max_tangent_angle_deg'],
                clearance_um=r.get('clearance',{}).get('clearance_lower_bound_um'),time_s=r['elapsed_seconds']))
    result=dict(status='ROM_screening_only_COMSOL_confirmation_required',evaluated=len(rows),
                feasible=len(acceptable),best_candidate=acceptable[0],top_ten=acceptable[:10],
                range_description='129 predeclared 64-strand designs in DESIGN_SPACE.md; not a global optimum over strand count or arbitrary topology')
    (ROOT/'search_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result

if __name__=='__main__':
    result=run_search()
    print('SCREENING_COMPLETE',result['best_candidate']['id'],result['best_candidate']['ROM']['K'])
