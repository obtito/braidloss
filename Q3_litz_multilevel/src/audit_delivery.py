"""Check the delivered values against raw logs and the preserved evidence.

This checks data integrity, not a new proof of electromagnetic accuracy.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import csv
import hashlib
import json
import re

import numpy as np
from finalize_study import file_hash
from summarize_pilots import load_result
from optimize_design import designs,design_id
from validate_rom import source_digest

ROOT=Path(__file__).resolve().parents[1]


def audit():
    d=json.loads((ROOT/'final_results.json').read_text(encoding='utf-8'))
    checks=[]
    def require(condition,description):
        if not condition:raise AssertionError(description)
        checks.append(description)
    cases={d['baseline']['case']:d['baseline'],d['recommended']['case']:d['recommended']}
    cases.update(d['verification']['cases'])
    cases.update({x['COMSOL']['case']:x['COMSOL'] for x in d['confirmation_cases']})
    for name,expected in cases.items():
        folder=ROOT/'data/raw'/name
        run=json.loads((folder/'run.json').read_text(encoding='utf-8'))
        cfg=json.loads((folder/'inputs.json').read_text(encoding='utf-8'))
        raw=load_result(folder);m=raw['METRIC']
        require(run['status']=='solved' and raw['complete'],name+': solver completion recorded')
        require(file_hash(folder/run['java'])==run['sha256'],name+': original Java matches run hash')
        for file,digest in expected['evidence_sha256'].items():
            require(file_hash(folder/file)==digest,name+': '+file+' matches final evidence hash')
        require(cfg['total_strands']==64 and cfg['copper_area_normal_mm2']==6
                and cfg['gap_design_um']==40,name+': geometric control variables unchanged')
        require(cfg['electrical']=='full_A_phi_periodic' and cfg['mapping']['geometry_mapping_pass'],
                name+': expected periodic field formulation and geometry map')
        actual={'Rdc':m[4,0].real,'Rac':m[4,-1].real,
                'K':m[4,-1].real/m[4,0].real,'Pcu':m[6,-1].real}
        require(all(np.isclose(expected[k],x,rtol=1e-11,atol=1e-13) for k,x in actual.items()),
                name+': final numerical values reproduced from raw METRIC rows')
        require(np.allclose(m[0].real,[10,200000]),name+': 10 Hz / 200 kHz frequencies')
        require(abs(raw['STRAND'][::3,-1].sum()-20)<1e-7,name+': strand phasors sum to 20 A RMS')
        require(np.isclose(actual['Pcu'],400*actual['Rac'],rtol=1e-12),name+': P = I_RMS^2 Rac')
    records=list(csv.DictReader((ROOT/'parameter_search.csv').open(encoding='utf-8-sig')))
    require(len(records)==129 and {x['id'] for x in records}=={design_id(f,p) for f,p in designs()},
            'All 129 predeclared candidates appear exactly once')
    require(dict(Counter(x['status'] for x in records))==d['search_counts'],
            'Candidate counts match final report')
    predictions=list((ROOT/'data/search_64').glob('*.json'))
    require(all(json.loads(p.read_text())['rom_sha256']==source_digest() for p in predictions),
            'All retained predictions correspond to the recorded ROM implementation')
    for key in ['Rac','K','Pcu']:
        value=100*(1-d['recommended'][key]/d['baseline'][key])
        require(np.isclose(value,d['improvement_pct'][key],rtol=1e-12),key+': improvement recomputed')
    require(all(abs(x)<1 for x in d['verification']['Rac_change_pct'].values()),
            'Declared successive integral-change checks pass (not an error bound)')
    require(d['verification']['mesh_history'][0]['status']=='exceeded_1_pct_further_refinement_required',
            'Original failed mesh comparison preserved')
    require(d['post_search_ROM_assessment']['status']=='restricted_after_new_FEM_evidence',
            'Post-search ROM limitation is explicit')
    field_summary=[]
    for name in [d['baseline']['case'],d['recommended']['case']]:
        folder=ROOT/'data/raw'/name
        sample=np.load(folder/'sampling_points.npz')
        lines=(folder/'stdout.log').read_text(encoding='utf-8',errors='replace').splitlines()
        fields=np.array([[float(v) for v in line.split('|')[1:]] for line in lines if line.startswith('FIELD|')])
        require(len(fields)==len(sample['points_m']) and np.isfinite(fields).all(),name+': field sample count and finiteness')
        require(np.array_equal(fields[:,0],np.arange(len(fields))),name+': field sample indices ordered')
        require(np.allclose(fields[:,1:4],sample['points_m'],rtol=1e-12,atol=1e-14),name+': field sample coordinates match source')
        require(fields[:,-1].min()>=0,name+': RMS field magnitudes are nonnegative')
        # Field coordinates are actual common-z cuts, not projected normal disks.
        cross=sample['cross_m']
        require(np.ptp(cross[...,2])<1e-10,name+': common-z cut geometry verified')
        field_summary.append(dict(case=name,count=len(fields),J_RMS_min_A_per_mm2=fields[:,-1].min()/1e6,
                                  J_RMS_max_A_per_mm2=fields[:,-1].max()/1e6))
    for filename in ['README.md','进度与控制变量.md','问题三_结果与控制变量.md']:
        text=(ROOT/filename).read_text(encoding='utf-8')
        links=re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',text)
        missing=[p for p in links if not p.startswith(('https://','http://','#')) and not (ROOT/p).exists()]
        require(not missing,filename+': all local links resolve')
    sources={str(p.relative_to(ROOT)).replace('\\','/'):file_hash(p) for p in (ROOT/'src').rglob('*')
             if p.is_file() and p.suffix in ['.py','.txt']}
    outputs={p.name:file_hash(p) for p in ROOT.iterdir() if p.is_file() and p.suffix in ['.md','.pdf','.csv']}
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),status='data_integrity_checks_passed',
        passed_count=len(checks),checks=checks,field_samples=field_summary,source_sha256=sources,
        output_sha256=outputs,
        limitations='Numerical acceptance is conditional on the reported discretization checks; integrity checks do not prove a unique optimum or validate all 129 predictions.')
    (ROOT/'delivery_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    return result


if __name__=='__main__':
    result=audit()
    print(json.dumps({'status':result['status'],'passed_count':result['passed_count']},indent=2))
