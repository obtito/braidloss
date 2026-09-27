"""Freeze a measured baseline only after its predeclared checks pass."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
from summarize_pilots import load_result
from recursive_geometry import RecursiveCable

ROOT=Path(__file__).resolve().parents[1]

def freeze():
    destination=ROOT/'baseline_frozen.json'
    if destination.exists():
        return json.loads(destination.read_text(encoding='utf-8'))
    names=['baseline_64_av_m2','baseline_64_av_m3','baseline_64_av_air10','baseline_64_av_length1']
    cases={}
    for name in names:
        folder=ROOT/'data/raw'/name
        run=json.loads((folder/'run.json').read_text(encoding='utf-8'))
        result=load_result(folder)
        if run['status']!='solved' or not result['complete']:
            raise RuntimeError('Unsolved verification: '+name)
        m=result['METRIC']
        power_error=float(100*max(abs(m[4].real/m[3].real-1)))
        if power_error>.1:raise RuntimeError('Power verification failed: '+name)
        cases[name]=dict(Rdc=float(m[4,0].real),Rac=float(m[4,-1].real),
                        K=float(m[4,-1].real/m[4,0].real),Pcu=float(m[6,-1].real),
                        power_error_pct=power_error,source_sha256=run['sha256'])
    base=cases[names[0]]
    changes={name:100*(cases[name]['Rac']/base['Rac']-1) for name in names[1:]}
    if any(abs(x)>=1 for x in changes.values()):
        raise RuntimeError('One-factor discretization/domain changes exceed 1%: '+str(changes))
    c=RecursiveCable((64,),(96.,),copper_area_mm2=6.,gap_um=40.)
    chosen=names[1]
    cfg=json.loads((ROOT/'data/raw'/chosen/'inputs.json').read_text(encoding='utf-8'))
    cfg.update(formal_baseline=True,validation_status='accepted_prespecified_baseline_checks')
    cfg['mapping']['electrical_mapping_status']='Full A-phi field continuity implemented and checked; see validation evidence'
    dc_error=100*(cases[chosen]['Rdc']/c.diagnostics()['Rdc_length_parallel_ohm_per_m']-1)
    if abs(dc_error)>.5:raise RuntimeError('DC length reference failed')
    mph=ROOT/'data/raw'/chosen/'model_fem.mph'
    manifest=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),case=chosen,
       meaning='Formal baseline frozen before design optimization; not an optimized design',
       fixed=dict(copper_area_normal_mm2=6,frequency_Hz=200000,current_RMS_A=20,
                  sigma_S_per_m=5.8e7,relative_permeability=1,temperature_C=20,
                  design_gap_um=40,magnetic_outer_boundary_radius_mm=6),
       geometry=c.diagnostics(),configuration=cfg,results=cases[chosen],
       validation=dict(one_factor_Rac_change_pct=changes,dc_reference_error_pct=dc_error,
                       accepted_limits=dict(Rac_change_pct=1,power_error_pct=.1,dc_reference_error_pct=.5),
                       verification_cases=cases),
       artifacts=dict(model=str(mph.relative_to(ROOT)).replace('\\','/'),
                      model_sha256=hashlib.sha256(mph.read_bytes()).hexdigest()))
    destination.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    return manifest

if __name__=='__main__':
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    print(json.dumps(freeze(),indent=2,ensure_ascii=False))
