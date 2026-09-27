"""Record blind ROM predictions and reconcile independent 3D results."""
from pathlib import Path
import hashlib
import json
import itertools
import numpy as np
from recursive_geometry import RecursiveCable
from multipole_rom import evaluate
from run_rom_validation import CASES
from summarize_pilots import load_result

ROOT=Path(__file__).resolve().parents[1]
FEM_ALIASES={'validate_two_64_16_96':'validate_two_64_16_96_curvature_m3'}

def source_digest():
    return hashlib.sha256((Path(__file__).with_name('multipole_rom.py')).read_bytes()).hexdigest()

def check():
    baseline=json.loads((ROOT/'baseline_frozen.json').read_text(encoding='utf-8'))
    parameters=[(baseline['case'],dict(factors=(64,),pitches=(96.,),length_mm=.5))]+CASES
    prediction_file=ROOT/'data/rom_blind_predictions_v2.json'
    if prediction_file.exists():
        stored=json.loads(prediction_file.read_text(encoding='utf-8'))
        if stored['source_sha256']!=source_digest():
            raise RuntimeError('ROM changed after blind predictions; retain old evidence and create a new validation version')
        predictions=stored['predictions']
    else:
        predictions={}
        for name,params in parameters:
            c=RecursiveCable(params['factors'],params['pitches'],gap_um=40.)
            r=evaluate(c,params['length_mm']*1e-3,slices=16,modes=3)
            predictions[name]={k:r[k] for k in ['Rac','Rdc','K','imbalance_complex','passivity_min_eigenvalue','slices','modes']}
        stored=dict(source_sha256=source_digest(),fit_coefficients=None,
                    note='No empirical fitting. Baseline and positive-pitch two-level case are development checks; two three-level cases and the new negative-pitch two-level case are holdouts.',
                    predictions=predictions)
        prediction_file.write_text(json.dumps(stored,indent=2),encoding='utf-8')
    rows=[]
    for name,params in parameters:
        folder=ROOT/'data/raw'/FEM_ALIASES.get(name,name)
        role='development' if name in [baseline['case'],'validate_two_64_16_96'] else 'independent'
        row=dict(case=name,FEM_folder=folder.name,role=role,prediction=predictions[name])
        if not (folder/'run.json').exists():row['status']='pending_COMSOL'
        else:
            run=json.loads((folder/'run.json').read_text(encoding='utf-8'))
            result=load_result(folder)
            if run['status']!='solved' or not result['complete']:row['status']='failed_COMSOL'
            else:
                m=result['METRIC'];power=100*np.max(abs(m[4].real/m[3].real-1))
                f=dict(Rac=float(m[4,-1].real),Rdc=float(m[4,0].real),K=float(m[4,-1].real/m[4,0].real))
                error={key:100*(predictions[name][key]/f[key]-1) for key in f}
                row.update(COMSol=f,error_pct=error,power_error_pct=float(power))
                row['status']='pass' if max(abs(error['Rac']),abs(error['K']))<=5 and power<.1 else 'rejected'
        rows.append(row)
    rank=[]
    independent=[r for r in rows if r['role']=='independent' and 'COMSol' in r]
    for a,b in itertools.combinations(independent,2):
        contrast=(a['COMSol']['K']-b['COMSol']['K'])/max(a['COMSol']['K'],b['COMSol']['K'])
        predicted=a['prediction']['K']-b['prediction']['K']
        rank.append(dict(cases=[a['case'],b['case']],relative_K_difference=float(contrast),
            considered_separable=bool(abs(contrast)>.02),
            rank_agrees=bool(contrast*predicted>0)))
    accepted=all(r['status']=='pass' for r in rows) and all(r['rank_agrees'] for r in rank if r['considered_separable'])
    result=dict(status='screening_accepted' if accepted else 'not_accepted',source_sha256=source_digest(),
                tolerance_pct=5,fit_coefficients=None,cases=rows,ranking=rank,
                limitations='Curved-wire inductance remains an approximation. Screening acceptance is not final 3D verification of the optimum.')
    (ROOT/'rom_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result

if __name__=='__main__':print(json.dumps(check(),indent=2))
