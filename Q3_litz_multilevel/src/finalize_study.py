"""Compile measured results, acceptance evidence, and matched scientific plots."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from summarize_pilots import load_result
from recursive_geometry import RecursiveCable
from visualize_results import field_dashboard

ROOT=Path(__file__).resolve().parents[1]
BASE='baseline_64_av_m3'
BEST='recommended_m325'


def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def measure(name):
    folder=ROOT/'data/raw'/name
    run=json.loads((folder/'run.json').read_text(encoding='utf-8'))
    raw=load_result(folder)
    if run['status']!='solved' or not raw['complete']:
        raise RuntimeError('Unsolved case '+name)
    if file_hash(folder/run['java'])!=run['sha256']:
        raise RuntimeError('Source differs from the recorded solver run: '+name)
    if not (folder/'model_fem.mph').is_file():
        raise RuntimeError('Missing saved FEM model: '+name)
    m=raw['METRIC']; current=raw['STRAND'][::3,-1]
    if not np.isfinite(m).all() or not np.isfinite(current).all():
        raise RuntimeError('Nonfinite measured result: '+name)
    if not np.allclose(m[0].real,[10.,200000.]):
        raise RuntimeError('Unexpected frequency convention: '+name)
    if abs(current.sum()-20)>1e-7 or not np.isclose(m[6,-1].real,400*m[4,-1].real):
        raise RuntimeError('RMS normalization or loss identity failed: '+name)
    cfg=json.loads((folder/'inputs.json').read_text(encoding='utf-8'))
    cut=raw.get('CUT')
    result=dict(case=name,Rdc=float(m[4,0].real),Rac=float(m[4,-1].real),
        K=float(m[4,-1].real/m[4,0].real),Pcu=float(m[6,-1].real),
        power_error_pct=float(100*np.max(abs(m[4].real/m[3].real-1))),
        current_sum_RMS_A=float(abs(current.sum())),
        current_amplitude_min_RMS_A=float(abs(current).min()),
        current_amplitude_max_RMS_A=float(abs(current).max()),
        current_amplitude_CV=float(abs(current).std()/abs(current).mean()),
        complex_current_imbalance_max=float(abs(current-current.sum()/len(current)).max()/abs(current.sum()/len(current))),
        source_sha256=run['sha256'],elapsed_seconds=run['elapsed_seconds'],
        currents_real_RMS_A=current.real.tolist(),currents_imag_RMS_A=current.imag.tolist())
    if cut is not None:
        cb=cut[::2,-1];ct=cut[1::2,-1]
        result['raw_gradient_cut_mismatch_pct_of_mean_amplitude']=float(100*max(abs(cb-ct))/abs(current).mean())
    lines=(folder/'stdout.log').read_text(encoding='utf-8',errors='replace').splitlines()
    rx={}
    for line in lines:
        if line.startswith('REACTION|'):
            _,j,side,s,re,im=line.split('|')
            if s=='1':rx[int(j),side]=complex(float(re),float(im))
    if len(rx)==2*len(current):
        bottom=np.array([rx[j,'bottom'] for j in range(len(current))])
        top=np.array([rx[j,'top'] for j in range(len(current))])
        perm=cfg['mapping']['strand_destination_to_source']
        result['reaction_flux_diagnostic']=dict(
            bottom_sum_real=float(bottom.sum().real),bottom_sum_imag=float(bottom.sum().imag),
            top_sum_real=float(top.sum().real),top_sum_imag=float(top.sum().imag),
            same_strand_conservation_pct=float(100*max(abs(bottom+top))/abs(current).mean()),
            screw_mapping_conservation_pct=float(100*max(abs(top+bottom[perm]))/abs(current).mean()),
            bottom_real=bottom.real.tolist(),bottom_imag=bottom.imag.tolist(),
            top_real=top.real.tolist(),top_imag=top.imag.tolist())
    result['evidence_sha256']={leaf:file_hash(folder/leaf) for leaf in
        ['inputs.json','run.json','stdout.log','model_fem.mph']}
    return result


def samples(name):
    folder=ROOT/'data/raw'/name
    geom=np.load(folder/'sampling_points.npz')
    lines=(folder/'stdout.log').read_text(encoding='utf-8',errors='replace').splitlines()
    values=np.array([float(x.split('|')[-1]) for x in lines if x.startswith('FIELD|')])
    if len(values)!=len(geom['points_m']) or not np.isfinite(values).all():
        raise RuntimeError('Incomplete field samples '+name)
    count=np.prod(geom['surface_m'].shape[:-1])
    return geom['surface_m']*1000,geom['cross_m']*1000,values[:count].reshape(geom['surface_m'].shape[:-1])/1e6,values[count:].reshape(geom['cross_m'].shape[:-1])/1e6


def plot_comparison(data, output_dir=None):
    output_dir=Path(output_dir) if output_dir else ROOT/'figures'
    output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,
        'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13.2,8.4),facecolor='#f7fafc')
    gs=fig.add_gridspec(1,2,left=.055,right=.91,top=.8,bottom=.28,wspace=.21)
    norm=Normalize(0,35);cmap=plt.get_cmap('turbo')
    for col,(name,title) in enumerate([(BASE,'基准：64 股单级，节距 96 mm'),(BEST,'推荐：4×16 两级，节距 −16 / 128 mm')]):
        _,cross,_,val=samples(name)
        ax=fig.add_subplot(gs[0,col]);ax.set_facecolor('white')
        for j in range(len(cross)):
            xy=np.vstack([cross[j,0,:1,:2],cross[j,1:,:,:2].reshape(-1,2)])
            v=np.r_[val[j,0,0],val[j,1:].ravel()]
            ax.tripcolor(mtri.Triangulation(xy[:,0],xy[:,1]),v,shading='gouraud',norm=norm,cmap=cmap,rasterized=True)
        ax.set_aspect('equal');ax.set(xlim=(-2.55,2.55),ylim=(-2.55,2.55),xlabel='x / mm',ylabel='y / mm',title=title)
        ax.grid(alpha=.1)
        row=data['baseline'] if col==0 else data['recommended']
        x=.065+col*.468
        fig.text(x,.213,f"Rac/Rdc  {row['K']:.4f}",fontsize=16,weight='bold',color='#142a43')
        fig.text(x,.17,f"交流电阻  {row['Rac']*1000:.4f} mΩ/m",fontsize=11,color='#334e68')
        fig.text(x,.128,f"铜损  {row['Pcu']:.4f} W/m",fontsize=14,weight='bold',color='#147d92')
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.932,.32,.017,.43]))
    cb.set_label('|J| RMS / (A/mm²)')
    fig.text(.055,.942,'COMSOL 实际求解：基准与推荐方案',fontsize=22,weight='bold',color='#142a43')
    fig.text(.055,.887,f"同样 6 mm²、64 股、200 kHz、20 A RMS；每米铜损降低 {data['improvement_pct']['Pcu']:.2f}%",fontsize=13,color='#147d92')
    fig.text(.055,.068,'两图使用相同尺寸与色标，截取各自周期单元的中部。各股电流由三维耦合场求解；未预设为 I/N。',fontsize=9,color='#526478')
    fig.text(.055,.031,'改善同时包含轨迹、外径及空隙率的影响。推荐限于本轮搜索范围；附近参数的微小差别不足以证明唯一最优。',fontsize=9,color='#526478')
    fig.savefig(output_dir/'final_comparison.png',dpi=175,facecolor=fig.get_facecolor());plt.close(fig)

    fig=plt.figure(figsize=(12.8,7.8),facecolor='#f7fafc')
    for col,(name,title) in enumerate([(BASE,'基准：实际求解单元 0.5 mm'),(BEST,'推荐：实际求解单元 4 mm')]):
        surface,_,values,_=samples(name)
        ax=fig.add_axes([.025+col*.445,.215,.405,.555],projection='3d',facecolor='#f7fafc')
        for j in range(len(surface)):
            pos=np.concatenate([surface[j],surface[j,:,:1]],axis=1)
            val=np.concatenate([values[j],values[j,:,:1]],axis=1)
            ax.plot_surface(pos[:,:,2],pos[:,:,0],pos[:,:,1],facecolors=cmap(norm(val)),
                rstride=1,cstride=1,linewidth=0,antialiased=True,shade=False)
        ax.set(xlabel='轴向 z / mm',ylabel='x / mm',zlabel='y / mm',title=title)
        ax.set_box_aspect((1.6,1,1));ax.view_init(24,-65)
        ax.tick_params(labelsize=8)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.936,.28,.016,.41]))
    cb.set_label('|J| RMS / (A/mm²)')
    fig.text(.055,.94,'真实三维铜实体上的电流密度',fontsize=22,weight='bold',color='#142a43')
    fig.text(.055,.875,'颜色来自 COMSOL 场解；显示铜半径 98% 处的采样面。不是几何颜色，也不是 Python 预测。',fontsize=11,color='#526478')
    fig.text(.055,.11,'周期单元沿导线重复，并按各自端面映射旋转、连接。三维视图调整轴向显示比例，长度请以坐标为准。',fontsize=10,color='#526478')
    fig.text(.055,.06,'局部颜色用于展示集肤和邻近效应；整体损耗采用全铜体积积分。点值与端面原始梯度通量的精度低于整体电阻。',fontsize=9,color='#526478')
    fig.savefig(output_dir/'final_3d_field.png',dpi=175,facecolor=fig.get_facecolor());plt.close(fig)


def finalize():
    names=['confirm_two_m16_128_m2','recommended_m3_v2',BEST,'recommended_length8','recommended_air10']
    checks={name:measure(name) for name in names}
    coarse=checks[names[0]]
    intermediate=checks['recommended_m3_v2']
    changes={
        'mesh_divisor_3_to_3.25':100*(checks[BEST]['Rac']/intermediate['Rac']-1),
        'cell_length_4_to_8_mm_same_mesh_settings':100*(checks['recommended_length8']['Rac']/coarse['Rac']-1),
        'air_radius_6_to_10_mm_same_mesh_settings':100*(checks['recommended_air10']['Rac']/coarse['Rac']-1)}
    mesh_history=[dict(from_divisor=2,to_divisor=3,
        Rac_change_pct=100*(intermediate['Rac']/coarse['Rac']-1),
        status='exceeded_1_pct_further_refinement_required'),
        dict(from_divisor=3,to_divisor=3.25,
             Rac_change_pct=changes['mesh_divisor_3_to_3.25'],status='check_against_unchanged_1_pct_limit')]
    if any(abs(v)>=1 for v in changes.values()):
        raise RuntimeError('Declared 1% Rac verification failed: '+str(changes))
    if any(row['power_error_pct']>=.1 for row in checks.values()):
        raise RuntimeError('Declared power check failed')
    geometry=RecursiveCable((4,16),(-16.,128.),copper_area_mm2=6.,gap_um=40.).diagnostics()
    search_records=[json.loads(p.read_text(encoding='utf-8')) for p in
                    (ROOT/'data/search_64').glob('*.json')]
    if len(search_records)!=129 or len({x['id'] for x in search_records})!=129:
        raise RuntimeError('The predeclared search record is incomplete')
    candidate=next(x for x in search_records if x['id']=='g4x16_pm16_128')
    dc_error=100*(checks[BEST]['Rdc']/geometry['Rdc_length_parallel_ohm_per_m']-1)
    if abs(dc_error)>=.5:raise RuntimeError('Declared DC check failed')
    baseline=measure(BASE);recommended=checks[BEST]
    frozen=json.loads((ROOT/'baseline_frozen.json').read_text(encoding='utf-8'))
    if baseline['source_sha256']!=frozen['results']['source_sha256']:
        raise RuntimeError('Frozen baseline source changed')
    for key in ['Rdc','Rac','K','Pcu']:
        if not np.isclose(baseline[key],frozen['results'][key],rtol=1e-10,atol=1e-12):
            raise RuntimeError('Frozen baseline value changed: '+key)
    current0=np.array(intermediate['currents_real_RMS_A'])+1j*np.array(intermediate['currents_imag_RMS_A'])
    current1=np.array(recommended['currents_real_RMS_A'])+1j*np.array(recommended['currents_imag_RMS_A'])
    improvement={k:100*(1-recommended[k]/baseline[k]) for k in ['Rac','K','Pcu']}
    cases=[]
    for name,id_ in [('confirm_two_m16_128_m2','g4x16_pm16_128'),
                     ('confirm_two_p16_128_m2','g4x16_p16_128'),
                     ('confirm_two_m32_128_m2','g4x16_pm32_128'),
                     ('confirm_three_m32_m16_64_m2','g4x4x4_pm32_m16_64')]:
        prediction=json.loads((ROOT/'data/search_64'/f'{id_}.json').read_text(encoding='utf-8'))
        actual=measure(name)
        cases.append(dict(candidate=id_,COMSOL=actual,ROM=prediction['ROM'],
            ROM_Rac_error_pct=100*(prediction['ROM']['Rac']/actual['Rac']-1),
            ROM_K_error_pct=100*(prediction['ROM']['K']/actual['K']-1)))
    exceeded=[row['candidate'] for row in cases
              if max(abs(row['ROM_Rac_error_pct']),abs(row['ROM_K_error_pct']))>5]
    assessment=dict(status='restricted_after_new_FEM_evidence' if exceeded else 'confirmed_on_selected_cases',
        threshold_pct=5,exceeded_candidates=exceeded,
        max_absolute_Rac_error_pct=max(abs(row['ROM_Rac_error_pct']) for row in cases),
        meaning='Historical screening retained. New designs require further validation; no reliable exhaustive ranking claim.',
        confirmed_cases=[dict(candidate=row['candidate'],Rac_error_pct=row['ROM_Rac_error_pct'],
                              K_error_pct=row['ROM_K_error_pct']) for row in cases])
    model=ROOT/'data/raw'/BEST/'model_fem.mph'
    data=dict(completed_utc=datetime.now(timezone.utc).isoformat(),
        status='recommended_design_passes_declared_integral_FEM_checks',
        scope='129 predeclared fixed-N=64 discrete candidates; screening and selected 3D confirmation; no global or unique optimum claim',
        baseline=baseline,recommended=recommended,recommended_geometry=geometry,
        search_counts=dict(Counter(x['status'] for x in search_records)),
        recommended_clearance=candidate['clearance'],
        post_search_ROM_assessment=assessment,
        improvement_pct=improvement,verification=dict(Rac_change_pct=changes,mesh_history=mesh_history,
            cumulative_mesh_2_to_3_25_Rac_change_pct=100*(recommended['Rac']/coarse['Rac']-1),
            interpretation='Successive changes satisfy declared checks; a 7.69% final h reduction is not a rigorous discretization-error bound.',
            dc_length_reference_error_pct=dc_error,
            volume_averaged_current_vector_mesh_change_pct=float(100*np.linalg.norm(current1-current0)/np.linalg.norm(current1)),
            accepted_limits=dict(Rac_change_pct=1,power_error_pct=.1,dc_error_pct=.5),cases=checks),
        confirmation_cases=cases,
        limitations=['N and diameter not searched; best inner and outer pitches touch search bounds',
            'Nearby opposite handedness is not uniquely distinguishable at the current numerical precision',
            'New three-level confirmation exceeds the initial 5% ROM target; ROM cannot certify all uncomputed candidates',
            'Current measurements are axial volume averages; raw boundary gradient fluxes are less accurate',
            'No thermal feedback, dielectric loss/capacitance, real leads/contacts, manufacturing costs or fatigue'],
        artifacts=dict(model=str(model.relative_to(ROOT)).replace('\\','/'),model_sha256=file_hash(model)))
    (ROOT/'final_results.json').write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf-8')
    (ROOT/'data/rom_post_confirmation.json').write_text(json.dumps(assessment,indent=2,ensure_ascii=False),encoding='utf-8')
    plot_comparison(data)
    field_dashboard(BEST,final=True)
    return data


if __name__=='__main__':
    result=finalize()
    print(json.dumps({k:result[k] for k in ['status','improvement_pct']},indent=2))
