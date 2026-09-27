"""Read raw solver evidence; solver completion is not physical acceptance."""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.special import jv

ROOT=Path(__file__).resolve().parents[1]


def load_result(folder):
    blocks={}
    log=folder/'stdout.log'
    text=log.read_text(encoding='utf-8',errors='replace') if log.exists() else ''
    for tag in ['METRIC','STRAND','CUT']:
        rows=[line.split('|') for line in text.splitlines() if line.startswith(tag+'|')]
        if not rows:continue
        a=np.zeros((max(int(r[1]) for r in rows)+1,max(int(r[2]) for r in rows)+1),complex)
        for _,i,j,re,im in rows:a[int(i),int(j)]=complex(float(re),float(im))
        blocks[tag]=a
    blocks['complete']='SOLVE_COMPLETE' in text
    return blocks


def audit():
    rows=[]
    frozen=json.loads((ROOT/'baseline_frozen.json').read_text(encoding='utf-8')) if (ROOT/'baseline_frozen.json').exists() else None
    final=json.loads((ROOT/'final_results.json').read_text(encoding='utf-8')) if (ROOT/'final_results.json').exists() else None
    for p in sorted((ROOT/'data/raw').iterdir()):
        if not (p/'run.json').exists():continue
        r=json.loads((p/'run.json').read_text(encoding='utf-8'))
        row=dict(case=p.name,solver_status=r['status'],physical_status='not_accepted',
          time_s=r.get('elapsed_seconds'),peak_rss_GiB=r.get('peak_family_rss_GiB'),
          peak_private_GiB=r.get('peak_family_private_GiB'),minimum_available_GiB=r.get('minimum_available_GiB'))
        a=load_result(p)
        if a.get('complete') and 'METRIC' in a:
            m=a['METRIC'];single='pilot_single' in p.name
            if single:
                radius=np.sqrt(6e-6/343/np.pi)
                frequency=m[0].real
                k=(1-1j)*np.sqrt(np.pi*frequency*4*np.pi*1e-7*5.8e7)
                z=k*jv(0,k*radius)/(2*np.pi*radius*5.8e7*jv(1,k*radius))
                err=100*(m[3].real/z.real-1)
                row['analytic_Rac_max_error_pct']=float(abs(err).max())
                row['physical_status']='single_wire_analytic_pass' if abs(err).max()<1 else 'rejected_analytic_error'
                row['power_error_pct']=float(100*np.max(abs(m[4].real/(.5*np.real(m[2]*m[1].conj())/.001)-1)))
            else:
                err=100*(m[4].real/m[3].real-1)
                row.update(power_error_pct=float(abs(err).max()),Rdc_limit_ohm_per_m=float(m[4,0].real),
                           Rac_ohm_per_m=float(m[4,-1].real),K=float(m[4,-1].real/m[4,0].real))
                row['physical_status']='power_pass_other_checks_pending' if abs(err).max()<.1 else 'rejected_power_balance'
                if p.name.startswith('pilot_7_helix'):
                    row['physical_status']='provisional_end_geometry_not_used_as_reference'
                if p.name in ['pilot_7_round_periodic_sign']:
                    row['physical_status']='power_pass_rejected_twisted_extrude_area_error'
                if (p/'inputs.json').exists() and 'STRAND' in a:
                    cfg=json.loads((p/'inputs.json').read_text(encoding='utf-8'))
                    perm=cfg.get('mapping',{}).get('strand_destination_to_source')
                    if perm is not None:
                        currents=a['STRAND'][::3]
                        mismatch=abs(currents-currents[perm]).max(axis=0)/abs(currents).mean(axis=0)
                        row['current_identity_mapping_error_pct']=float(100*mismatch.max())
                        if mismatch.max()>.001:
                            row['physical_status']='rejected_periodic_current_continuity'
                        if 'CUT' in a:
                            cb=a['CUT'][::2];ct=a['CUT'][1::2]
                            row['cut_same_strand_error_pct']=float(100*abs(cb-ct).max()/abs(currents).mean())
                            row['cut_periodic_mapping_error_pct']=float(100*abs(ct-cb[perm]).max()/abs(currents).mean())
        if frozen and p.name==frozen['case'] and r.get('sha256')==frozen['results']['source_sha256']:
            row['physical_status']='accepted_formal_baseline'
        if final and p.name==final['recommended']['case'] and r.get('sha256')==final['recommended']['source_sha256']:
            row['physical_status']='accepted_recommendation_integral_checks'
        rows.append(row)
        (p/'validation.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (ROOT/'pilot_evidence.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    (ROOT/'pilot_evidence.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    return rows


def field_figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,'font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    folder=ROOT/'data/raw/pilot_9_recursive_medium_fields'
    raw=(folder/'stdout.log').read_text(encoding='utf-8',errors='replace')
    rows=np.array([[float(x) for x in line.split('|')[1:]] for line in raw.splitlines() if line.startswith('FIELD|')])
    q=np.load(folder/'sampling_points.npz')['points']
    if len(rows)!=np.prod(q.shape[:-1]) or not np.isfinite(rows).all():raise ValueError('Field export is incomplete or invalid')
    value=rows[:,-1].reshape(q.shape[:-1])/1e6
    m=load_result(folder)['METRIC'];currents=load_result(folder)['STRAND'][::3,-1]
    fig=plt.figure(figsize=(13,7.7),facecolor='#f7fafc')
    ax=fig.add_axes([.035,.16,.51,.67],projection='3d',facecolor='#f7fafc')
    norm=Normalize(1.5,5.5);cmap=plt.get_cmap('turbo')
    for j in range(9):
        pos=np.concatenate([q[j],q[j,:,:1]],axis=1)*1000
        v=np.concatenate([value[j],value[j,:,:1]],axis=1)
        ax.plot_surface(pos[:,:,2],pos[:,:,0],pos[:,:,1],facecolors=cmap(norm(v)),rstride=1,cstride=1,
                        linewidth=0,antialiased=True,shade=False)
    ax.set(xlabel='轴向 z / mm',ylabel='x / mm',zlabel='y / mm',title='COMSOL 求解的铜内电流密度有效值')
    ax.set_box_aspect((2.1,1,1));ax.view_init(25,-65)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.03,pad=.06,shrink=.7)
    cb.set_label('|J| RMS / (A/mm²)')
    bx=fig.add_axes([.65,.48,.3,.32],facecolor='white')
    x=np.arange(1,10);colors=plt.get_cmap('tab10')(np.repeat(np.arange(3),3))
    bx.bar(x,abs(currents)*1000,color=colors,width=.65)
    bx.axhline(abs(currents.sum())/9*1000,ls='--',color='#273f57',lw=1,label='总电流幅值 ÷ 9，仅作参考')
    bx.set(xlabel='股线编号',ylabel='电流幅值 RMS / mA',title='各股电流由耦合场求解')
    bx.set_xticks(x);bx.legend(fontsize=8,loc='lower right');bx.grid(axis='y',alpha=.15)
    fig.text(.645,.4,'有限样品试算参数',weight='bold',fontsize=13,color='#142a43')
    info=[f'3×3 两级绞合；单丝直径 149.24 μm',
          f'轴向长度 1.333 mm；节距 4 / 12 mm',
          f'铜面积 0.1574 mm²；总电流 {20*9/343:.4f} A RMS',
          f'200 kHz；Rac/Rdc ≈ {m[4,-1].real/m[4,0].real:.4f}',
          f'输入功率与铜损相差 {abs(m[4,-1].real/m[3,-1].real-1)*100:.5f}%']
    for k,line in enumerate(info):fig.text(.645,.348-k*.045,line,fontsize=10,color='#334e68')
    fig.text(.055,.94,'3×3 多级绞合：真实三维实体场试算',fontsize=22,weight='bold',color='#142a43')
    fig.text(.055,.886,'颜色来自 COMSOL 数值结果；未给各股强制分配 I/N。此图用于验证建模方法。',fontsize=11,color='#526478')
    fig.text(.055,.065,'边界：有限样品两端共同接线，外侧采用磁绝缘。此结果不代表 343 股正式基准，也不能直接当作无限长周期导线的参数。',fontsize=9,color='#526478')
    fig.text(.055,.032,'表面颜色在 0.99 倍铜半径处采样；三维显示调整了坐标轴比例。直流参考取 10 Hz 极限，网格验证尚在进行。',fontsize=9,color='#526478')
    folder=ROOT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/'pilot_9_actual_3d_field.png',dpi=180,facecolor=fig.get_facecolor())
    plt.close(fig)
    np.savez_compressed(ROOT/'data/pilot_9_field.npz',coordinates_m=q,J_rms_A_per_m2=value*1e6,currents_RMS_A=currents)


if __name__=='__main__':
    rows=audit();field_figure()
    print(json.dumps([r for r in rows if r['physical_status']!='not_accepted'],indent=2))
