"""Readable Chinese report assembled only from retained numerical evidence.

Run with a Python environment containing ReportLab. No solver is launched.
"""
from pathlib import Path
from xml.sax.saxutils import escape
import json
import math

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak,
    KeepTogether,
)

ROOT=Path(__file__).resolve().parents[1]
INK=colors.HexColor('#19324b')
TEAL=colors.HexColor('#147d92')
GREY=colors.HexColor('#526478')
PALE=colors.HexColor('#eef4f8')
WIDTH=A4[0]-96


def build():
    pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\msyh.ttc',subfontIndex=0))
    pdfmetrics.registerFont(TTFont('CNBold',r'C:\Windows\Fonts\msyhbd.ttc',subfontIndex=0))
    pdfmetrics.registerFont(TTFont('Math',r'C:\Windows\Fonts\cambria.ttc',subfontIndex=0))
    pdfmetrics.registerFontFamily('CN',normal='CN',bold='CNBold',italic='CN',boldItalic='CNBold')
    styles={
        'title':ParagraphStyle('title',fontName='CNBold',fontSize=22,leading=31,textColor=INK,spaceAfter=13,wordWrap='CJK'),
        'h':ParagraphStyle('h',fontName='CNBold',fontSize=14,leading=21,textColor=INK,spaceBefore=12,spaceAfter=8,keepWithNext=True,wordWrap='CJK'),
        'body':ParagraphStyle('body',fontName='CN',fontSize=10,leading=16.2,textColor=INK,spaceAfter=8,wordWrap='CJK'),
        'small':ParagraphStyle('small',fontName='CN',fontSize=8.4,leading=13,textColor=GREY,spaceAfter=7,wordWrap='CJK'),
        'cell':ParagraphStyle('cell',fontName='CN',fontSize=8.6,leading=13,textColor=INK,wordWrap='CJK'),
        'headcell':ParagraphStyle('headcell',fontName='CNBold',fontSize=8.7,leading=13,textColor=colors.white,wordWrap='CJK'),
    }
    story=[]
    def para(text,style='body'):
        # YaHei omits nabla and several Unicode subscripts. Use a verified
        # embedded math font for those glyphs instead of emitting empty boxes.
        base=pdfmetrics.getFont('CN').face.charToGlyph
        math_chars=pdfmetrics.getFont('Math').face.charToGlyph
        for char in set(text):
            if ord(char)>32 and ord(char) not in base:
                if ord(char) not in math_chars:
                    raise ValueError('No embedded font covers '+repr(char))
                text=text.replace(char,'<font name="Math">'+char+'</font>')
        return Paragraph(text,styles[style])
    def p(text,style='body'):
        story.append(para(text,style))
    def heading(text):p(text,'h')
    def table(headers,rows,widths):
        data=[[para(escape(str(x)),'headcell') for x in headers]]
        data.extend([[para(escape(str(x)),'cell') for x in row] for row in rows])
        t=Table(data,colWidths=[WIDTH*x for x in widths],repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),INK),('VALIGN',(0,0),(-1,-1),'TOP'),
            ('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),
            ('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,PALE]),
            ('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#cad7df')),
        ]))
        story.extend([t,Spacer(1,9)])
    def figure(filename,caption,width=WIDTH):
        im=Image(str(ROOT/'figures'/filename))
        im.drawHeight=width*im.imageHeight/im.imageWidth
        im.drawWidth=width
        story.append(KeepTogether([im,Spacer(1,5),para(caption,'small')]))
    def page():story.append(PageBreak())
    def read(rel):return json.loads((ROOT/rel).read_text(encoding='utf-8'))

    d=read('final_results.json');b=d['baseline'];r=d['recommended'];g=d['recommended_geometry']
    v=d['verification'];frozen=read('baseline_frozen.json');bg=frozen['geometry']
    rom=read('rom_validation.json');post=d['post_search_ROM_assessment']
    coarse=v['cases']['confirm_two_m16_128_m2'];mid=v['cases']['recommended_m3_v2']
    delta=math.sqrt(1/(math.pi*200000*(4*math.pi*1e-7)*5.8e7))*1000

    p('传统多级绞合 Litz 线\n问题三续作报告'.replace('\n','<br/>'),'title')
    p('2026-09-27  ·  三维有限元验证与受控参数筛选','small')
    p(f'<b>推荐：4×16 两级绞合，内级节距 -16 mm，外级 +128 mm。</b>'
      f'在固定 64 股、6 mm² 铜面积、200 kHz、20 A RMS 下，三维 COMSOL 结果相对冻结基准的每米铜损下降 '
      f'<b>{d["improvement_pct"]["Pcu"]:.2f}%</b>。结论是已验证的结构改进，不能证明唯一或全局最优。')
    table(['指标','冻结基准：单级 64','推荐：两级 4×16'],[
        ['Rdc / (mΩ/m)',f'{b["Rdc"]*1000:.6f}',f'{r["Rdc"]*1000:.6f}'],
        ['Rac / (mΩ/m)',f'{b["Rac"]*1000:.6f}',f'{r["Rac"]*1000:.6f}'],
        ['Rac/Rdc',f'{b["K"]:.6f}',f'{r["K"]:.6f}'],
        ['铜损 / (W/m)',f'{b["Pcu"]:.6f}',f'{r["Pcu"]:.6f}'],
    ],[.36,.32,.32])
    figure('final_comparison.png','图 1  实际场解的截面对比。两图使用相同坐标范围与色标；所有电阻由全铜体积损耗积分提取。')
    p('本报告聚焦既有工作目录中的问题三。赛题规定的四题总论文、完美换位新拓扑与公开网站不属于本次已完成成果。','small')

    page();p('1  设计范围与结构定义','title')
    p(f'20 A / 6 mm² = 3.333 A/mm²，满足题目给定的名义直流电流密度上限 4 A/mm²。'
      f'由 Nπd²/4 = 6 mm² 得 d = {g["diameter_um"]/1000:.6f} mm。200 kHz 的铜趋肤深度为 {delta:.6f} mm；'
      '64 股是便于控制分组变量的结构研究起点，未证明其丝径最优。')
    table(['类型','内 → 外分组','主动变量'],[
        ['单级','64','外级节距 64、96、128 mm'],
        ['两级','4×16','内级 ±16、±32、±64 mm；外级同上'],
        ['三级','4×4×4','内、中级各 ±16、±32、±64 mm；外级同上'],
    ],[.17,.22,.61])
    p('本轮固定材料、电流、频率、64 股、单丝直径及 40 μm 名义间隙，枚举 3+18+108=129 个候选。'
      '绞角、丝长、外径、间隙和弯曲半径随轨迹联动计算，不作为彼此独立的输入。')
    figure('recommended_structure.png','图 2  先由 4 根绝缘单丝形成子束，再将 16 个子束绞合。颜色标识股线或子束，属于几何示意。')
    p('各级节距是在父束材料坐标中，相对转动一周对应的全局轴向参数增量。它不同于沿弯曲父轴测量的局部弧长节距。','small')

    page();p('2  场方程与统一计量','title')
    heading('递归轨迹与真实铜实体')
    p('Cchild(t) = Cparent(t) + r[cos φ(t)e₁(t) + sin φ(t)e₂(t)]，其中 φ(t)=2πt/p+φ₀。'
      '父轴法向基随曲线转动；每股铜丝截面垂直于自身局部轴线，随后扫掠成独立实体。'
      '总铜面积指各股法向面积之和。')
    table(['几何量','基准','推荐'],[
        ['外径保守上界 / mm',f'{bg["outer_diameter_bound_mm"]:.4f}',f'{g["outer_diameter_bound_mm"]:.4f}'],
        ['实际平均丝长增长',f'{100*(bg["growth_mean"]-1):.4f}%',f'{100*(g["growth_mean"]-1):.4f}%'],
        ['最大局部绞角',f'{bg["max_tangent_angle_deg"]:.3f}°',f'{g["max_tangent_angle_deg"]:.3f}°'],
        ['完整股线身份重复长度 / mm','96','128'],
        ['计算单元长度 / mm','0.5','4'],
    ],[.5,.25,.25])
    p(f'推荐的股间保守间隙下界 {d["recommended_clearance"]["clearance_lower_bound_um"]:.3f} μm，'
      f'最小弯曲半径 {g["min_bend_radius_mm"]:.3f} mm。间隙检查包含相邻周期样条副本；不代表制造公差认证。'
      f'推荐外径相对基准增加 {100*(g["outer_diameter_bound_mm"]/bg["outer_diameter_bound_mm"]-1):.1f}%，'
      '损耗变化同时包含分组、径向运动及空隙率的影响。')
    heading('磁准静态 A–φ 联合求解')
    p('J = σ[-∇φ - iωA + (V₀/L)e<sub>z</sub>]；∇·J=0；∇×(μ⁻¹∇×A)=J。'
      '使用 COMSOL 磁场与电流接口，静态形式 Coulomb 规范；铜丝间电绝缘。'
      '上下端的 A 和 φ 采用同一个经验证的旋转周期映射，每个股线映射闭环只设一个电势参考点。'
      '<b>各股电流由耦合场求得，没有给定为 I/N。</b>')
    p('外部磁绝缘圆柱半径统一为 6 mm，并另以 10 mm 检查边界影响。周期模型描述长导线内部，'
      '不包含实际引线、接头或端部损耗。推荐 4 mm 单元的旋转角为 11.25°；几何短周期并不等于完整股线身份重复周期。')
    heading('有效值与每米电阻')
    p('原始求解采用峰值复相量，Iraw=∫Jz dV/L。各股电流统一乘以 20/Iraw，场幅值乘以 20/|Iraw|，'
      '从而归一化到 20 A RMS。Rac=∫(|J|²/σ)dV/(|Iraw|²L)，Pcu=20²Rac。'
      '另用 Re(V₀/Iraw)/L 核对实功率；Rdc 取 10 Hz 极限，并与实际丝长并联公式比较。')

    page();p('3  数值检查与精度边界','title')
    p('原中断点位于补充网格加密之前。首次从 a/2 加密至 a/3 时，Rac 变化约 -1.0063%，'
      '略超预先设定的 1% 门槛；本次完成 a/3.25 的真实三维求解，并保留首次未达标的记录。')
    table(['推荐方案网格','Rac / (mΩ/m)','相对上一档'],[
        ['a/2',f'{coarse["Rac"]*1000:.6f}','初算参考'],
        ['a/3',f'{mid["Rac"]*1000:.6f}',f'{100*(mid["Rac"]/coarse["Rac"]-1):+.4f}%'],
        ['a/3.25',f'{r["Rac"]*1000:.6f}',f'{100*(r["Rac"]/mid["Rac"]-1):+.4f}%'],
    ],[.4,.3,.3])
    table(['独立检查','测得变化或误差','约定限值'],[
        ['单元长度 4 → 8 mm（a/2）',f'{v["Rac_change_pct"]["cell_length_4_to_8_mm_same_mesh_settings"]:+.4f}%','|ΔRac| < 1%'],
        ['空气半径 6 → 10 mm（a/2）',f'{v["Rac_change_pct"]["air_radius_6_to_10_mm_same_mesh_settings"]:+.4f}%','|ΔRac| < 1%'],
        ['最终模型输入功率与铜损',f'{r["power_error_pct"]:.6f}%','< 0.1%'],
        ['最终 Rdc 与丝长公式',f'{v["dc_length_reference_error_pct"]:+.6f}%','|误差| < 0.5%'],
        ['体积平均股电流向量网格变化',f'{v["volume_averaged_current_vector_mesh_change_pct"]:.4f}%','辅助诊断'],
    ],[.49,.28,.23])
    heading('这些检查支持什么')
    p(f'从 a/2 到 a/3.25，Rac 累计变化 {v["cumulative_mesh_2_to_3_25_Rac_change_pct"]:+.4f}%。'
      '最后一次单元尺寸仅减小 7.69%，相邻网格变化通过 1% 检查不等于已经证明真实离散误差小于 1%。'
      '空气与长度检查采用同一较粗网格作单因素比较，也不等同于更细网格下的联合误差评估。')
    p('冻结基准的网格、空气域及长度变化分别约 -0.6894%、+0.1132%、+0.0184%；'
      '基准输入功率与铜损差约 0.00590%。完整原始证据见 baseline_frozen.json。')
    heading('场点值与局部通量')
    p(f'推荐模型直接由电势梯度计算的两端股电流不平衡，仍达平均股电流幅值的 '
      f'{r["raw_gradient_cut_mismatch_pct_of_mean_amplitude"]:.3f}%。因此整体损耗的稳定程度不能直接套用到端面梯度或热点点值。')
    rx=r.get('reaction_flux_diagnostic')
    if rx:
        p(f'约束反力通量的同股守恒偏差为 {rx["same_strand_conservation_pct"]:.3g}%，'
          f'旋转周期配对偏差为 {rx["screw_mapping_conservation_pct"]:.3g}%。'
          '这说明离散约束满足良好，不构成对局部梯度精度或连续问题误差的独立证明。')

    page();p('4  参数筛选与独立复算','title')
    c=d['search_counts']
    p(f'129 组 Python 筛选中，{c["excluded_geometry"]} 组被几何条件排除，{c["excluded_loss"]} 组预测铜损高于基准，'
      f'{c["feasible_screening_candidate"]} 组进入可行筛选集。Python 使用物理多极近似，计入横向及轴向磁场响应，无拟合损耗系数。')
    labels=['两级 -16 / +128','两级 +16 / +128','两级 -32 / +128','三级 -32 / -16 / +64']
    table(['候选（节距 / mm）','COMSOL K','铜损 / (W/m)','Python Rac 偏差'],[
        [labels[i],f'{x["COMSOL"]["K"]:.5f}',f'{x["COMSOL"]["Pcu"]:.5f}',f'{x["ROM_Rac_error_pct"]:+.3f}%']
        for i,x in enumerate(d['confirmation_cases'])
    ],[.38,.18,.21,.23])
    p('表中四组统一使用候选筛选网格，故推荐行与首页最终细网格结果不同。'
      '三级候选只作初算比较，未通过最终推荐模型的全套收敛检查。','small')
    p(f'初始独立验证集的最大 Rac 偏差为 '
      f'{max(abs(x["error_pct"]["Rac"]) for x in rom["cases"] if x["role"]=="independent"):.2f}%，排序一致。'
      f'新增三级候选的偏差绝对值达到 {post["max_absolute_Rac_error_pct"]:.3f}%，超过原定 5% 目标。'
      '它同时含近似模型与筛选网格误差；本次保留该失败证据并限制模型后续使用，未调整系数以强行过关。')
    figure('parameter_controls.png','图 3  已保存的 Python 控制实验。小于模型误差的趋势只能作为待验证线索。')
    table(['两级内 / 外节距 / mm','Python K','丝长增长'],[
        [label,f'{read("data/search_64/"+name+".json")["ROM"]["K"]:.5f}',
         f'{100*(read("data/search_64/"+name+".json")["geometry"]["growth_mean"]-1):.4f}%']
        for name,label in [('g4x16_pm16_64','-16 / +64'),('g4x16_pm16_96','-16 / +96'),
                           ('g4x16_pm16_128','-16 / +128'),('g4x16_pm32_128','-32 / +128'),
                           ('g4x16_pm64_128','-64 / +128')]
    ],[.48,.26,.26])
    p('外级节距增大时，本组离散点的预测 K 与丝长均下降；内级节距影响更小。推荐参数触及搜索边界，'
      '不能据此证明连续最优。正反绞向候选的初算 K 仅相差约 0.17%，当前精度不足以可靠区分。','small')

    page();p('5  真实场图与电流分布','title')
    figure('final_3d_field.png','图 4  场值直接来自 COMSOL 解，在铜半径 98% 处取样。各模型实际单元长度不同，三维显示的轴向比例经过调整。')
    table(['股电流指标（20 A RMS）','基准','推荐'],[
        ['最小股电流幅值 / A',f'{b["current_amplitude_min_RMS_A"]:.6f}',f'{r["current_amplitude_min_RMS_A"]:.6f}'],
        ['最大股电流幅值 / A',f'{b["current_amplitude_max_RMS_A"]:.6f}',f'{r["current_amplitude_max_RMS_A"]:.6f}'],
        ['股电流幅值变异系数',f'{b["current_amplitude_CV"]:.5f}',f'{r["current_amplitude_CV"]:.5f}'],
        ['最大复数均流偏差',f'{100*b["complex_current_imbalance_max"]:.2f}%',f'{100*r["complex_current_imbalance_max"]:.2f}%'],
    ],[.5,.25,.25])
    p('复数均流偏差定义为 max|Iⱼ-I总/N| / |I总/N|，同时含幅值与相位，可能超过 100%。'
      '各股净电流取轴向体积平均值；局部 |J| 云图与每股净电流代表不同物理量。均流参考值 20/64=0.3125 A 仅为比较基准，未用于施加各股激励。')
    p('同一股内部仍可存在趋肤与邻近涡流。推荐损耗降低并不意味着处处电流均匀，'
      '也不意味着已经得到问题四所要求的完整径向遍历或“完美换位”。')

    page();p('6  证据索引与复现','title')
    p('工作成果保存在用户指定的 litz_multilevel 文件夹。目录中的失败试算与早期几何探索继续保留，'
      '以 final_results.json、更新后的进度文档及本报告为本轮结论入口。')
    table(['文件 / 目录','用途'],[
        ['final_results.json','最终指标、独立检查、源与模型 SHA-256'],
        ['baseline_frozen.json','优化前冻结的正式单级基准'],
        ['parameter_search.csv；data/search_64/','129 组参数、预测、几何及排除原因'],
        ['data/raw/baseline_64_av_m3/','基准 MPH、Java、输入、原始结果和资源日志'],
        ['data/raw/recommended_m325/','本次补充加密的推荐模型与完整原始证据'],
        ['data/rom_post_confirmation.json','追加验证超标记录与后续模型使用限制'],
        ['MODEL_METHOD.md；DESIGN_SPACE.md','场方程、几何定义、先验搜索范围'],
        ['src/；figures/','可复现程序与实际场图、几何图'],
    ],[.49,.51])
    heading('本机工具与复现步骤')
    p('COMSOL：D:/tools/comosol/COMSOL62/Multiphysics/bin/win64。'
      '科学计算 Python：D:/tools/anaconda/python.exe。主要库为 NumPy、SciPy、Matplotlib、psutil、threadpoolctl；'
      'PDF 生成使用 ReportLab 和微软雅黑字体。')
    p('在研究目录依次运行 python src/finalize_study.py、python src/write_report.py，'
      '可由已有解重建结果与图文。python src/build_report_pdf.py 生成本 PDF；'
      'python src/audit_delivery.py 核验数据与证据。若环境未安装 ReportLab，可使用 Codex 随附的文档 Python 运行 PDF 脚本。')
    p('用 COMSOL 打开 model_fem.mph 可查看几何、网格、场与设置。需要重新求解时，从已保存输入生成新案例名，'
      '再调用 src/run_java.py；不覆盖已记录的案例。最多使用 4 个求解线程，连续低内存时只终止本次计算。')
    heading('范围与后续研究')
    p('本次完成固定 64 股的既定研究闭环。未优化股数与联动丝径，未进行扫频、真实温升反馈、介质损耗、电容、接头、成本或疲劳验证。'
      '若拓展节距或股数范围，应重新建立控制组并补充三维验证；整套赛题的问题一、二、四和公开网站需在其他工作中完成。')
    heading('AI 工具使用说明')
    p('本次续作由 Codex（当前会话基于 GPT-6）辅助阅读、遗留资料审查、程序修订、计算续跑、数据核验及报告制作。'
      '数值由实际 Python 与 COMSOL 运行产生。历史目录未完整保存此前 AI 模型标识，故不补写未经证实的版本；总论文须合并队伍完整使用记录。','small')

    out=ROOT/'问题三_续作验证报告.pdf'
    def decoration(canvas,doc):
        canvas.saveState()
        canvas.setFillColor(TEAL);canvas.rect(48,A4[1]-35,WIDTH,2,fill=1,stroke=0)
        canvas.setFont('CN',8);canvas.setFillColor(GREY)
        canvas.drawString(48,A4[1]-25,'LITZ  /  问题三 · 受控研究与三维验证')
        canvas.drawString(48,27,'6 mm²  ·  64 股  ·  200 kHz  ·  20 A RMS')
        canvas.drawRightString(A4[0]-48,27,str(doc.page))
        canvas.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,leftMargin=48,rightMargin=48,
        topMargin=53,bottomMargin=47,title='问题三：传统多级绞合 Litz 线续作验证报告',
        author='研究工作记录 / Codex 辅助整理')
    doc.build(story,onFirstPage=decoration,onLaterPages=decoration)
    return out


if __name__=='__main__':
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    print(build())
