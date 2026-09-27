"""Re-render existing geometry; write only this request's explanatory artifacts."""
from pathlib import Path
import json, importlib.util, urllib.request
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE=Path(__file__).resolve().parent
DATA=json.loads((HERE.parent/'render_data.json').read_text(encoding='utf-8'))
VIS=Path('C:/Users/zq257/.codex/visualizations/2026/09/26/01a0dc56-200c-71b3-913a-f642f822c58b')
SKILL=Path('C:/Users/zq257/.codex/plugins/cache/openai-bundled/visualize/1.0.39/skills/visualize')

# Preserve the pinned library locally, including its original license notice.
# The result does not require the CDN when the user opens the model.
lib=HERE/'three-0.128.0.min.js'
if not lib.exists():
    request=urllib.request.Request('https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js',headers={'User-Agent':'Mozilla/5.0'})
    lib.write_bytes(urllib.request.urlopen(request,timeout=30).read())
source=(VIS/'litz-transposition-volume.html').read_text(encoding='utf-8')
tag='<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>'
assert tag in source
source=source.replace(tag,'<script>\n'+lib.read_text(encoding='utf-8').replace('</script','<\\/script')+'\n</script>')
source=source.replace("let mode='principle', phase=0", "let mode='bundle', phase=0")
fragment=VIS/'litz-transposition-rerender.html'
fragment.write_text(source,encoding='utf-8')
assert fragment.stat().st_size<1_000_000
spec=importlib.util.spec_from_file_location('viz_render',SKILL/'scripts/render.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.export_html(fragment,HERE/'litz-interactive.html',title='Litz 完整换位：可旋转三维模型',force=True)

# A new scientific figure, rendered directly from the same spline coefficients.
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei'],'axes.unicode_minus':False,
                     'font.size':13,'savefig.facecolor':'white'})
C=np.array(DATA['coeff16'])
def xy(u):
    u=np.atleast_1d(u)%16;i=np.floor(u).astype(int);t=(u-i)[:,None];a=C[i]
    return ((a[:,3]*t+a[:,2])*t+a[:,1])*t+a[:,0]
blue='#2478b4';orange='#d88231';grey='#b5bdc6';ink='#253440'
fig,axes=plt.subplots(2,2,figsize=(11,10.5),dpi=180)
fig.subplots_adjust(left=.03,right=.97,bottom=.09,top=.88,wspace=.07,hspace=.25)
fig.suptitle('同一根单丝，依次从内部走向边缘，再返回',fontsize=21,x=.5,y=.972,color=ink)
fig.text(.5,.925,'16 根缩例｜四幅图对应不同轴向位置｜蓝色始终是单丝 A',ha='center',fontsize=13,color='#536273')
labels=['① 起点：A 在内部','② 走过 1/4 圈：A 到边缘','③ 走过 1/2 圈：A 到角部','④ 完整一圈：A 回到原位']
for ax,phase,title in zip(axes.flat,[0,4,8,16],labels):
    route=xy(np.linspace(0,16,641));ax.plot(route[:,0],route[:,1],color='#8c99a8',lw=1.3,zorder=1)
    if phase:
        trace=xy(np.linspace(5,5+phase,phase*50+1));ax.plot(trace[:,0],trace[:,1],color=blue,lw=3,alpha=.75,zorder=2)
    for u in [1.3,4.3,8.3,12.3]:
        a,b=xy([u,u+.45]);ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=16,color='#637283',lw=1.2,zorder=3))
    for j in range(16):
        p=xy(j)[0];ax.add_patch(Circle(p,.26,fill=False,color='#dae0e6',lw=1))
        ax.text(p[0]+.29,p[1]+.29,str(j+1),fontsize=11,ha='center',va='center',color='#566576')
    for j in [i for i in range(16) if i not in [5,13]]+[13,5]:
        p=xy(j+phase)[0];fc=blue if j==5 else orange if j==13 else grey
        ax.add_patch(Circle(p,.225,facecolor=fc,edgecolor='white',lw=1.2,zorder=4))
        if j in [5,13]:ax.text(p[0],p[1],'A' if j==5 else 'B',ha='center',va='center',fontsize=14,color='white',weight='bold',zorder=5)
    ax.set_title(title,fontsize=15,pad=12,color=ink)
    ax.set_xlim(-2.02,2.02);ax.set_ylim(-1.93,2.02);ax.set_aspect('equal');ax.axis('off')
fig.text(.5,.047,'沿闭合路线连续前进；编号随单丝移动。此图为原理缩例，非 484 根实际尺寸截面。',ha='center',fontsize=12,color='#536273')
fig.savefig(HERE/'01-transposition-steps.png',dpi=180)
fig.savefig(HERE/'01-transposition-steps.svg')
plt.close(fig)
print(json.dumps({'fragment_bytes':fragment.stat().st_size,'interactive':str(HERE/'litz-interactive.html'),'static_steps':str(HERE/'01-transposition-steps.png')},ensure_ascii=False))
