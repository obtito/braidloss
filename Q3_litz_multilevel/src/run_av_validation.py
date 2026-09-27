"""Sequential pilot checks. No optimization or final acceptance is inferred."""
from pathlib import Path
from build_av_periodic import build_av
from run_java import run
from summarize_pilots import audit

if __name__=='__main__':
    cases=[
        ('pilot_9_av_m3',dict(mesh_div=3,av_ooc=False)),
        ('pilot_9_av_doublecell',dict(mesh_div=2,length_mm=8/3,av_ooc=False)),
        ('pilot_9_av_air15',dict(mesh_div=2,air_radius_mm=1.5,av_ooc=False)),
    ]
    for name,kw in cases:
        p=Path(__file__).resolve().parents[1]/'data/raw'/name
        if (p/'run.json').exists():raise RuntimeError('Refusing to overwrite a recorded case: '+name)
        result=run(build_av(name=name,**kw))
        audit()
        if result['status']!='solved':break
