"""Verify the formal-area baseline before freezing it or starting optimization."""
from pathlib import Path
from build_av_periodic import build_av
from export_field_samples import add_field_samples
from run_java import run
from summarize_pilots import audit

if __name__=='__main__':
    shared=dict(factors=(64,),pitches=(96.,),copper_area=6.,length_mm=.5,
                mesh_div=2,air_radius_mm=6.,gap_um=40.,mesh_kind='sweep',axial_layers=3,av_ooc=False)
    for name,changes in [
       ('baseline_64_av_m3',dict(mesh_div=3)),
       ('baseline_64_av_air10',dict(air_radius_mm=10.)),
       ('baseline_64_av_length1',dict(length_mm=1.,axial_layers=6)),
    ]:
        folder=Path(__file__).resolve().parents[1]/'data/raw'/name
        if (folder/'run.json').exists():raise RuntimeError('Refusing to overwrite recorded case '+name)
        p=build_av(name=name,**(shared|changes))
        if changes.get('mesh_div')==3:p=add_field_samples(p,nz=21)
        result=run(p);audit()
        if result['status']!='solved':break
