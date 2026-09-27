"""Additional refinement after the 2-to-3 mesh change exceeded 1%.

Only the copper/face size divisor changes. 3.25 is a further 7.69% reduction
in maximum transverse element size, chosen within this machine's RAM budget.
The original 1.0063% failed check remains in the final convergence history.
"""
from pathlib import Path
from build_av_periodic import build_av
from export_field_samples import add_field_samples
from verify_recommended import add_reaction_currents
from run_java import run
from summarize_pilots import audit

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    name='recommended_m325'
    if (root/'data/raw'/name/'run.json').exists():
        raise RuntimeError('Refusing to overwrite recorded refinement')
    source=build_av(name=name,factors=(4,16),pitches=(-16.,128.),length_mm=4.,
        copper_area=6.,mesh_div=3.25,air_radius_mm=6.,gap_um=40.,zscale=.1,
        mesh_kind='tet',av_ooc=True,profile_frame='adaptive')
    source=add_field_samples(source,nz=21)
    source=add_reaction_currents(source)
    result=run(source)
    audit()
