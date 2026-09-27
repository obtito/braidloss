"""Refine the two-level verification and run the remaining blind cases."""
from pathlib import Path
from build_av_periodic import build_av
from run_java import run
from summarize_pilots import audit

ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    shared=dict(copper_area=6.,mesh_div=2,zscale=.1,air_radius_mm=6.,gap_um=40.,
                mesh_kind='tet',av_ooc=True,profile_frame='curvature')
    cases=[
        ('validate_two_64_16_96_curvature_m3',dict(factors=(4,16),pitches=(16.,96.),length_mm=4.,mesh_div=3)),
        ('validate_three_64_32_32_96',dict(factors=(4,4,4),pitches=(32.,32.,96.),length_mm=8.)),
        ('holdout_three_64_m16_16_96',dict(factors=(4,4,4),pitches=(-16.,16.,96.),length_mm=4.)),
    ]
    for name,changes in cases:
        if (ROOT/'data/raw'/name/'run.json').exists():raise RuntimeError('Refusing to overwrite '+name)
        source=build_av(name=name,**(shared|changes))
        result=run(source);audit()
        if result['status']!='solved':break
