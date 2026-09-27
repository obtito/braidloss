"""Independent full-area 3D validation cases, selected before ROM searching.

No empirical coefficients are fitted. These runs test different recursive
structures and an opposite-chirality holdout. Results must be reviewed before
the experimental impedance ROM is enabled for any design selection.
"""
from pathlib import Path
from build_av_periodic import build_av
from run_java import run
from summarize_pilots import audit

ROOT=Path(__file__).resolve().parents[1]

CASES=[
    ('validate_two_64_16_96',dict(factors=(4,16),pitches=(16.,96.),length_mm=4.)),
    ('validate_three_64_32_32_96',dict(factors=(4,4,4),pitches=(32.,32.,96.),length_mm=8.)),
    ('holdout_three_64_m16_16_96',dict(factors=(4,4,4),pitches=(-16.,16.,96.),length_mm=4.)),
    ('holdout_two_64_m32_96',dict(factors=(4,16),pitches=(-32.,96.),length_mm=8.)),
]

if __name__=='__main__':
    if not (ROOT/'baseline_frozen.json').exists():
        raise RuntimeError('The full-area baseline must be frozen first')
    shared=dict(copper_area=6.,mesh_div=2,zscale=.08,air_radius_mm=6.,
                gap_um=40.,mesh_kind='tet',av_ooc=True)
    for name,changes in CASES:
        folder=ROOT/'data/raw'/name
        if (folder/'run.json').exists():raise RuntimeError('Refusing to overwrite '+name)
        source=build_av(name=name,**(shared|changes))
        result=run(source);audit()
        if result['status']!='solved':break
