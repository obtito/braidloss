"""Actual 3D confirmation of the screened winner and nearby alternatives."""
from pathlib import Path
import json
from build_av_periodic import build_av
from export_field_samples import add_field_samples
from run_java import run
from summarize_pilots import audit

ROOT=Path(__file__).resolve().parents[1]
CASES=[
    ('confirm_two_m16_128_m2',dict(factors=(4,16),pitches=(-16.,128.),length_mm=4.,zscale=.1)),
    ('confirm_two_p16_128_m2',dict(factors=(4,16),pitches=(16.,128.),length_mm=4.,zscale=.1)),
    ('confirm_two_m32_128_m2',dict(factors=(4,16),pitches=(-32.,128.),length_mm=8.,zscale=.05)),
    ('confirm_three_m32_m16_64_m2',dict(factors=(4,4,4),pitches=(-32.,-16.,64.),length_mm=8.,zscale=.1)),
]

if __name__=='__main__':
    search=json.loads((ROOT/'search_result.json').read_text(encoding='utf-8'))
    if search['best_candidate']['id']!='g4x16_pm16_128':raise RuntimeError('Review candidate list after a changed search')
    shared=dict(copper_area=6.,mesh_div=2,air_radius_mm=6.,gap_um=40.,
                mesh_kind='tet',av_ooc=True,profile_frame='adaptive')
    for name,changes in CASES:
        folder=ROOT/'data/raw'/name
        if (folder/'run.json').exists():raise RuntimeError('Refusing to overwrite '+name)
        source=build_av(name=name,**(shared|changes))
        if name==CASES[0][0]:source=add_field_samples(source,nz=25)
        result=run(source);audit()
        if result['status']!='solved':break
