"""One-factor FEM checks of the Python-screened recommendation.

The primary acceptance limits remain those declared for the frozen baseline.
Reaction currents are an additional postprocessing diagnostic; they do not
change either field equations or the original volume/power measurements.
"""
from pathlib import Path
import json
from build_av_periodic import build_av
from export_field_samples import add_field_samples
from run_java import run
from summarize_pilots import audit

ROOT = Path(__file__).resolve().parents[1]


def add_reaction_currents(source):
    code = source.read_text()
    extra = '''
    try {
      m.result().numerical().create("reactionflux","IntSurface");
      m.result().numerical("reactionflux").set("expr",new String[]{"reacf(V)*I0/Itot"});
      for(int j=0;j<64;j++)for(String side:new String[]{"bottom","top"}) {
        m.result().numerical("reactionflux").selection().named(side+j);
        double[][] rxRealData=m.result().numerical("reactionflux").getReal();
        double[][] rxImagData=m.result().numerical("reactionflux").getImag();
        for(int s=0;s<rxRealData[0].length;s++)
          System.out.println("REACTION|"+j+"|"+side+"|"+s+"|"+rxRealData[0][s]+"|"+rxImagData[0][s]);
      }
    } catch(Exception e) {
      System.out.println("REACTION_UNAVAILABLE|"+e.getMessage());
    }
'''
    code = code.replace('    System.out.println("SOLVE_COMPLETE");',
                        extra+'    System.out.println("SOLVE_COMPLETE");')
    source.write_text(code)
    return source


if __name__ == '__main__':
    shared = dict(factors=(4,16),pitches=(-16.,128.),length_mm=4.,
                  copper_area=6.,mesh_div=2,air_radius_mm=6.,gap_um=40.,
                  zscale=.1,mesh_kind='tet',av_ooc=True,profile_frame='adaptive')
    cases = [
        ('recommended_m3_v2',dict(mesh_div=3)),
        ('recommended_length8',dict(length_mm=8.)),
        ('recommended_air10',dict(air_radius_mm=10.)),
    ]
    for name,changes in cases:
        folder=ROOT/'data/raw'/name
        if (folder/'run.json').exists():
            recorded=json.loads((folder/'run.json').read_text(encoding='utf-8'))
            if recorded.get('status')=='solved':
                print('Keeping completed verification: '+name,flush=True)
                continue
            raise RuntimeError('Refusing to overwrite '+name)
        source=build_av(name=name,**(shared|changes))
        if name=='recommended_m3_v2':
            source=add_field_samples(source,nz=31)
        source=add_reaction_currents(source)
        result=run(source)
        audit()
        if result['status']!='solved':
            break
