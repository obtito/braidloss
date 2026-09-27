"""Build real round swept copper solids about recursively stranded axes."""
from pathlib import Path
import json
import numpy as np
from recursive_geometry import RecursiveCable
from run_java import run

ROOT=Path(__file__).resolve().parents[1]


def build(name="pilot_9_recursive_v1", length_mm=1.33333333333333, layers=12):
    c=RecursiveCable((3,3),(4.,12.),copper_area_mm2=6*9/343)
    length=length_mm*1e-3
    pts,_,_=c.at_physical_z(np.linspace(-2*c.a,length+2*c.a,97))
    paths=[','.join(format(v,'.16g') for v in path.ravel()) for path in pts]
    code=(ROOT/'data/raw/pilot_7_helix_v2/PilotHelix.java').read_text()
    code=code.replace('PilotHelix','PilotRecursive')
    code=code.replace('public class PilotRecursive {','public class PilotRecursive {\n  static String[] paths={'+','.join('"'+s+'"' for s in paths)+'};')
    code=code.replace('m.param().set("L", "0.5[mm]");',f'm.param().set("L", "{length_mm}[mm]");')
    code=code.replace('"20[A]*7/343"','"20[A]*9/343"')
    code=code.replace('"2*a+16[um]"',f'"{c.envelope_radius:.16g}[m]"')
    i=code.index('    g.create("cu0", "Cylinder");');j=code.index('    g.create("air", "Cylinder");',i)
    geo='''    for(int j=0;j<paths.length;j++) {
      String wp="wp"+j, ic="ic"+j, sw="sw"+j, clip="clip"+j, cu="cu"+j;
      String[] raw=paths[j].split(",");
      double[][] pts=new double[raw.length/3][3];
      for(int k=0;k<raw.length;k++)pts[k/3][k%3]=Double.parseDouble(raw[k]);
      g.create(wp,"WorkPlane");
      g.feature(wp).geom().create("circle","Circle");
      g.feature(wp).geom().feature("circle").set("r","a");
      g.run(wp);
      g.create(ic,"InterpolationCurve");
      g.feature(ic).set("table",pts);
      g.feature(ic).set("rtol",1e-7);
      g.run(ic);
      String curveObj="";
      for(String ob:g.objectNames())if(ob.startsWith(ic))curveObj=ob;
      System.out.println("STAGE|curve|"+j+"|"+curveObj);
      g.create(sw,"Sweep");
      g.feature(sw).selection("enttosweep").init(2);
      g.feature(sw).selection("enttosweep").set(wp+".circle",new int[]{1});
      g.feature(sw).selection("edge").set(curveObj,new int[]{1});
      g.feature(sw).set("movetospine",true);
      g.feature(sw).set("rtol",1e-6);
      g.create(clip,"Block");
      g.feature(clip).set("size",new String[]{"2*Rair","2*Rair","L"});
      g.feature(clip).set("pos",new String[]{"-Rair","-Rair","0"});
      g.create(cu,"Intersection");
      g.feature(cu).selection("input").set(new String[]{sw,clip});
      g.feature(cu).set("selresult",true);
      g.run(cu);
      System.out.println("STAGE|solid|"+j);
    }
'''
    code=code[:i]+geo+code[j:]
    code=code.replace('new String[7]','new String[9]').replace('j<7','j<9')
    code=code.replace('mf.ICoil_6+mf.ICoil_7','mf.ICoil_6+mf.ICoil_7+mf.ICoil_8+mf.ICoil_9')
    code=code.replace('new String[21]','new String[27]')
    code=code.replace('set("numelem",6)',f'set("numelem",{layers})')
    code=code.replace('3D solid seven-wire helix pilot - common voltage','3D round swept solids - recursive 3 x 3 stranding pilot')
    folder=ROOT/'data/raw'/name;folder.mkdir(parents=True,exist_ok=True)
    java=folder/'PilotRecursive.java';java.write_text(code)
    np.savez_compressed(folder/'paths.npz',xyz=pts,normal_radius_m=c.a)
    info=c.diagnostics();info.update(formal_baseline=False,axial_sample_length_mm=length_mm,
      geometry_type="Round normal-plane sweeps, extended and then cut at common physical z planes",phase_convention="Peak COMSOL phasors; report normalized RMS results")
    (folder/'inputs.json').write_text(json.dumps(info,indent=2))
    return java


if __name__=='__main__':
    run(build())
