"""Seven round helical solids with explicitly paired periodic end meshes."""
from pathlib import Path
from run_java import run

ROOT=Path(__file__).resolve().parents[1]


def build(name="pilot_7_round_periodic_v1"):
    code=(ROOT/'data/raw/pilot_7_periodic_conforming_v2/PilotHelix.java').read_text()
    start=code.index('    g.create("cu0", "Cylinder");')
    stop=code.index('    g.create("wpair","WorkPlane");',start)
    geo='''    g.create("wpc","WorkPlane");
    g.feature("wpc").geom().create("circle","Circle");
    g.feature("wpc").geom().feature("circle").set("r","a");
    g.run("wpc");
    g.create("cu0","Extrude");
    g.feature("cu0").selection("input").set(new String[]{"wpc"});
    g.feature("cu0").set("distance","L");
    g.feature("cu0").set("twist","-360[deg]*L/pitch");
    g.feature("cu0").set("selresult",true);
    for(int j=1;j<7;j++) {
      String h="hx"+j, clip="clip"+j, cu="cu"+j;
      g.create(h,"Helix");
      g.feature(h).set("rmin","a");
      g.feature(h).set("rmaj","spacing");
      g.feature(h).set("axialpitch","pitch");
      g.feature(h).set("turns","(L+4*a)/pitch");
      g.feature(h).set("rot",(60*(j-1))+"[deg]-720[deg]*a/pitch");
      g.feature(h).set("pos",new String[]{"0","0","-2*a"});
      g.feature(h).set("endcaps","perpspine");
      g.feature(h).set("twistcomp",false);
      g.feature(h).set("rtol",1e-7);
      g.create(clip,"Block");
      g.feature(clip).set("size",new String[]{"2*Rair","2*Rair","L"});
      g.feature(clip).set("pos",new String[]{"-Rair","-Rair","0"});
      g.create(cu,"Intersection");
      g.feature(cu).selection("input").set(new String[]{h,clip});
      g.feature(cu).set("selresult",true);
    }
'''
    code=code[:start]+geo+code[stop:]
    # COMSOL Extrude's positive twist is clockwise; native Helix right is CCW.
    # The Rotated coordinate system for the destination uses the same passive convention.
    code=code.replace('g.feature("air").set("twist","360[deg]*L/pitch")','g.feature("air").set("twist","-360[deg]*L/pitch")')
    code=code.replace('new String[]{"360[deg]*L/pitch","0","0"}','new String[]{"-360[deg]*L/pitch","0","0"}')
    i=code.index('    m.component("c").mesh().create("mesh");')
    setup='''    java.util.Set<Integer> bcu=new java.util.HashSet<Integer>();
    java.util.Set<Integer> tcu=new java.util.HashSet<Integer>();
    for(int j=0;j<7;j++)for(String side:new String[]{"bottom","top"}) {
      String tag=side+j;
      m.component("c").selection().create(tag,"Intersection");
      m.component("c").selection(tag).set("entitydim",2);
      m.component("c").selection(tag).set("input",new String[]{side,"g_cu"+j+"_bnd"});
      for(int e:m.component("c").selection(tag).entities()) {
        if(side.equals("bottom"))bcu.add(e);else tcu.add(e);
      }
    }
    for(String side:new String[]{"bottom","top"}) {
      java.util.ArrayList<Integer> vals=new java.util.ArrayList<Integer>();
      java.util.Set<Integer> used=side.equals("bottom")?bcu:tcu;
      for(int e:m.component("c").selection(side).entities())if(!used.contains(e))vals.add(e);
      int[] a=new int[vals.size()];for(int k=0;k<a.length;k++)a[k]=vals.get(k);
      m.component("c").selection().create(side+"air","Explicit");
      m.component("c").selection(side+"air").geom("g",2);
      m.component("c").selection(side+"air").set(a);
    }
'''
    code=code[:i]+setup+code[i:]
    i=code.index('    m.component("c").mesh("mesh").create("copy","CopyFace");')
    j=code.index('    m.component("c").mesh("mesh").create("sweep","Sweep");',i)
    copies='''    for(String part:new String[]{"0","1","2","3","4","5","6","air"}) {
      String tag="copy"+part;
      m.component("c").mesh("mesh").create(tag,"CopyFace");
      m.component("c").mesh("mesh").feature(tag).selection("source").named("bottom"+part);
      m.component("c").mesh("mesh").feature(tag).selection("destination").named("top"+part);
    }
'''
    code=code[:i]+copies+code[j:]
    p=ROOT/'data/raw'/name/'PilotHelix.java';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(code)
    return p


if __name__=='__main__':run(build())
