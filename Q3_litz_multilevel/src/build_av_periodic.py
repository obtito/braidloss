"""Full periodic A-phi field formulation with no strand-current prescription.

J = sigma * (-grad(phi) - i*omega*A + (V0/L)*ez).
Magnetic curl-curl and electric current conservation are solved together.
Both A and phi have the same verified screw map. A potential reference is set
at one point per connected permutation cycle; it is not a physical busbar.
All cases require field, power, mesh and cell-length verification.
"""
from pathlib import Path
import json
import numpy as np
from build_recursive_periodic import build, cycles, rows_java
from recursive_geometry import RecursiveCable
from run_java import run


def build_av(name='pilot_9_av_v1', av_ooc=True, **kwargs):
    kwargs.update(electrical='uniform_voltage_diagnostic',solver='ooc')
    p=build(name=name,**kwargs)
    cfg=json.loads((p.parent/'inputs.json').read_text(encoding='utf-8'))
    c=RecursiveCable(tuple(cfg['grouping_inner_to_outer']),tuple(cfg['signed_relative_pitches_mm']),
       copper_area_mm2=cfg['copper_area_normal_mm2'],gap_um=cfg['gap_design_um'])
    n=c.n;perm=cfg['mapping']['strand_destination_to_source'];groups=cycles(perm)
    centers=c.at_physical_z([0.])[0][:,0]
    code=p.read_text()
    # MF's default frequency-domain gauge feature imposes current conservation,
    # which would duplicate EC's equation. The stationary form imposes div(A)=0
    # and is the required Coulomb gauge for this separately coupled A-phi system.
    code=code.replace('m.component("c").physics("mf").feature("gfa1").selection().all();',
       'm.component("c").physics("mf").feature("gfa1").selection().all();\n'
       '    m.component("c").physics("mf").feature("gfa1").set("equationForm","Stationary");')
    i=code.index('    for(int j=0;j<'+str(n)+';j++) {\n      String t="coil"')
    j=code.index('    m.component("c").variable().create("vars");',i)
    setup='''
    m.component("c").physics("mf").create("als","AmperesLawSolid",3);
    m.component("c").physics("mf").feature("als").selection().named("copper");
    m.component("c").physics("mf").create("ext","ExternalCurrentDensity",3);
    m.component("c").physics("mf").feature("ext").selection().named("copper");
    m.component("c").physics("mf").feature("ext").set("Je",new String[]{"-sig*Vx","-sig*Vy","sig*(V0/L-Vz)"});
    m.component("c").physics().create("ec","ConductiveMedia","g");
    m.component("c").physics("ec").selection().named("copper");
    m.component("c").physics("ec").create("ext","ExternalCurrentDensity",3);
    m.component("c").physics("ec").feature("ext").selection().named("copper");
    m.component("c").physics("ec").feature("ext").set("Je",new String[]{"-sig*i*2*pi*freq*Ax","-sig*i*2*pi*freq*Ay","sig*(V0/L-i*2*pi*freq*Az)"});
    m.component("c").physics("ec").create("pc","PeriodicCondition",2);
    m.component("c").physics("ec").feature("pc").selection().named("caps");
    m.component("c").physics("ec").feature("pc").set("TransformationMethod","GlobalSystem");
    m.component("c").physics("ec").feature("pc").set("manualDestinationSelection",true);
    m.component("c").physics("ec").feature("pc").selection("destinationDomains").named("top");
    m.component("c").physics("ec").feature("pc").set("TransformationMethod_dst","rot");
    double[][] center0=CENTERS;
    int[] refs=REFS;
    double[][] vertex=g.getVertexCoord();
    int[] pointIDs=new int[refs.length];
    for(int j0=0;j0<refs.length;j0++) {
      int j=refs[j0];double best=1e99;int point=-1;
      int[] candidates=m.component("c").selection("g_cu"+j+"_pnt").entities(0);
      for(int pointCandidate:candidates)if(Math.abs(vertex[2][pointCandidate-1])<1e-9) {
        int k=pointCandidate-1;
        double dx=vertex[0][k]-center0[j][0],dy=vertex[1][k]-center0[j][1];
        double dd=dx*dx+dy*dy;if(dd<best){best=dd;point=k+1;}
      }
      if(point<1)throw new RuntimeException("No potential-reference vertex found");
      pointIDs[j0]=point;
      System.out.println("REFERENCE_POINT|"+j+"|"+point+"|"+Math.sqrt(best));
    }
    m.component("c").physics("ec").create("ref","Ground",0);
    m.component("c").physics("ec").feature("ref").selection().set(pointIDs);
    for(int j=0;j<NWIRE;j++) {
      m.component("c").cpl().create("intw"+j,"Integration");
      m.component("c").cpl("intw"+j).selection().named("g_cu"+j+"_dom");
    }
    System.out.println("STAGE|coupled_physics_built");
'''.replace('CENTERS',rows_java(centers)).replace('REFS','new int[]{'+','.join(str(g[0]) for g in groups)+'}').replace('NWIRE',str(n))
    code=code[:i]+setup+code[j:]
    i=code.index('    m.component("c").variable("vars").set("Itot"')
    j=code.index('    java.util.Set<Integer> bcu',i)
    setup='''
    m.component("c").variable("vars").set("Jxphys","-sig*(Vx+i*2*pi*freq*Ax)");
    m.component("c").variable("vars").set("Jyphys","-sig*(Vy+i*2*pi*freq*Ay)");
    m.component("c").variable("vars").set("Jzphys","sig*(V0/L-Vz-i*2*pi*freq*Az)");
    m.component("c").variable("vars").set("Jnormphys","sqrt(abs(Jxphys)^2+abs(Jyphys)^2+abs(Jzphys)^2)");
    m.component("c").variable("vars").set("Qphys","Jnormphys^2/(2*sig)");
    m.component("c").variable("vars").set("Itot","intcu(Jzphys)/L");
    m.component("c").variable("vars").set("Jrms","Jnormphys*I0/abs(Itot)");
'''
    code=code[:i]+setup+code[j:]
    i=code.index('    m.study().create("std");');j=code.index('    m.result().numerical().create("ev"',i)
    setup='''
    for(int j=0;j<NWIRE;j++)for(String side:new String[]{"bottom","top"}) {
      String t="int"+side+j;
      m.component("c").cpl().create(t,"Integration");
      m.component("c").cpl(t).selection().geom("g",2);
      m.component("c").cpl(t).selection().named(side+j);
    }
    m.study().create("std");
    m.study("std").create("freq","Frequency");
    m.study("std").feature("freq").set("plist","10 200000");
    m.study("std").createAutoSequences("all");
    m.sol("sol1").feature("s1").create("allCoupled","FullyCoupled");
    m.sol("sol1").feature("s1").create("dirAV","Direct");
    m.sol("sol1").feature("s1").feature("dirAV").set("linsolver","pardiso");
    m.sol("sol1").feature("s1").feature("dirAV").set("ooc","on");
    m.sol("sol1").feature("s1").feature("dirAV").set("incore","manual");
    m.sol("sol1").feature("s1").feature("dirAV").set("oocmemory",512);
    m.sol("sol1").feature("s1").feature("allCoupled").set("linsolver","dirAV");
    m.sol("sol1").feature("v1").feature("c_A").set("scalemethod","manual");
    m.sol("sol1").feature("v1").feature("c_A").set("scaleval",1e-7);
    m.sol("sol1").feature("v1").feature("c_V").set("scalemethod","manual");
    m.sol("sol1").feature("v1").feature("c_V").set("scaleval",1e-4);
    m.sol("sol1").runAll();
'''.replace('NWIRE',str(n))
    code=code[:i]+setup+code[j:]
    code=code.replace('mf.Qh','Qphys')
    code=code.replace('sx[3*j]="mf.ICoil_"+(j+1)+"*I0/Itot";',
                      'sx[3*j]="intw"+j+"(Jzphys)/L*I0/Itot";')
    i=code.index('    m.result().create("pg"')
    setup='''
    m.result().numerical().create("cuts","EvalGlobal");
    String[] cx=new String[2*NWIRE];
    for(int j=0;j<NWIRE;j++) {
      cx[2*j]="intbottom"+j+"(Jzphys)*I0/Itot";
      cx[2*j+1]="inttop"+j+"(Jzphys)*I0/Itot";
    }
    m.result().numerical("cuts").set("expr",cx);
    double[][] cr=m.result().numerical("cuts").getReal(),ci=m.result().numerical("cuts").getImag();
    for(int j=0;j<cr.length;j++)for(int s=0;s<cr[j].length;s++)
      System.out.println("CUT|"+j+"|"+s+"|"+cr[j][s]+"|"+ci[j][s]);
'''.replace('NWIRE',str(n))
    code=code[:i]+setup+code[i:]
    code=code.replace('3D round swept solids - recursive 3 x 3 stranding pilot','Full 3D periodic A-phi electromagnetic field pilot')
    if not av_ooc:
        code=code.replace('feature("dirAV").set("ooc","on")','feature("dirAV").set("ooc","off")')
    p.write_text(code)
    cfg.update(electrical='full_A_phi_periodic',av_ooc=av_ooc,formal_baseline=False,
               current_prescription='None; common axial voltage gradient; periodic potential; one gauge reference per permutation cycle')
    (p.parent/'inputs.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    return p


if __name__=='__main__':run(build_av())
