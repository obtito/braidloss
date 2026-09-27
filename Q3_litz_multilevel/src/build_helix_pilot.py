"""Small *solid* seven-wire native Helix pilot; not a formal 6 mm² baseline."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]


def build(name="pilot_7_helix_v1", length_mm=0.5, pitch_mm=5., layers=6, mesh_div=4):
    code = (ROOT / "data/raw/pilot_single_verified/PilotSingle.java").read_text()
    code = code.replace("PilotSingle", "PilotHelix")
    code = code.replace("3D solid single wire pilot - voltage driven", "3D solid seven-wire helix pilot - common voltage")
    code = code.replace('m.param().set("L", "1[mm]");', f'm.param().set("L", "{length_mm}[mm]");\n    m.param().set("pitch", "{pitch_mm}[mm]");\n    m.param().set("spacing", "2*a+16[um]");\n    m.param().set("I0", "20[A]*7/343");')
    i = code.index('    g.create("cu", "Cylinder");')
    j = code.index('    g.create("air", "Cylinder");', i)
    geom = '''    g.create("cu0", "Cylinder");
    g.feature("cu0").set("r", "a");
    g.feature("cu0").set("h", "L");
    g.feature("cu0").set("selresult", true);
    for(int j=1;j<7;j++) {
      String t="cu"+j;
      g.create(t,"Helix");
      g.feature(t).set("rmin","a");
      g.feature(t).set("rmaj","spacing");
      g.feature(t).set("axialpitch","pitch");
      g.feature(t).set("turns","L/pitch");
      g.feature(t).set("rot",60*(j-1));
      g.feature(t).set("endcaps","perpaxis");
      g.feature(t).set("selresult",true);
      g.feature(t).set("rtol",1e-6);
    }
'''
    code = code[:i] + geom + code[j:]
    code = code.replace('    System.out.println("STAGE|geometry_built");', '''    m.component("c").selection().create("copper","Union");
    String[] cuNames=new String[7];
    for(int j=0;j<7;j++)cuNames[j]="g_cu"+j+"_dom";
    m.component("c").selection("copper").set("input",cuNames);
    System.out.println("STAGE|geometry_built");''')
    code = code.replace('"g_cu_dom"', '"copper"')
    code = code.replace('"-a-1[nm]"', '"-spacing-a-1[nm]"').replace('"a+1[nm]"', '"spacing+a+1[nm]"')
    i = code.index('    m.component("c").physics("mf").create("coil1","Coil",3);')
    j = code.index('    m.component("c").mesh().create("mesh");', i)
    coils = '''    for(int j=0;j<7;j++) {
      String t="coil"+(j+1);
      m.component("c").physics("mf").create(t,"Coil",3);
      m.component("c").physics("mf").feature(t).selection().named("g_cu"+j+"_dom");
      m.component("c").physics("mf").feature(t).set("ConductorModel","Single");
      m.component("c").physics("mf").feature(t).set("CoilName",Integer.toString(j+1));
      m.component("c").physics("mf").feature(t).set("CoilExcitation","Voltage");
      m.component("c").physics("mf").feature(t).set("VCoil","V0");
      m.component("c").physics("mf").feature(t).feature("ccc1").feature("ct1").selection().named("bottom");
      m.component("c").physics("mf").feature(t).feature("ccc1").create("cg1","CoilGround",2);
      m.component("c").physics("mf").feature(t).feature("ccc1").feature("cg1").selection().named("top");
      m.component("c").cpl().create("intw"+j,"Integration");
      m.component("c").cpl("intw"+j).selection().named("g_cu"+j+"_dom");
    }
    m.component("c").variable().create("vars");
    m.component("c").variable("vars").set("Itot","mf.ICoil_1+mf.ICoil_2+mf.ICoil_3+mf.ICoil_4+mf.ICoil_5+mf.ICoil_6+mf.ICoil_7");
    m.component("c").variable("vars").set("Jrms","mf.normJ*I0/abs(Itot)");
'''
    code = code[:i] + coils + code[j:]
    code = code.replace('"a/6"', f'"a/{mesh_div}"').replace('set("numelem",3)', f'set("numelem",{layers})')
    code = code.replace('"1000 200000 2000000"', '"10 200000"')
    code = code.replace('String[] ex={"freq","mf.ICoil_1","mf.VCoil_1","mf.RCoil_1/L","intcu(mf.Qh)/L","intcu(1)","intcu(mf.Jz)/L"};',
                        'String[] ex={"freq","Itot","V0","real(V0/Itot)/L","2*intcu(mf.Qh)/abs(Itot)^2/L","intcu(1)","2*intcu(mf.Qh)*I0^2/abs(Itot)^2/L"};')
    code = code.replace('    m.result().create("pg","PlotGroup3D");', '''    m.result().numerical().create("strands","EvalGlobal");
    String[] sx=new String[21];
    for(int j=0;j<7;j++) {
      sx[3*j]="mf.ICoil_"+(j+1)+"*I0/Itot";
      sx[3*j+1]="2*intw"+j+"(mf.Qh)*I0^2/abs(Itot)^2/L";
      sx[3*j+2]="intw"+j+"(1)";
    }
    m.result().numerical("strands").set("expr",sx);
    double[][] sr=m.result().numerical("strands").getReal();
    double[][] si=m.result().numerical("strands").getImag();
    for(int j=0;j<sr.length;j++) for(int s=0;s<sr[j].length;s++)
      System.out.println("STRAND|"+j+"|"+s+"|"+sr[j][s]+"|"+si[j][s]);
    m.result().create("pg","PlotGroup3D");''')
    code = code.replace('"mf.normJ"', '"Jrms"')
    job = ROOT / "data/raw" / name
    job.mkdir(parents=True, exist_ok=True)
    path = job / "PilotHelix.java"
    path.write_text(code, encoding="utf-8")
    (job / "inputs.json").write_text(json.dumps(dict(strands=7, formal_baseline=False,
      copper_area_mm2=6*7/343, diameter_from_N=343, frequency_Hz=[10,200000],
      target_RMS_A=20*7/343, pitch_mm=pitch_mm, axial_length_mm=length_mm,
      mesh_axial_layers=layers, mesh_radius_divisor=mesh_div,
      magnetic_boundary="AxA=0 on outer cylinder and physical end faces",
      electrical_termination="Common voltage excitation per single-conductor coil; only physical ends",
      geometry="native solid Helix, plane end cuts normal to global z"), indent=2), encoding="utf-8")
    return path


if __name__ == "__main__":
    from run_java import run
    run(build(*sys.argv[1:]))
