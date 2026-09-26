import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.Locale;

public class Q1Post {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q1_COMSOL/";
 static Model m;
 public static Model run() throws Exception {
  m=ModelUtil.load("Model",OUT+"Q1_SolidWire_200kHz.mph");
  double[][] check=m.result().numerical("gev").getReal();
  if(Math.abs(check[6][0]-20)>1e-5 || Math.abs(check[3][0]-0.0200162)>1e-4) throw new RuntimeException("Validation failed: I="+check[6][0]+" R="+check[3][0]);
  m.result().dataset().create("radline","CutLine2D");
  m.result().dataset("radline").label("Radius at wire midspan z = 0.5 m");
  m.result().dataset("radline").set("data","dset1");
  m.result().dataset("radline").set("genpoints",new double[][]{{0,0.5},{0.001,0.5}});
  m.result().create("pgRad","PlotGroup1D");
  m.result("pgRad").label("Radial current density - RMS");
  m.result("pgRad").set("data","radline");
  m.result("pgRad").create("line1","LineGraph");
  m.result("pgRad").feature("line1").set("expr","Jrms");
  m.result("pgRad").feature("line1").set("unit","A/mm^2");
  m.result("pgRad").feature("line1").set("xdata","expr");
  m.result("pgRad").feature("line1").set("xdataexpr","r");
  m.result("pgRad").feature("line1").set("xdataunit","mm");
  m.result("pgRad").feature("line1").set("resolution","extrafine");
  m.result().dataset().create("revWire","Revolve2D");
  m.result().dataset("revWire").set("data","dset1");
  m.result().dataset("revWire").set("revangle",360);
  m.result().dataset("revWire").set("layermethod","custom");
  m.result().dataset("revWire").set("revlayers",160);
  m.result().dataset("revWire").label("Full cylindrical wire, D=2 mm, L=1 m");
  m.result().dataset().create("cross","CutPlane");
  m.result().dataset("cross").set("data","revWire");
  m.result().dataset("cross").set("quickplane","xy");
  m.result().dataset("cross").set("quickz","L/2");
  m.result().dataset("cross").label("Wire transverse section at z = L/2");
  m.result().create("pgCross","PlotGroup2D");
  m.result("pgCross").label("Cross section - current density RMS (A/mm^2)");
  m.result("pgCross").set("data","cross");
  m.result("pgCross").set("titletype","manual");
  m.result("pgCross").set("title","200 kHz | Current density RMS (A/mm^2) | z = 0.5 m");
  m.result("pgCross").create("surf1","Surface");
  m.result("pgCross").feature("surf1").set("expr","Jrms");
  m.result("pgCross").feature("surf1").set("unit","A/mm^2");
  m.result("pgCross").feature("surf1").set("colortable","ThermalLight");
  m.result("pgCross").feature("surf1").set("resolution","fine");
  m.result().create("pgWire","PlotGroup3D");
  m.result("pgWire").label("Full 1 m copper wire (true aspect ratio)");
  m.result("pgWire").set("data","revWire");
  m.result("pgWire").create("surf1","Surface");
  m.result("pgWire").feature("surf1").set("expr","Jrms");
  m.result("pgWire").feature("surf1").set("unit","A/mm^2");
  m.result("pgWire").feature("surf1").set("colortable","ThermalLight");
  m.result("pgCross").run();
  m.result("pgRad").run();
  m.result().table().create("tblMesh","Table");
  m.result().table("tblMesh").label("Mesh convergence: 8 / 16 / 32 / 64 / 128 radial elements");
  m.result().table("tblMesh").importData(OUT+"mesh_convergence.csv");
  m.result().numerical("gev").setResult();
  m.save(OUT+"Q1_SolidWire_200kHz.mph");
  m.save(OUT+"Q1_SolidWire_generated","java");
  image("pgCross","comsol_cross_section.png");
  image("pgRad","comsol_radial_curve.png");
  // Verify z-invariance with twice as many axial elements.
  m.param().set("Nz","8");
  m.component("comp1").mesh("mesh1").run();
  m.study("std1").run();
  writeResults("axial_mesh_check.csv");
  m.save(OUT+"q1_axial_check_Nz8.mph");
  m.param().set("Nz","4");
  m.component("comp1").mesh("mesh1").run();
  m.study("std1").label("Frequency sweep: 50 kHz to 1 MHz");
  m.study("std1").feature("freq").set("plist","50000 100000 200000 300000 500000 750000 1000000");
  m.result().numerical("gev").setIndex("expr","sqrt(1/(pi*freq*muCu*sigmaCu))",1);
  m.study("std1").run();
  writeResults("frequency_sweep.csv");
  m.result().numerical("gev").setResult();
  m.label("Q1 - Solid copper wire - frequency sweep 50 kHz to 1 MHz");
  m.result("pgCross").set("titletype","auto");
  m.result("pgRad").label("Radial current density - frequency sweep");
  m.save(OUT+"Q1_FrequencySweep.mph");
  return m;
 }
 static void writeResults(String filename) throws Exception {
  double[][] v=m.result().numerical("gev").getReal();
  PrintWriter w=new PrintWriter(OUT+filename,"UTF-8");
  w.println("frequency_Hz,delta_m,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A");
  for(int j=0;j<v[0].length;j++){
   for(int k=0;k<v.length;k++){if(k>0)w.print(",");w.printf(Locale.US,"%.15g",v[k][j]);}w.println();
  }w.close();
 }
 static void image(String pg,String file) {
  try {
   String tag="img"+pg;
   m.result().export().create(tag,pg,"Image");
   m.result().export(tag).set("pngfilename",OUT+file);
   m.result().export(tag).set("size","manualweb");
   m.result().export(tag).set("unit","px");
   m.result().export(tag).set("width",1400);
   m.result().export(tag).set("height",1050);
   m.result().export(tag).set("resolution",160);
   m.result().export(tag).set("background","color");
   m.result().export(tag).set("customcolor",new double[]{1,1,1});
   m.result().export(tag).set("zoomextents",true);
   m.result().export(tag).set("options2d","on");
   m.result().export(tag).set("legend2d","on");
   m.result().export(tag).set("axes2d","on");
   m.result().export(tag).set("title2d","on");
   m.result().export(tag).run();
   System.out.println("IMAGE_SUCCESS "+file);
  }catch(Exception ex){System.out.println("IMAGE_ERROR "+file+": "+ex.getMessage());}
 }
 public static void main(String[] args) throws Exception {run();}
}
