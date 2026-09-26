import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

/** Explicitly resolved, insulated parallel strands. All share one voltage.
 * Coil group is OFF: only TOTAL current is imposed, never equal strand currents.
 */
public class Q2CircularBundle {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_Circular_COMSOL/";
 static Model m;
 static ArrayList<double[]> pts=new ArrayList<double[]>();
 static int[] cu;
 public static Model run() throws Exception {
  new File(OUT).mkdirs();
  new File(OUT+"verification").mkdirs();
  for(int ring=0;ring<=10;ring++){
   int count=ring==0?1:6*ring;
   for(int j=0;j<count;j++){
    double theta=2*Math.PI*j/count;
    pts.add(new double[]{ring*Math.cos(theta),ring*Math.sin(theta),ring,j});
   }
  }
  m=ModelUtil.create("Model");
  m.label("Q2 - 331 straight insulated strands in circular concentric rings - 20 A RMS - 200 kHz");
  m.param().set("N",Integer.toString(pts.size()),"Number of individually insulated copper strands");
  m.param().set("dwire","0.14[mm]","Bare copper strand diameter");
  m.param().set("pitch","0.15[mm]","Radial ring spacing and minimum center-to-center clearance");
  m.param().set("L","1[m]","Evaluation length and out-of-plane thickness");
  m.param().set("sigmaCu","5.8e7[S/m]","Copper conductivity, 20 degC");
  m.param().set("muCu","4*pi*1e-7[H/m]","Copper and surrounding medium permeability");
  m.param().set("Irms","20[A]");
  m.param().set("Ipk","sqrt(2)*Irms");
  m.param().set("f0","200[kHz]");
  m.param().set("Rair","20[mm]","Centered circular exterior magnetic boundary radius");
  m.param().set("Acu","N*pi*dwire^2/4","Total copper area");
  m.param().set("Rdc","L/(sigmaCu*Acu)");
  m.param().set("delta0","sqrt(1/(pi*f0*muCu*sigmaCu))");
  m.param().set("beta","dwire/delta0");
  m.param().set("meshDiv","3","Maximum strand element size: dwire/meshDiv");
  m.component().create("comp1",true);
  m.component("comp1").geom().create("geom1",2);
  m.component("comp1").geom("geom1").create("outer","Circle");
  m.component("comp1").geom("geom1").feature("outer").set("r","Rair");
  for(int i=0;i<pts.size();i++){
   double[] p=pts.get(i);String tag="w"+(i+1);
   m.component("comp1").geom("geom1").create(tag,"Circle");
   m.component("comp1").geom("geom1").feature(tag).set("r","dwire/2");
   m.component("comp1").geom("geom1").feature(tag).set("pos",new String[]{Double.toString(p[0])+"*pitch",Double.toString(p[1])+"*pitch"});
   m.component("comp1").geom("geom1").feature(tag).set("selresult",true);
  }
  m.component("comp1").geom("geom1").run();
  cu=new int[pts.size()];
  for(int i=0;i<cu.length;i++){
   int[] ds=m.component("comp1").selection("geom1_w"+(i+1)+"_dom").entities(2);
   if(ds.length!=1)throw new RuntimeException("Expected one copper domain per strand");
   cu[i]=ds[0];
  }
  m.component("comp1").selection().create("selCu","Explicit");
  m.component("comp1").selection("selCu").geom("geom1",2);
  m.component("comp1").selection("selCu").set(cu);
  m.component("comp1").selection("selCu").label("All 331 copper strands, separated by insulating gaps");
  m.component("comp1").material().create("air","Common");
  m.component("comp1").material("air").label("Nonconducting surrounding medium / insulation gaps");
  m.component("comp1").material("air").selection().all();
  material("air","0[S/m]");
  m.component("comp1").material().create("copper","Common");
  m.component("comp1").material("copper").label("Oxygen-free copper at 20 degC");
  m.component("comp1").material("copper").selection().named("selCu");
  material("copper","sigmaCu");
  m.component("comp1").physics().create("mf","InductionCurrents","geom1");
  m.component("comp1").physics("mf").label("Resolved strands: Az magnetic diffusion");
  m.component("comp1").physics("mf").prop("d").set("d","L");
  m.component("comp1").physics("mf").create("coil","Coil",2);
  m.component("comp1").physics("mf").feature("coil").selection().named("selCu");
  m.component("comp1").physics("mf").feature("coil").set("ConductorModel","Single");
  m.component("comp1").physics("mf").feature("coil").set("coilGroup",false);
  m.component("comp1").physics("mf").feature("coil").set("CoilExcitation","Current");
  m.component("comp1").physics("mf").feature("coil").set("ICoil","Ipk");
  m.component("comp1").physics("mf").feature("coil").label("Parallel terminal connection: TOTAL 20 A RMS, shared voltage");
  m.component("comp1").cpl().create("intCu","Integration");
  m.component("comp1").cpl("intCu").selection().named("selCu");
  m.component("comp1").variable().create("v1");
  m.component("comp1").variable("v1").set("Jrms","mf.normJ/sqrt(2)");
  m.component("comp1").variable("v1").set("qloss","mf.normJ^2/(2*sigmaCu)");
  m.component("comp1").variable("v1").set("Pcu","L*intCu(qloss)");
  m.component("comp1").variable("v1").set("Rac","Pcu/Irms^2");
  m.component("comp1").variable("v1").set("ratio","Rac/Rdc");
  m.component("comp1").variable("v1").set("Icheck","intCu(mf.Jz)/sqrt(2)");
  m.component("comp1").mesh().create("mesh1");
  m.component("comp1").mesh("mesh1").feature("size").set("custom",true);
  m.component("comp1").mesh("mesh1").feature("size").set("hmax","Rair/6");
  m.component("comp1").mesh("mesh1").feature("size").set("hmin","(pitch-dwire)/4");
  m.component("comp1").mesh("mesh1").feature("size").set("hgrad",1.25);
  m.component("comp1").mesh("mesh1").create("tri","FreeTri");
  m.component("comp1").mesh("mesh1").feature("tri").create("szCu","Size");
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").selection().geom("geom1",2);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").selection().set(cu);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("custom",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmaxactive",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmax","dwire/meshDiv");
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hminactive",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmin","(pitch-dwire)/4");
  m.study().create("std1");
  m.study("std1").create("freq","Frequency");
  m.study("std1").feature("freq").set("plist","f0");
  m.study("std1").label("200 kHz: naturally distributed strand currents");
  m.result().numerical().create("global","EvalGlobal");
  m.result().numerical("global").set("expr",new String[]{"freq","Acu","Rdc","comp1.Rac","comp1.ratio","comp1.Pcu","real(comp1.Icheck)","imag(comp1.Icheck)","real(comp1.mf.VCoil_1/Ipk)","real(comp1.mf.VCoil_1)/sqrt(2)","imag(comp1.mf.VCoil_1)/sqrt(2)"});
  m.result().numerical("global").set("unit",new String[]{"Hz","m^2","ohm","ohm","1","W","A","A","ohm","V","V"});
  m.result().table().create("results","Table");
  m.result().table("results").label("RMS current, resistance, copper loss, terminal voltage");
  m.result().numerical("global").set("table","results");
  PrintWriter w=new PrintWriter(OUT+"mesh_convergence.csv","UTF-8");
  w.println("mesh_div,strand_hmax_um,frequency_Hz,copper_area_m2,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A,Rac_voltage_ohm,V_real_rms_V,V_imag_rms_V");
  for(int div:new int[]{3,5,8}){
   m.param().set("meshDiv",Integer.toString(div));
   m.component("comp1").mesh("mesh1").run();
   m.study("std1").run();
   m.result().numerical("global").set("data","dset1");
   double[][] v=m.result().numerical("global").getReal();
   w.printf(Locale.US,"%d,%.12g",div,m.param().evaluate("dwire")*1e6/div);
   for(double[] row:v)w.printf(Locale.US,",%.15g",row[0]);
   w.println();w.flush();
   System.out.println("Q2 MESH "+div+" Rac="+v[3][0]+" ratio="+v[4][0]+" I="+v[6][0]);
   if(Math.abs(v[6][0]-20)>1e-4 || Math.abs(v[3][0]/v[8][0]-1)>1e-3)throw new RuntimeException("Total current or power balance failed");
   m.result().numerical("global").setResult();
   m.save(OUT+"verification/q2_mesh_"+div+".mph");
  }
  w.close();
  globalExport("bundle_200kHz_results.csv");
  strandExport("strand_currents_200kHz.csv");
  plots();
  m.result().table().create("meshTable","Table");
  m.result().table("meshTable").importData(OUT+"mesh_convergence.csv");
  m.result().table("meshTable").label("200 kHz mesh convergence");
  m.result().table().create("strandTable","Table");
  m.result().table("strandTable").importData(OUT+"strand_currents_200kHz.csv");
  m.result().table("strandTable").label("200 kHz: per-strand complex RMS currents and copper loss");
  image("pgJ","comsol_current_density.png");
  image("pgLoss","comsol_loss_density.png");
  m.save(OUT+"Q2_CircularBundle_200kHz.mph");
  m.study("std1").feature("freq").set("plist","10");
  m.study("std1").label("Low-frequency validation: 10 Hz");
  m.study("std1").run();
  globalExport("low_frequency_check.csv");
  strandExport("strand_currents_10Hz.csv");
  m.result().numerical("global").setResult();
  m.result().table().remove("strandTable");
  m.result().table().remove("meshTable");
  m.label("Q2 - low-frequency validation - 10 Hz");
  m.save(OUT+"verification/Q2_LowFrequencyCheck.mph");
  m.study("std1").feature("freq").set("plist","50000 100000 200000 300000 500000 750000 1000000");
  m.study("std1").label("Frequency sweep: 50 kHz to 1 MHz");
  m.study("std1").run();
  globalExport("frequency_sweep.csv");
  m.result().numerical("global").setResult();
  m.label("Q2 - circular bundle frequency sweep - 50 kHz to 1 MHz");
  m.save(OUT+"verification/Q2_FrequencySweep.mph");
  return m;
 }
 static void material(String tag,String sig){
  m.component("comp1").material(tag).propertyGroup("def").set("electricconductivity",new String[]{sig});
  m.component("comp1").material(tag).propertyGroup("def").set("relpermeability",new String[]{"muCu/mu0_const"});
  m.component("comp1").material(tag).propertyGroup("def").set("relpermittivity",new String[]{"1"});
 }
 static void globalExport(String fn)throws Exception{
  double[][] v=m.result().numerical("global").getReal();
  PrintWriter w=new PrintWriter(OUT+fn,"UTF-8");
  w.println("frequency_Hz,copper_area_m2,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A,Rac_voltage_ohm,V_real_rms_V,V_imag_rms_V");
  for(int j=0;j<v[0].length;j++){for(int k=0;k<v.length;k++){if(k>0)w.print(",");w.printf(Locale.US,"%.15g",v[k][j]);}w.println();}w.close();
 }
 static void strandExport(String fn)throws Exception{
  if(!Arrays.asList(m.result().numerical().tags()).contains("strand"))m.result().numerical().create("strand","IntSurface");
  m.result().numerical("strand").set("data","dset1");
  m.result().numerical("strand").set("expr",new String[]{"1","real(mf.Jz)/sqrt(2)","imag(mf.Jz)/sqrt(2)","L*qloss"});
  PrintWriter w=new PrintWriter(OUT+fn,"UTF-8");
  w.println("strand_id,domain_id,radial_ring,ring_index,x_mm,y_mm,radius_mm,area_m2,I_real_rms_A,I_imag_rms_A,loss_W_per_m");
  for(int i=0;i<cu.length;i++){
   double[] p=pts.get(i);double spacing=m.param().evaluate("pitch")*1000;
   double x=p[0]*spacing,y=p[1]*spacing;
   m.result().numerical("strand").selection().set(cu[i]);
   double[][] v=m.result().numerical("strand").getReal();
   w.printf(Locale.US,"%d,%d,%d,%d,%.12g,%.12g,%.12g",i+1,cu[i],(int)p[2],(int)p[3],x,y,Math.hypot(x,y));
   for(double[] row:v)w.printf(Locale.US,",%.15g",row[0]);w.println();
  }w.close();
 }
 static void plots(){
  m.result().dataset().create("copperData","Solution");
  m.result().dataset("copperData").set("solution","sol1");
  m.result().dataset("copperData").selection().geom("geom1",2);
  m.result().dataset("copperData").selection().set(cu);
  m.result().create("pgJ","PlotGroup2D");
  m.result("pgJ").label("331-strand current density RMS (A/mm^2)");
  m.result("pgJ").set("data","copperData");
  m.result("pgJ").create("s1","Surface");
  m.result("pgJ").feature("s1").set("expr","Jrms");
  m.result("pgJ").feature("s1").set("unit","A/mm^2");
  m.result("pgJ").feature("s1").set("colortable","ThermalLight");
  m.result("pgJ").feature("s1").set("resolution","fine");
  m.result().create("pgLoss","PlotGroup2D");
  m.result("pgLoss").label("Copper Joule loss density (W/m^3)");
  m.result("pgLoss").set("data","copperData");
  m.result("pgLoss").create("s1","Surface");
  m.result("pgLoss").feature("s1").set("expr","qloss");
  m.result("pgLoss").feature("s1").set("colortable","ThermalLight");
  m.result("pgJ").set("titletype","manual");
  m.result("pgJ").set("title","Current density RMS (A/mm^2) | 331 straight insulated strands in circular concentric rings");
  m.result("pgLoss").set("titletype","manual");
  m.result("pgLoss").set("title","Time-average copper Joule loss density (W/m^3)");
  m.result("pgJ").run();
  m.result("pgLoss").run();
 }
 static void image(String pg,String file){
  try{
   String tag="img"+pg;m.result().export().create(tag,pg,"Image");
   m.result().export(tag).set("pngfilename",OUT+file);
   m.result().export(tag).set("size","manualweb");
   m.result().export(tag).set("unit","px");
   m.result().export(tag).set("width",1600);
   m.result().export(tag).set("height",1200);
   m.result().export(tag).set("options2d","on");
   m.result().export(tag).set("legend2d","on");
   m.result().export(tag).set("axes2d","on");
   m.result().export(tag).set("title2d","on");
   m.result().export(tag).set("zoomextents",true);
   m.result().export(tag).run();
  }catch(Exception ex){throw new RuntimeException("IMAGE_ERROR "+file,ex);}
 }
 public static void main(String[] args)throws Exception{run();}
}
