import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

/** Helicoidal three-component Maxwell solve of 331 solid round copper helices.
 * Exact transverse tube sections, transformed permeability and conductivity.
 * Shared complex axial voltage gradient; total current normalized to 20 A RMS.
 * Zero-conductivity exterior. GMRES handles the harmless magnetic-potential gauge nullspace.
 */
public class Q2RegularTwist {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_RegularTwist_COMSOL/";
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
  m.label("Q2 - 331 regularly twisted insulated strands, 40 mm lay length - 20 A RMS - 200 kHz");
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
  m.param().set("twistPitch","40[mm]");
  m.param().set("Edrive","1[V/m]");
  m.param().set("sigmaAir","0[S/m]","Exactly insulating strands and nonconducting exterior");
  m.param().set("alphaTw","2*pi/twistPitch");
  m.param().set("Rdc","L/(sigmaCu*pi*dwire^2/4*(1+6/sqrt(1+(alphaTw*1*pitch)^2)+12/sqrt(1+(alphaTw*2*pitch)^2)+18/sqrt(1+(alphaTw*3*pitch)^2)+24/sqrt(1+(alphaTw*4*pitch)^2)+30/sqrt(1+(alphaTw*5*pitch)^2)+36/sqrt(1+(alphaTw*6*pitch)^2)+42/sqrt(1+(alphaTw*7*pitch)^2)+48/sqrt(1+(alphaTw*8*pitch)^2)+54/sqrt(1+(alphaTw*9*pitch)^2)+60/sqrt(1+(alphaTw*10*pitch)^2)))");
  m.component().create("comp1",true);
  m.component("comp1").geom().create("geom1",2);
  m.component("comp1").geom("geom1").create("outer","Circle");
  m.component("comp1").geom("geom1").feature("outer").set("r","Rair");
  for(int i=0;i<pts.size();i++){
   double[] p=pts.get(i);String tag="w"+(i+1);
   double rr=p[2]; double th=rr==0?0:2*Math.PI*p[3]/(6*rr);
   String rho="("+rr+"*pitch)", h="(alphaTw*"+rho+")", sl="sqrt(1+"+h+"^2)";
   String theta="("+th+"+alphaTw*dwire/2*"+h+"/"+sl+"*sin(s))";
   String rad="("+rho+"+dwire/2*cos(s))", tang="(dwire/2/"+sl+"*sin(s))";
   String pc="pc"+(i+1);
   m.component("comp1").geom("geom1").create(pc,"ParametricCurve");
   m.component("comp1").geom("geom1").feature(pc).set("parmax","2*pi");
   m.component("comp1").geom("geom1").feature(pc).set("coord",new String[]{rad+"*cos("+theta+")-"+tang+"*sin("+theta+")",rad+"*sin("+theta+")+"+tang+"*cos("+theta+")"});
   m.component("comp1").geom("geom1").create(tag,"ConvertToSolid");
   m.component("comp1").geom("geom1").feature(tag).selection("input").set(pc);
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
  m.component("comp1").selection().create("selAir","Complement");
  m.component("comp1").selection("selAir").set("input",new String[]{"selCu"});
  m.component("comp1").material().create("air","Common");
  m.component("comp1").material("air").label("Perfect electric insulation and nonconducting air");
  m.component("comp1").material("air").selection().all();
  material("air","sigmaAir");
  m.component("comp1").material().create("copper","Common");
  m.component("comp1").material("copper").label("Oxygen-free copper at 20 degC");
  m.component("comp1").material("copper").selection().named("selCu");
  material("copper","sigmaCu");
  m.component("comp1").physics().create("mf","InductionCurrents","geom1");
  m.component("comp1").physics("mf").prop("components").set("components","threecomponent");
  m.component("comp1").physics("mf").create("amp","AmperesLawSolid",2);
  m.component("comp1").physics("mf").feature("amp").selection().all();
  m.component("comp1").physics("mf").label("Helicoidal reduction: three-component vector potential");
  m.component("comp1").physics("mf").prop("d").set("d","L");
  m.component("comp1").physics("mf").create("drive","ExternalCurrentDensity",2);
  m.component("comp1").physics("mf").feature("drive").selection().named("selCu");
  m.component("comp1").physics("mf").feature("drive").set("Je",new String[]{"sigmaCu*alphaTw*y*Edrive","-sigmaCu*alphaTw*x*Edrive","sigmaCu*Edrive"});
  m.component("comp1").cpl().create("intCu","Integration");
  m.component("comp1").cpl("intCu").selection().named("selCu");
  m.component("comp1").cpl().create("intAir","Integration");
  m.component("comp1").cpl("intAir").selection().named("selAir");
  m.component("comp1").variable().create("v1");
  m.component("comp1").variable("v1").set("Jrms","sqrt(abs(mf.Jx-alphaTw*y*mf.Jz)^2+abs(mf.Jy+alphaTw*x*mf.Jz)^2+abs(mf.Jz)^2)/sqrt(2)");
  m.component("comp1").variable("v1").set("qloss","Jrms^2/sigmaCu");
  m.component("comp1").variable("v1").set("Pcu","L*intCu(qloss)");
  m.component("comp1").variable("v1").set("Rac","Pcu/abs(Icheck)^2");
  m.component("comp1").variable("v1").set("ratio","Rac/Rdc");
  m.component("comp1").variable("v1").set("Icheck","intCu(mf.Jz)/sqrt(2)");
  m.component("comp1").variable("v1").set("Pair","0[W]");
  m.component("comp1").variable("v1").set("Pterminal","real(Edrive*L*conj(sqrt(2)*Icheck))/2");
  m.component("comp1").variable("v1").set("powerResidual","(Pterminal-Pcu-Pair)/Pterminal");
  m.component("comp1").mesh().create("mesh1");
  m.component("comp1").mesh("mesh1").feature("size").set("custom",true);
  m.component("comp1").mesh("mesh1").feature("size").set("hmax","Rair/6");
  m.component("comp1").mesh("mesh1").feature("size").set("hmin","2.5[um]");
  m.component("comp1").mesh("mesh1").feature("size").set("hgrad",1.25);
  m.component("comp1").mesh("mesh1").create("tri","FreeTri");
  m.component("comp1").mesh("mesh1").feature("tri").create("szCu","Size");
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").selection().geom("geom1",2);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").selection().set(cu);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("custom",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmaxactive",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmax","dwire/meshDiv");
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hminactive",true);
  m.component("comp1").mesh("mesh1").feature("tri").feature("szCu").set("hmin","2.5[um]");
  m.study().create("std1");
  m.study("std1").create("freq","Frequency");
  m.study("std1").feature("freq").set("plist","f0");
  m.study("std1").label("200 kHz: naturally distributed strand currents");
  m.result().numerical().create("global","EvalGlobal");
  m.result().numerical("global").set("expr",new String[]{"freq","Acu","Rdc","comp1.Rac","comp1.ratio","comp1.Pcu","real(comp1.Icheck)","imag(comp1.Icheck)","real(Edrive*L/(sqrt(2)*comp1.Icheck))","real(Edrive*L)/sqrt(2)","imag(Edrive*L)/sqrt(2)","comp1.Pair","comp1.Pterminal","comp1.powerResidual","comp1.intCu(1)","sigmaAir"});
  m.result().numerical("global").set("unit",new String[]{"Hz","m^2","ohm","ohm","1","W","A","A","ohm","V","V","W","W","1","m^2","S/m"});
  m.result().table().create("results","Table");
  m.result().table("results").label("RMS current, resistance, copper loss, terminal voltage");
  m.result().numerical("global").set("table","results");
  m.param().set("meshDiv","3");
  m.component("comp1").mesh("mesh1").run();
  m.study("std1").createAutoSequences("all");
  configureSolver();
  m.result().numerical("global").set("data","dset1");
  m.save(OUT+"verification/ready_to_solve.mph");
  for(String div:new String[]{"3","5","8"}){
   m.param().set("meshDiv",div);
   m.component("comp1").mesh("mesh1").run();
   solveCurrent("200000","mesh_"+div);
   globalExport("verification/mesh_"+div+".csv");
   m.save(OUT+"verification/mesh_"+div+".mph");
  }
  plots();
  globalExport("regular_twist_200kHz_results.csv");
  strandExport("strand_currents_200kHz.csv");
  image("pgJ","comsol_current_density.png");
  image("pgLoss","comsol_loss_density.png");
  m.save(OUT+"Q2_RegularTwist_200kHz.mph");
  m.param().set("Rair","40[mm]");
  m.component("comp1").geom("geom1").run();
  m.component("comp1").mesh("mesh1").run();
  solveCurrent("200000","outer_boundary_40mm");
  globalExport("verification/outer_boundary_40mm.csv");
  m.save(OUT+"verification/outer_boundary_40mm.mph");
  m.param().set("Rair","20[mm]");
  m.component("comp1").geom("geom1").run();
  m.component("comp1").mesh("mesh1").run();
  solveCurrent("10","low_frequency_10Hz");
  globalExport("verification/low_frequency_10Hz.csv");
  strandExport("verification/strand_currents_10Hz.csv");
  m.save(OUT+"verification/low_frequency_10Hz.mph");
  m.param().set("alphaTw","0[1/m]");
  m.component("comp1").geom("geom1").run();
  m.component("comp1").mesh("mesh1").run();
  solveCurrent("200000","zero_twist_limit");
  globalExport("verification/zero_twist_limit.csv");
  m.save(OUT+"verification/zero_twist_limit.mph");
  System.out.println("REGULAR_TWIST_COMPLETE");
  return m;
 }
 static void configureSolver(){
  m.sol("sol1").feature("v1").feature("comp1_A").set("scalemethod","manual");
  m.sol("sol1").feature("v1").feature("comp1_A").set("scaleval","1e-4");
  m.sol("sol1").feature("s1").set("stol","1e-5");
  m.sol("sol1").feature("s1").create("i1","Iterative");
  m.sol("sol1").feature("s1").feature("i1").set("linsolver","gmres");
  m.sol("sol1").feature("s1").feature("i1").set("maxlinit",200);
  m.sol("sol1").feature("s1").feature("i1").set("prefuntype","right");
  m.sol("sol1").feature("s1").feature("i1").set("rhob",20);
  m.sol("sol1").feature("s1").feature("i1").create("d1","DirectPreconditioner");
  m.sol("sol1").feature("s1").feature("i1").feature("d1").set("linsolver","mumps");
  m.sol("sol1").feature("s1").feature("i1").feature("d1").set("pivotperturb","1e-8");
  m.sol("sol1").feature("s1").feature("fc1").set("linsolver","i1");
 }
 static void solveCurrent(String f,String caseName)throws Exception{
  m.study("std1").feature("freq").set("plist",f);
  m.sol("sol1").feature("s1").feature("p1").set("plistarr",new String[]{f});
  m.param().set("Edrive","1[V/m]");
  m.sol("sol1").runAll();
  double[][] first=m.result().numerical("global").getReal();
  double ir=first[6][0],ii=first[7][0],den=ir*ir+ii*ii;
  double target=m.param().evaluate("Irms");
  double er=target*ir/den,ei=-target*ii/den;
  m.param().set("Edrive","("+er+"+i*("+ei+"))[V/m]");
  m.sol("sol1").runAll();
  double[][] v=m.result().numerical("global").getReal();
  if(Math.abs(v[6][0]-target)>target*1e-5||Math.abs(v[7][0])>target*1e-5)
   throw new RuntimeException("Current normalization failed "+caseName+" "+Arrays.deepToString(v));
  m.result().numerical("global").setResult();
  System.out.println("CASE "+caseName+" "+Arrays.deepToString(v));
 }
 static void material(String tag,String sig){
  String[] T={"1+alphaTw^2*y^2","-alphaTw^2*x*y","alphaTw*y","-alphaTw^2*x*y","1+alphaTw^2*x^2","-alphaTw*x","alphaTw*y","-alphaTw*x","1"};
  String[] sc=new String[9],mu=new String[9];
  for(int k=0;k<9;k++){sc[k]="("+sig+")*("+T[k]+")";mu[k]="muCu/mu0_const*("+T[k]+")";}
  m.component("comp1").material(tag).propertyGroup("def").set("electricconductivity",sc);
  m.component("comp1").material(tag).propertyGroup("def").set("relpermeability",mu);
  m.component("comp1").material(tag).propertyGroup("def").set("relpermittivity",new String[]{"1"});
 }
 static void globalExport(String fn)throws Exception{
  double[][] v=m.result().numerical("global").getReal();
  PrintWriter w=new PrintWriter(OUT+fn,"UTF-8");
  w.println("frequency_Hz,copper_area_m2,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A,Rac_voltage_ohm,V_real_rms_V,V_imag_rms_V,artificial_air_loss_W,terminal_real_power_W,power_balance_relative_residual,copper_transverse_area_m2,air_regularization_S_per_m");
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
  m.result("pgJ").set("title","Current density RMS (A/mm^2) | 331 regularly twisted insulated strands, 40 mm lay length");
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
