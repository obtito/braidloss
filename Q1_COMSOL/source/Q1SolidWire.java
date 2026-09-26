import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.Locale;

/** Q1: native COMSOL AC/DC H-formulation, 2D axisymmetric solid copper wire. */
public class Q1SolidWire {
  static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q1_COMSOL/";
  static Model m;
  public static Model run() throws Exception {

    new File(OUT).mkdirs();
    m=ModelUtil.create("Model");
    m.label("Q1 - Solid copper wire - 200 kHz - 20 A RMS");
    m.param().set("a","1[mm]","Copper radius");
    m.param().set("L","1[m]","Evaluated wire length");
    m.param().set("sigmaCu","5.8e7[S/m]","Copper conductivity at 20 degC");
    m.param().set("muCu","4*pi*1e-7[H/m]","Permeability (mu_r=1)");
    m.param().set("Irms","20[A]","RMS transport current");
    m.param().set("Ipk","sqrt(2)*Irms","Peak phasor current for frequency domain");
    m.param().set("f0","200[kHz]","Operating frequency");
    m.param().set("delta0","sqrt(1/(pi*f0*muCu*sigmaCu))","Theoretical skin depth");
    m.param().set("Rdc","L/(sigmaCu*pi*a^2)","DC resistance");
    m.param().set("Nr","16","Radial mapped elements");
    m.param().set("Nz","4","Axial mapped elements (z invariant solution)");
    m.component().create("comp1",true);
    m.component("comp1").geom().create("geom1",2);
    m.component("comp1").geom("geom1").axisymmetric(true);
    m.component("comp1").geom("geom1").create("r1","Rectangle");
    m.component("comp1").geom("geom1").feature("r1").set("size",new String[]{"a","L"});
    m.component("comp1").geom("geom1").run();
    m.component("comp1").material().create("mat1","Common");
    m.component("comp1").material("mat1").label("Oxygen-free copper, 20 degC");
    m.component("comp1").material("mat1").propertyGroup("def").set("electricconductivity",new String[]{"sigmaCu"});
    m.component("comp1").material("mat1").propertyGroup("def").set("relpermeability",new String[]{"muCu/mu0_const"});
    m.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity",new String[]{"1"});
    m.component("comp1").physics().create("mfh","MagneticFieldFormulation","geom1");
    m.component("comp1").physics("mfh").label("H-phi formulation: axial transport current");
    // COMSOL API 'components' names the current direction: inplane => Hphi only.
    m.component("comp1").physics("mfh").prop("components").set("components","inplane");
    m.component("comp1").physics("mfh").prop("DivergenceConstraint").set("DivergenceConstraint",false);
    m.component("comp1").physics("mfh").create("mfb1","MagneticFieldBoundary",1);
    m.component("comp1").physics("mfh").feature("mfb1").selection().set(4);
    m.component("comp1").physics("mfh").feature("mfb1").label("Surface Hphi = Ipeak / (2 pi a)");
    m.component("comp1").physics("mfh").feature("mfb1").set("H0",new String[]{"0","Ipk/(2*pi*a)","0"});
    m.component("comp1").physics("mfh").create("axisH","PointwiseConstraint",1);
    m.component("comp1").physics("mfh").feature("axisH").selection().set(1);
    m.component("comp1").physics("mfh").feature("axisH").label("Axis regularity: Hphi = 0 at r = 0");
    m.component("comp1").physics("mfh").feature("axisH").set("constraintExpression","Hphi");
    // On this very slender mapped mesh COMSOL's default 0.001*h_spatial
    // axis regularization extends too far radially. Use the exact cylindrical
    // curl, with its regular limit only on the mathematical axis.
    m.component("comp1").physics("mfh").feature("fl1").featureInfo("info").set("mfh.curlHz",new String[]{"Hphir+if(abs(r)<1e-12[m],Hphir,Hphi/r)"});
    m.component("comp1").cpl().create("intV","Integration");
    m.component("comp1").cpl("intV").selection().all();
    m.component("comp1").variable().create("var1");
    m.component("comp1").variable("var1").set("Jrms","mfh.normJ/sqrt(2)","RMS magnitude of current density");
    m.component("comp1").variable("var1").set("qloss","mfh.normJ^2/(2*sigmaCu)","Time-average Joule loss density");
    m.component("comp1").cpl("intV").label("Meridian area integral; use explicit 2*pi*r Jacobian");
    m.component("comp1").variable("var1").set("Pcu","intV(2*pi*r*qloss)","Copper loss over 1 m");
    m.component("comp1").variable("var1").set("Rac","Pcu/Irms^2","AC resistance from Joule losses");
    m.component("comp1").variable("var1").set("ratio","Rac/Rdc","AC/DC resistance ratio");
    m.component("comp1").variable("var1").set("Icheck","intV(2*pi*r*mfh.Jz)/L/sqrt(2)","Complex RMS current, axial average");
    m.component("comp1").mesh().create("mesh1");
    m.component("comp1").mesh("mesh1").create("map1","Map");
    m.component("comp1").mesh("mesh1").feature("map1").create("disr","Distribution");
    m.component("comp1").mesh("mesh1").feature("map1").feature("disr").selection().set(2,3);
    m.component("comp1").mesh("mesh1").feature("map1").feature("disr").set("numelem","Nr");
    m.component("comp1").mesh("mesh1").feature("map1").create("disz","Distribution");
    m.component("comp1").mesh("mesh1").feature("map1").feature("disz").selection().set(1,4);
    m.component("comp1").mesh("mesh1").feature("map1").feature("disz").set("numelem","Nz");
    m.study().create("std1");
    m.study("std1").label("200 kHz frequency-domain solution");
    m.study("std1").create("freq","Frequency");
    m.study("std1").feature("freq").set("plist","f0");
    m.result().numerical().create("gev","EvalGlobal");
    m.result().numerical("gev").set("expr",new String[]{"freq","delta0","Rdc","comp1.Rac","comp1.ratio","comp1.Pcu","real(comp1.Icheck)","imag(comp1.Icheck)"});
    m.result().numerical("gev").set("unit",new String[]{"Hz","m","ohm","ohm","1","W","A","A"});
    m.result().table().create("tbl1","Table");
    m.result().table("tbl1").label("Electrical results (RMS convention)");
    m.result().numerical("gev").set("table","tbl1");
    PrintWriter w=new PrintWriter(OUT+"mesh_convergence.csv","UTF-8");
    try {
      w.println("radial_elements,axial_elements,quad_elements,radial_h_um,frequency_Hz,delta_m,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A");
      for(int nr:new int[]{8,16,32,64,128}){
        m.param().set("Nr",Integer.toString(nr));
        m.component("comp1").mesh("mesh1").run();
        m.study("std1").run();
        m.result().numerical("gev").set("data","dset1");
        double[][] v=m.result().numerical("gev").getReal();
        w.printf(Locale.US,"%d,4,%d,%.12g",nr,nr*4,1000.0/nr);
        for(double[] row:v) w.printf(Locale.US,",%.15g",row[0]);
        w.println(); w.flush();
        System.out.println("MESH nr="+nr+" Rac="+v[3][0]+" ratio="+v[4][0]+" I="+v[6][0]);
        m.result().numerical("gev").setResult();
        m.save(OUT+"q1_mesh_"+nr+".mph");
      }
    } finally {w.close();}
    exportRadial("radial_200kHz.csv");
    m.save(OUT+"Q1_SolidWire_200kHz.mph");
    m.save(OUT+"Q1_SolidWire_generated","java");
    return m;
  }
  static void exportRadial(String name) throws Exception {
    m.result().numerical().create("rad","Interp");
    m.result().numerical("rad").set("data","dset1");
    m.result().numerical("rad").set("expr",new String[]{"r","z","Jrms","real(mfh.Jz)/sqrt(2)","imag(mfh.Jz)/sqrt(2)","qloss","abs(mfh.Jr)/sqrt(2)"});
    double[][] coords=new double[2][2001];
    for(int j=0;j<2001;j++){coords[0][j]=0.001*j/2000;coords[1][j]=0.5;}
    m.result().numerical("rad").setInterpolationCoordinates(coords);
    double[][][] v=m.result().numerical("rad").getData();
    PrintWriter w=new PrintWriter(OUT+name,"UTF-8");
    try {
      w.println("r_m,z_m,J_rms_A_m2,Jz_real_rms_A_m2,Jz_imag_rms_A_m2,loss_W_m3,Jr_rms_A_m2");
      for(int j=0;j<2001;j++){
        for(int k=0;k<v.length;k++){if(k>0)w.print(",");w.printf(Locale.US,"%.15g",v[k][0][j]);}w.println();
      }
    } finally {w.close();}
  }
  public static void main(String[] args) throws Exception {run();}
}

