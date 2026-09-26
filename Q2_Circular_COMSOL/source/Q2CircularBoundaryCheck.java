import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.Locale;
public class Q2CircularBoundaryCheck {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_Circular_COMSOL/";
 public static Model run()throws Exception{
  Model m=ModelUtil.load("Model",OUT+"Q2_CircularBundle_200kHz.mph");
  m.param().set("Rair","40[mm]");
  m.component("comp1").geom("geom1").run();
  m.component("comp1").mesh("mesh1").run();
  m.study("std1").run();
  double[][] v=m.result().numerical("global").getReal();
  PrintWriter w=new PrintWriter(OUT+"outer_boundary_check.csv","UTF-8");
  w.println("frequency_Hz,copper_area_m2,Rdc_ohm,Rac_ohm,Rac_over_Rdc,loss_W,I_real_rms_A,I_imag_rms_A,Rac_voltage_ohm,V_real_rms_V,V_imag_rms_V");
  for(int k=0;k<v.length;k++){if(k>0)w.print(",");w.printf(Locale.US,"%.15g",v[k][0]);}w.println();w.close();
  m.result().numerical("global").setResult();
  m.label("Q2 validation - outer magnetic boundary radius 40 mm");
  m.result().table().remove("strandTable");
  m.result().table().remove("meshTable");
  m.save(OUT+"verification/Q2_OuterBoundaryCheck.mph");
  return m;
 }
 public static void main(String[] args)throws Exception{run();}
}
