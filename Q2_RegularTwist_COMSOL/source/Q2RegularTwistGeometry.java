import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
public class Q2RegularTwistGeometry {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_RegularTwist_COMSOL/";
 public static Model run()throws Exception{
  new File(OUT).mkdirs();
  Model m=ModelUtil.create("Model");
  m.label("Q2 regular twist: 331 round wires, one 40 mm pitch, repeated 25 times per meter");
  m.param().set("dwire","0.14[mm]","Diameter measured normal to each wire centerline");
  m.param().set("ringStep","0.15[mm]");
  m.param().set("twistPitch","40[mm]");
  m.param().set("L_eval","1[m]","Axial evaluation length, not individual wire length");
  m.param().set("cellTurns","1","One helical period represented by this geometry");
  m.param().set("sigmaCu","5.8e7[S/m]");
  m.component().create("comp3d",true);
  m.component("comp3d").geom().create("geom1",3);
  int id=0;
  for(int ring=0;ring<=10;ring++){
   int count=ring==0?1:6*ring;
   for(int j=0;j<count;j++){
    id++;String tag="strand"+id;
    if(ring==0){
     m.component("comp3d").geom("geom1").create(tag,"Cylinder");
     m.component("comp3d").geom("geom1").feature(tag).set("r","dwire/2");
     m.component("comp3d").geom("geom1").feature(tag).set("h","cellTurns*twistPitch");
    }else{
     m.component("comp3d").geom("geom1").create(tag,"Helix");
     m.component("comp3d").geom("geom1").feature(tag).set("rmaj",ring+"*ringStep");
     m.component("comp3d").geom("geom1").feature(tag).set("rmin","dwire/2");
     m.component("comp3d").geom("geom1").feature(tag).set("axialpitch","twistPitch");
     m.component("comp3d").geom("geom1").feature(tag).set("turns","cellTurns");
     m.component("comp3d").geom("geom1").feature(tag).set("rot",360.0*j/count);
     m.component("comp3d").geom("geom1").feature(tag).set("endcaps","perpaxis");
     m.component("comp3d").geom("geom1").feature(tag).set("rtol",1e-5);
    }
    m.component("comp3d").geom("geom1").feature(tag).set("selresult",true);
   }
  }
  // All tubes are mathematically verified disjoint. Assembly avoids an
  // unnecessary Boolean union of 331 long, mutually insulated helical solids.
  m.component("comp3d").geom("geom1").feature("fin").set("action","assembly");
  m.component("comp3d").geom("geom1").feature("fin").set("createpairs",false);
  m.component("comp3d").geom("geom1").run();
  m.component("comp3d").material().create("cu","Common");
  m.component("comp3d").material("cu").label("331 individually insulated round copper strands");
  m.component("comp3d").material("cu").propertyGroup("def").set("electricconductivity",new String[]{"sigmaCu"});
  m.component("comp3d").material("cu").propertyGroup("def").set("relpermeability",new String[]{"1"});
  m.save(OUT+"Q2_RegularTwist3D_OnePitch.mph");
  m.component("comp3d").geom("geom1").image().set("pngfilename",OUT+"comsol_regular_twist_3d.png");
  m.component("comp3d").geom("geom1").image().set("size","manualweb");
  m.component("comp3d").geom("geom1").image().set("unit","px");
  m.component("comp3d").geom("geom1").image().set("width",1500);
  m.component("comp3d").geom("geom1").image().set("height",1200);
  m.component("comp3d").geom("geom1").image().set("options3d","on");
  m.component("comp3d").geom("geom1").image().set("axisorientation","on");
  m.component("comp3d").geom("geom1").image().set("grid","on");
  m.component("comp3d").geom("geom1").image().set("zoomextents",true);
  m.component("comp3d").geom("geom1").image().export();
  System.out.println("REGULAR_TWIST_GEOMETRY_COMPLETE strands="+id+" pitch=40mm axial_cell_length=40mm");
  return m;
 }
 public static void main(String[] args)throws Exception{run();}
}
