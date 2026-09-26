import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
public class Q2CircularGeometry {
 static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_Circular_COMSOL/";
 public static Model run() throws Exception {
  new File(OUT).mkdirs();
  Model m=ModelUtil.create("Model");
  m.label("Q2 - 3D circular bundle - 331 insulated straight strands - 1 m");
  m.param().set("dwire","0.14[mm]");
  m.param().set("pitch","0.15[mm]");
  m.param().set("L","1[m]");
  m.param().set("sigmaCu","5.8e7[S/m]");
  m.component().create("comp3d",true);
  m.component("comp3d").geom().create("geom1",3);
  int id=0;
  for(int ring=0;ring<=10;ring++){
   int count=ring==0?1:6*ring;
   for(int j=0;j<count;j++){
    double theta=2*Math.PI*j/count;
    id++;String tag="strand"+id;
    m.component("comp3d").geom("geom1").create(tag,"Cylinder");
    m.component("comp3d").geom("geom1").feature(tag).set("r","dwire/2");
    m.component("comp3d").geom("geom1").feature(tag).set("h","L");
    m.component("comp3d").geom("geom1").feature(tag).set("pos",new String[]{Double.toString(ring*Math.cos(theta))+"*pitch",Double.toString(ring*Math.sin(theta))+"*pitch","0"});
   }
  }
  m.component("comp3d").geom("geom1").run();
  m.component("comp3d").material().create("cu","Common");
  m.component("comp3d").material("cu").label("Copper; strands separated by 10 micrometer insulation gaps");
  m.component("comp3d").material("cu").propertyGroup("def").set("electricconductivity",new String[]{"sigmaCu"});
  m.component("comp3d").material("cu").propertyGroup("def").set("relpermeability",new String[]{"1"});
  m.save(OUT+"Q2_CircularBundle3D_Geometry.mph");
  System.out.println("3D_GEOMETRY_COMPLETE "+id+" isolated copper cylinders, 1 m each");
  return m;
 }
 public static void main(String[] args)throws Exception{run();}
}
