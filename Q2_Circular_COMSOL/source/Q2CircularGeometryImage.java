import com.comsol.model.*;
import com.comsol.model.util.*;
public class Q2CircularGeometryImage {
 public static Model run()throws Exception{
  String OUT="C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q2_Circular_COMSOL/";
  Model m=ModelUtil.load("Model",OUT+"Q2_CircularBundle3D_Geometry.mph");
  m.param().set("L","8[mm]","Short visual excerpt ONLY; saved full geometry has length 1 m");
  m.component("comp3d").geom("geom1").run();
  m.component("comp3d").geom("geom1").image().set("pngfilename",OUT+"comsol_3d_geometry_excerpt.png");
  m.component("comp3d").geom("geom1").image().set("size","manualweb");
  m.component("comp3d").geom("geom1").image().set("unit","px");
  m.component("comp3d").geom("geom1").image().set("width",1400);
  m.component("comp3d").geom("geom1").image().set("height",1100);
  m.component("comp3d").geom("geom1").image().set("options3d","on");
  m.component("comp3d").geom("geom1").image().set("axisorientation","on");
  m.component("comp3d").geom("geom1").image().set("grid","on");
  m.component("comp3d").geom("geom1").image().set("zoomextents",true);
  m.component("comp3d").geom("geom1").image().export();
  return m; // Do not overwrite the full-length geometry file.
 }
 public static void main(String[] args)throws Exception{run();}
}
