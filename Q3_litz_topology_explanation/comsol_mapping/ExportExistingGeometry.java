import com.comsol.model.*;
import com.comsol.model.util.*;
import java.util.Arrays;

/** Read the existing MPH and export its native geometry; never overwrite it. */
public class ExportExistingGeometry {
  static final String OUT="C:/Users/zq257/Documents/Codex/2026-09-26/sa/outputs/litz_topology_explanation/comsol_mapping/";
  static void export(Model m,String name,String view) {
    m.component("paths").geom("g3").image().set("view",view);
    m.component("paths").geom("g3").image().set("imagetype","png");
    m.component("paths").geom("g3").image().set("pngfilename",OUT+name+".png");
    m.component("paths").geom("g3").image().set("size","manualweb");
    m.component("paths").geom("g3").image().set("unit","px");
    m.component("paths").geom("g3").image().set("width",1500);
    m.component("paths").geom("g3").image().set("height",1100);
    m.component("paths").geom("g3").image().set("background","color");
    m.component("paths").geom("g3").image().set("customcolor",new double[]{1,1,1});
    m.component("paths").geom("g3").image().set("options3d","on");
    m.component("paths").geom("g3").image().set("axisorientation","on");
    m.component("paths").geom("g3").image().set("grid","on");
    m.component("paths").geom("g3").image().set("logo3d","on");
    m.component("paths").geom("g3").image().set("fontsize",15);
    m.component("paths").geom("g3").image().set("zoomextents","on");
    m.component("paths").geom("g3").image().export();
    System.out.println("EXPORTED|"+name);
  }
  public static Model run() throws Exception {
    Model m=ModelUtil.load("geometryReadOnly", "C:/Users/zq257/Documents/Codex/2026-09-26/sa/outputs/litz_q3/models/optimized_full_484_paths_fem.mph");
    System.out.println("SOURCE_LOADED|"+m.label());
    System.out.println("GEOMETRY_FEATURES|"+m.component("paths").geom("g3").feature().tags().length);
    String v="mappingview";
    m.component("paths").view().create(v,"g3");
    m.component("paths").view(v).camera().set("projection","orthographic");
    m.component("paths").view(v).camera().set("position",new double[]{.6,-.8,.5});
    m.component("paths").view(v).camera().set("target",new double[]{0,0,.25});
    m.component("paths").view(v).camera().set("up",new double[]{0,0,1});
    m.component("paths").view(v).set("showunits",true);
    export(m,"comsol_484_original_scale",v);
    m.component("paths").view(v).camera().set("viewscaletype","manual");
    m.component("paths").view(v).set("xscale",1);
    m.component("paths").view(v).set("yscale",1);
    m.component("paths").view(v).set("zscale",.02);
    export(m,"comsol_484_axial_compressed_50x",v);
    // Strand 44 is path43 (grid ID 43), corresponding to loop index 462.
    double[][] table=m.component("paths").geom("g3").feature("path43").getDoubleMatrix("table");
    for(int k=0;k<table.length;k++)System.out.println("STRAND44_POINT|"+k+"|"+table[k][0]+"|"+table[k][1]+"|"+table[k][2]);
    for(String f:m.component("paths").geom("g3").feature().tags())
      if(f.startsWith("path")&&!f.equals("path43"))m.component("paths").geom("g3").feature().remove(f);
    m.component("paths").geom("g3").run();
    export(m,"comsol_strand44_axial_compressed_50x",v);
    System.out.println("EXPORT_COMPLETE|original file not saved or modified");
    return m;
  }
  public static void main(String[]args) throws Exception {run();}
}
