from pathlib import Path
import zipfile, re, json, xml.etree.ElementTree as ET
import numpy as np

out=Path(__file__).resolve().parent
source=out.parent.parent/'litz_q3/models/optimized_full_484_paths_fem.mph'
s=zipfile.ZipFile(source).read('dmodel.xml').decode('utf-8')
selected={}
for gid in [43,21]:
    i=s.index('tag="path'+str(gid)+'"')
    a=s.rfind('<GeomFeature ',0,i);b=s.index('</GeomFeature>',i)+len('</GeomFeature>')
    block=s[a:b];root=ET.fromstring(block)
    prop=next(e for e in root.findall('propertyValue') if e.get('name')=='p:table')
    pts=np.array([float(x) for x in re.findall(r"'([^']*)'",prop.get('valueMatrix'))]).reshape(-1,3)
    selected[str(gid)]={'label':root.get('name'),'points_m':pts.tolist()}
    print('EXTRACTED',root.get('name'),pts.shape,'first/last',pts[[0,-1]].tolist())
(out/'actual_mph_strands.json').write_text(json.dumps(selected,separators=(',',':')),encoding='utf-8')
table=';'.join(','.join(format(x,'.12g') for x in row) for row in selected['43']['points_m'])
assert len(table.encode())<65000
code='''import com.comsol.model.*;
import com.comsol.model.util.*;
public class NativeShapeExport {
 public static Model run() throws Exception {
  Model m=ModelUtil.create("nativeShape");m.label("Strand 44 coordinates extracted from existing MPH; geometry illustration");
  m.component().create("paths",true);m.component("paths").geom().create("g3",3);
  String data="__TABLE__";String[] rows=data.split(";");double[][] xyz=new double[rows.length][3];
  for(int k=0;k<rows.length;k++){String[] s=rows[k].split(",");for(int j=0;j<3;j++)xyz[k][j]=Double.parseDouble(s[j]);}
  m.component("paths").geom("g3").create("path43","InterpolationCurve");
  m.component("paths").geom("g3").feature("path43").label("Strand 44 - original coordinate table");
  m.component("paths").geom("g3").feature("path43").set("table",xyz);
  m.component("paths").geom("g3").run();
  m.component("paths").view().create("vmap","g3");
  m.component("paths").view("vmap").camera().set("projection","orthographic");
  m.component("paths").view("vmap").camera().set("viewscaletype","automatic");
  m.component("paths").view("vmap").camera().set("autocontext","anisotropic");
  m.component("paths").view("vmap").camera().set("xweight",1);m.component("paths").view("vmap").camera().set("yweight",1);m.component("paths").view("vmap").camera().set("zweight",2);
  System.out.println("VIEW_PROPERTIES|"+java.util.Arrays.toString(m.component("paths").view("vmap").properties()));
  System.out.println("CAMERA_PROPERTIES|"+java.util.Arrays.toString(m.component("paths").view("vmap").camera().properties()));
  m.component("paths").geom("g3").image().set("view","vmap");
  m.component("paths").geom("g3").image().set("pngfilename","__OUT__/comsol_native_strand44.png");
  m.component("paths").geom("g3").image().set("size","manualweb");
  m.component("paths").geom("g3").image().set("width",1500);m.component("paths").geom("g3").image().set("height",1100);
  m.component("paths").geom("g3").image().set("options3d","on");m.component("paths").geom("g3").image().set("grid","on");
  m.component("paths").geom("g3").image().set("axisorientation","on");m.component("paths").geom("g3").image().set("logo3d","on");
  m.component("paths").geom("g3").image().set("zoomextents","on");
  m.component("paths").geom("g3").image().export();System.out.println("NATIVE_IMAGE_EXPORTED");
  return m;
 }
 public static void main(String[]args) throws Exception {run();}
}
'''.replace('__TABLE__',table).replace('__OUT__',out.as_posix())
(out/'NativeShapeExport.java').write_text(code,encoding='utf-8')
