from pathlib import Path
out=Path(__file__).resolve().parent
s=(out/'NativeShapeExport.java').read_text(encoding='utf-8')
s=s.replace('public class NativeShapeExport','public class NativeClipboardExport')
old='m.component("paths").geom("g3").image().export();System.out.println("NATIVE_IMAGE_EXPORTED");'
new='''m.result().export().create("img","Image");
  m.result().export("img").set("sourcetype","geometry");
  m.result().export("img").set("sourceobject","g3");
  m.result().export("img").set("view","vmap");
  m.result().export("img").set("size","manualweb");
  m.result().export("img").set("unit","px");
  m.result().export("img").set("width",1500);m.result().export("img").set("height",1100);
  m.result().export("img").set("target","clipboard");
  m.result().export("img").run();
  System.out.println("NATIVE_CLIPBOARD_EXPORTED");'''
assert old in s
s=s.replace(old,new)
(out/'NativeClipboardExport.java').write_text(s,encoding='utf-8')
