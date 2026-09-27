"""Prepare a small, separate COMSOL geometry for inspecting the plotted strand.

The original MPH and all parent-task files remain read-only. COMSOL's normal
batch output option saves the returned model; its Java code accesses no files.
"""
from pathlib import Path

out = Path(__file__).resolve().parent
code = (out / 'NativeShapeExport.java').read_text(encoding='utf-8')
code = code.replace('public class NativeShapeExport', 'public class InspectStrand44')
code = code.replace('ModelUtil.create("nativeShape")', 'ModelUtil.create("strand44")')
cut = code.index('  System.out.println("VIEW_PROPERTIES|')
code = code[:cut] + '''  m.component("paths").label("Strand 44: original centerline coordinates, geometry only");
  m.param().set("fullPeriod", "500[mm]", "Full 484-position transposition period");
  m.param().set("copperDiameter", "125.634236[um]", "Reference diameter; no swept solid in this model");
  System.out.println("INSPECTION_GEOMETRY_BUILT|969 original MPH points|no field solve");
  return m;
 }
 public static void main(String[]args) throws Exception {run();}
}
'''
(out / 'InspectStrand44.java').write_text(code, encoding='utf-8')
print('Prepared separate inspection geometry; original model unchanged.')
