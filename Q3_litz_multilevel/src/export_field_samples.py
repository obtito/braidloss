"""Add point evaluation to a freshly built COMSOL Java model, before solving.

Uses solver output on stdout, without changing COMSOL file-security settings.
All sampling coordinates and original field values are retained alongside MPH.
"""
from pathlib import Path
import json
import numpy as np
from recursive_geometry import RecursiveCable
from visual_geometry import tube_surface


def add_field_samples(java, nz=31, ntheta=12):
    java=Path(java);folder=java.parent
    cfg=json.loads((folder/'inputs.json').read_text(encoding='utf-8'))
    c=RecursiveCable(tuple(cfg['grouping_inner_to_outer']),tuple(cfg['signed_relative_pitches_mm']),
        copper_area_mm2=cfg['copper_area_normal_mm2'],gap_um=cfg['gap_design_um'])
    length=cfg['mapping']['cell_length_m']
    z=np.linspace(length*.01,length*.99,nz)
    surface=tube_surface(c,z,ntheta=ntheta,radius_fraction=.98)
    radii=np.array([0.,.2,.4,.6,.8,.98])
    cross=np.stack([tube_surface(c,[length/2],ntheta=24,radius_fraction=r)[:,0] for r in radii],axis=1)
    points=np.concatenate([surface.reshape(-1,3),cross.reshape(-1,3)])
    np.savez_compressed(folder/'sampling_points.npz',surface_m=surface,cross_m=cross,
                        points_m=points,radius_fractions=radii,axial_levels_m=z)
    strings=[]
    # Keep every Java string below the 65535-byte constant-pool limit.
    for part in np.array_split(points,max(1,int(np.ceil(len(points)/400)))):
        strings.append('"'+','.join(format(v,'.16g') for v in part.ravel())+'"')
    code=java.read_text()
    code=code.replace('public class PilotRecursive {','public class PilotRecursive {\n  static String[] fieldCoords={'+','.join(strings)+'};')
    extra='''
    int npoint=0;for(String row:fieldCoords)npoint+=row.split(",").length/3;
    double[][] coords=new double[3][npoint];int offset=0;
    for(String row:fieldCoords) {
      String[] a=row.split(",");
      for(int k=0;k<a.length;k++)coords[k%3][offset+k/3]=Double.parseDouble(a[k]);
      offset+=a.length/3;
    }
    m.result().numerical().create("fields","Interp");
    m.result().numerical("fields").set("coord",coords);
    m.result().numerical("fields").set("expr",new String[]{"Jrms"});
    double[][][] fv=m.result().numerical("fields").getData();
    int ss=fv[0].length-1;
    for(int k=0;k<npoint;k++)System.out.println("FIELD|"+k+"|"+coords[0][k]+"|"+coords[1][k]+"|"+coords[2][k]+"|"+fv[0][ss][k]);
'''
    code=code.replace('    System.out.println("SOLVE_COMPLETE");',extra+'    System.out.println("SOLVE_COMPLETE");')
    java.write_text(code)
    return java
