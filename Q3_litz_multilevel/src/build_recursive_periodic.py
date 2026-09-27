"""True swept solids, screw-periodic fields, and strand-identity circuit cycles.

This is a research builder. Every new case needs geometric, power, current
continuity, mesh and cell-length validation before it can be used for design.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from recursive_geometry import RecursiveCable
from run_java import run

ROOT=Path(__file__).resolve().parents[1]
TEMPLATES=Path(__file__).resolve().parent/'templates'


def cycles(permutation):
    todo=set(range(len(permutation)));ans=[]
    while todo:
        k=min(todo);group=[]
        while k in todo:
            todo.remove(k);group.append(k);k=int(permutation[k])
        if k!=group[0]:raise ValueError('Not a permutation cycle')
        ans.append(group)
    return ans


def rows_java(a):
    return 'new double[][]{'+','.join('{'+','.join(format(v,'.16g') for v in row)+'}' for row in a)+'}'


def choose_profile_frames(c,groups):
    """Choose the less rapidly spinning CAD frame, consistently per orbit.

    Both frames are equivariant under axial rotation. Changing the frame of
    a circular isotropic profile changes only its CAD parameterization.
    """
    d=c.evaluate(c.t,1);dd=c.evaluate(c.t,2)
    tangent=d/np.linalg.norm(d,axis=2,keepdims=True)
    guides=[np.array([0.,0.,1.])-tangent[:,:,2,None]*tangent,
            dd-np.sum(dd*tangent,axis=2,keepdims=True)*tangent]
    scores=[]
    for guide in guides:
        norm=np.linalg.norm(guide,axis=2,keepdims=True)
        frame=guide/np.maximum(norm,1e-30)
        derivative=np.gradient(frame,c.t,axis=1,edge_order=2)
        spin=abs(np.sum(derivative*np.cross(tangent,frame),axis=2)).max(axis=1)
        spin[norm[:,:,0].min(axis=1)<1e-12]=np.inf
        scores.append(spin)
    flags=np.zeros(c.n,bool)
    for group in groups:
        flags[group]=max(scores[1][group])<=max(scores[0][group])
    return flags,scores


def build(name='pilot_9_rve_v1', factors=(3,3), pitches=(4.,12.), copper_area=6*9/343,
          length_mm=4/3, mesh_div=3, zscale=.15, air_radius_mm=1., solver='ooc',
          gap_um=20., electrical='cycles', mesh_kind='tet', axial_layers=6,
          profile_frame='axis_projection'):
    c=RecursiveCable(tuple(factors),tuple(pitches),copper_area_mm2=copper_area,gap_um=gap_um)
    n=c.n;length=length_mm*1e-3
    mapping=c.screw_mapping(length,2*np.pi*length/(pitches[-1]*1e-3))
    if not mapping['geometry_mapping_pass']:raise ValueError('Screw periodic geometry did not match')
    perm=mapping['strand_destination_to_source'];groups=cycles(perm)
    if electrical not in ['cycles','identity_voltage','uniform_voltage_diagnostic']:raise ValueError(electrical)
    use_ge=electrical=='cycles'
    if electrical=='identity_voltage' and perm!=list(range(n)):
        raise ValueError('Common voltage per cell is valid here only with identity strand mapping')
    gids=np.zeros(n,int)
    for i,group in enumerate(groups):gids[group]=i
    if profile_frame in ['curvature','adaptive']:
        intervals=max(128,int(np.ceil(128*length/(min(abs(x) for x in pitches)*1e-3))))
        step=length/intervals
        pad=int(np.ceil(2*c.a/step))+2
        z=np.arange(-pad,intervals+pad+1)*step
    else:z=np.linspace(-2*c.a,length+2*c.a,97)
    pts,_,tt=c.at_physical_z(z)
    frames=[];straight=[]
    if profile_frame not in ['axis_projection','curvature','adaptive']:raise ValueError(profile_frame)
    profile_curvature=np.full(n,profile_frame=='curvature')
    if profile_frame=='adaptive':profile_curvature,_=choose_profile_frames(c,groups)
    for j in range(n):
        cs=CubicSpline(c.t,c.xyz[j]-c.t[:,None]*[0,0,1],bc_type='periodic',extrapolate='periodic')
        tang=cs(tt[j,0],1)+[0,0,1];tang/=np.linalg.norm(tang)
        guide=np.array([0.,0.,1.])-tang[2]*tang
        if profile_curvature[j]:
            second=cs(tt[j,0],2)
            guide=second-np.dot(second,tang)*tang
        isstraight=np.linalg.norm(np.ptp(pts[j,:,:2],axis=0))<1e-12
        if not isstraight and np.linalg.norm(guide)<1e-9:
            raise ValueError('Projection frame is singular at the sweep start; change the cut offset')
        straight.append(isstraight)
        if isstraight:e1=np.array([1.,0.,0.])
        else:e1=guide/np.linalg.norm(guide)
        e2=np.cross(tang,e1)
        frames.append(np.array([pts[j,0],pts[j,0]+c.a*e1,pts[j,0]+c.a*e2]))
    code=(TEMPLATES/'recursive_base.java.txt').read_text()
    i=code.index('  static String[] paths=');j=code.index('  public static Model run()',i)
    pathtext=','.join('"'+','.join(format(v,'.16g') for v in p.ravel())+'"' for p in pts)
    extra='  static String[] paths={'+pathtext+'};\n'
    extra+='  static double[][][] frames=new double[][][]{'+','.join(rows_java(f) for f in frames)+'};\n'
    extra+='  static int[] perm=new int[]{'+','.join(map(str,perm))+'};\n'
    extra+='  static int[] groupID=new int[]{'+','.join(map(str,gids))+'};\n'
    extra+='  static boolean[] straight=new boolean[]{'+','.join('true' if b else 'false' for b in straight)+'};\n'
    extra+='  static boolean[] profileCurvature=new boolean[]{'+','.join('true' if b else 'false' for b in profile_curvature)+'};\n'
    helper=(TEMPLATES/'twisted_cylinder_method.java.txt').read_text()
    extra+=helper
    code=code[:i]+extra+code[j:]
    code=code.replace('"sqrt(6[mm^2]/343/pi)"',f'"{c.a:.16g}[m]"')
    # Replace named parameter lines without modifying coordinates embedded above.
    import re
    values={'L':f'{length_mm:.16g}[mm]','pitch':f'{pitches[-1]:.16g}[mm]',
            'I0':f'{20*copper_area/6:.16g}[A]','Rair':f'{air_radius_mm}[mm]',
            'spacing':f'{c.envelope_radius:.16g}[m]'}
    for key,val in values.items():code=re.sub(r'm\.param\(\)\.set\("'+key+r'", "[^"]*"\);',f'm.param().set("{key}", "{val}");',code)
    code=code.replace('new String[9]',f'new String[{n}]').replace('j<9',f'j<{n}').replace('new String[27]',f'new String[{3*n}]')
    code=code.replace('      g.create(wp,"WorkPlane");','''      if(straight[j]) { twistedCylinder(g,cu,"a");continue; }
      g.create(wp,"WorkPlane");
      g.feature(wp).set("planetype","coordinates");
      g.feature(wp).set("genpoints",frames[j]);''')
    code=code.replace('g.feature(sw).set("movetospine",true);','g.feature(sw).set("movetospine",false);\n      g.feature(sw).set("twisting","projvector");\n      g.feature(sw).set("projvector",new String[]{"0","0","1"});')
    if profile_frame in ['curvature','adaptive']:
        code=code.replace('g.feature(sw).set("twisting","projvector");','g.feature(sw).set("twisting",profileCurvature[j]?"curvature":"projvector");')
    code=code.replace('g.feature(cu).set("selresult",true);','g.feature(cu).set("selresult",true);\n      g.feature(cu).set("selresultshow","all");')
    i=code.index('    g.create("air", "Cylinder");');j=code.index('    g.run();',i)
    code=code[:i]+'    twistedCylinder(g,"air","Rair");\n'+code[j:]
    i=code.index('    m.component("c").physics("mf").create("gfa1"')
    pc='''    m.component("c").selection().create("caps","Union");
    m.component("c").selection("caps").set("entitydim",2);
    m.component("c").selection("caps").set("input",new String[]{"bottom","top"});
    m.component("c").coordSystem().create("rot","Rotated");
    m.component("c").coordSystem("rot").set("angle",new String[]{"360[deg]*L/pitch","0","0"});
    m.component("c").physics("mf").create("pc","PeriodicCondition",2);
    m.component("c").physics("mf").feature("pc").selection().named("caps");
    m.component("c").physics("mf").feature("pc").set("TransformationMethod","GlobalSystem");
    m.component("c").physics("mf").feature("pc").set("manualDestinationSelection",true);
    m.component("c").physics("mf").feature("pc").selection("destinationDomains").named("top");
    m.component("c").physics("mf").feature("pc").set("TransformationMethod_dst","rot");
'''
    code=code[:i]+pc+code[i:]
    if use_ge:
        code=code.replace('set("CoilExcitation","Voltage")','set("CoilExcitation","Current")')
        code=code.replace('set("VCoil","V0")','set("ICoil","Is"+groupID[j])')
    for tag in ['ct1','cg1']:
        key=f'      m.component("c").physics("mf").feature(t).feature("ccc1").feature("{tag}").selection()'
        code=code.replace(key,f'      m.component("c").physics("mf").feature(t).feature("ccc1").feature("{tag}").set("SlantedCut",true);\n'+key)
    code=re.sub(r'set\("Itot","[^"]*"\)', 'set("Itot","'+'+'.join('mf.ICoil_'+str(k+1) for k in range(n))+'")',code)
    i=code.index('    m.component("c").mesh().create("mesh");')
    setup=(TEMPLATES/'periodic_cap_selections.java.txt').read_text().replace('j<7',f'j<{n}')
    ge='''    m.component("c").physics().create("ge","GlobalEquations","g");
    m.component("c").physics("ge").prop("EquationForm").set("form","Automatic");
'''
    names=['Is'+str(k) for k in range(len(groups))]
    eqs=['('+'+'.join('mf.VCoil_'+str(j+1) for j in group)+')/'+str(len(group))+'-V0' for group in groups]
    for prop,arr in [('name',names),('equation',eqs),('initialValueU',['0.05']*len(groups))]:
        ge+='    m.component("c").physics("ge").feature("ge1").set("'+prop+'",new String[]{'+','.join('"'+x+'"' for x in arr)+'});\n'
    for key,val in [('DependentVariableQuantity','none'),('SourceTermQuantity','none'),('CustomDependentVariableUnit','A'),('CustomSourceTermUnit','V')]:
        ge+=f'    m.component("c").physics("ge").feature("ge1").set("{key}","{val}");\n'
    code=code[:i]+setup+(ge if use_ge else '')+code[i:]
    i=code.index('    m.component("c").mesh("mesh").create("sweep","Sweep");')
    j=code.index('    m.component("c").mesh("mesh").run();',i)
    mesh='''    for(int j=0;j<perm.length;j++) {
      String tag="copy"+j;
      m.component("c").mesh("mesh").create(tag,"CopyFace");
      m.component("c").mesh("mesh").feature(tag).selection("source").named("bottom"+perm[j]);
      m.component("c").mesh("mesh").feature(tag).selection("destination").named("top"+j);
    }
    m.component("c").mesh("mesh").create("copyair","CopyFace");
    m.component("c").mesh("mesh").feature("copyair").selection("source").named("bottomair");
    m.component("c").mesh("mesh").feature("copyair").selection("destination").named("topair");
    m.component("c").mesh("mesh").create("tet","FreeTet");
    m.component("c").mesh("mesh").feature("tet").set("zscale",ZSCALE);
    m.component("c").mesh("mesh").feature("tet").create("cuSize","Size");
    m.component("c").mesh("mesh").feature("tet").feature("cuSize").selection().named("copper");
    m.component("c").mesh("mesh").feature("tet").feature("cuSize").set("custom",true);
    m.component("c").mesh("mesh").feature("tet").feature("cuSize").set("hmaxactive",true);
    m.component("c").mesh("mesh").feature("tet").feature("cuSize").set("hmax","a/DIV");
'''.replace('ZSCALE',str(zscale)).replace('DIV',str(mesh_div))
    if mesh_kind=='sweep':
        if perm!=list(range(n)):raise ValueError('Swept mesh needs identity face mapping in this builder')
        mesh=mesh[:mesh.index('    m.component("c").mesh("mesh").create("tet"')]+'''
    m.component("c").mesh("mesh").create("sweepAll","Sweep");
    m.component("c").mesh("mesh").feature("sweepAll").selection("sourceface").named("bottom");
    m.component("c").mesh("mesh").feature("sweepAll").selection("targetface").named("top");
    m.component("c").mesh("mesh").feature("sweepAll").create("dist","Distribution");
    m.component("c").mesh("mesh").feature("sweepAll").feature("dist").set("numelem",LAYERS);
'''.replace('LAYERS',str(axial_layers))
    elif mesh_kind!='tet':raise ValueError(mesh_kind)
    code=code[:i]+mesh+code[j:]
    code=code.replace('"a/4"',f'"a/{mesh_div}"')
    if use_ge:
        code=code.replace('m.study("std").create("ccc","CoilCurrentCalculation");','m.study("std").create("ccc","CoilCurrentCalculation");\n    m.study("std").feature("ccc").activate("ge",false);')
    # Circuit currents and field-derived coil voltages must be solved together.
    # The automatic segregated GE block has no self derivative and is singular.
    sol='    m.study("std").createAutoSequences("all");\n'
    if use_ge:sol+='''
    m.sol("sol1").feature("s2").create("fc1","FullyCoupled");
    m.sol("sol1").feature("s2").feature().remove("se1");
    m.sol("sol1").feature("s2").feature("fc1").set("linsolver","d1");
'''
    # Explicit physical scales avoid an automatic 1e2 T*m A-field estimate
    # caused by near-zero diagonal entries of a periodic curl formulation.
    sol+='''    m.sol("sol1").feature("v2").feature("c_A").set("scalemethod","manual");
    m.sol("sol1").feature("v2").feature("c_A").set("scaleval",1e-7);
'''
    if use_ge:sol+='''
    m.sol("sol1").feature("v2").feature("c_ODE1").set("scalemethod","manual");
    m.sol("sol1").feature("v2").feature("c_ODE1").set("scaleval",0.05);
'''
    for j in range(n) if use_ge else []:
        sol+=f'    m.sol("sol1").feature("v2").feature("c_mf_coil{j+1}_VCoil_ode").set("scalemethod","manual");\n'
        sol+=f'    m.sol("sol1").feature("v2").feature("c_mf_coil{j+1}_VCoil_ode").set("scaleval",1e-4);\n'
    if solver in ['ooc','mumps']:
        sol+='''    m.sol("sol1").feature("s2").feature("d1").set("ooc","on");
    m.sol("sol1").feature("s2").feature("d1").set("incore","manual");
    m.sol("sol1").feature("s2").feature("d1").set("oocmemory",512);
'''
        if solver=='mumps':
            sol+='    m.sol("sol1").feature("s2").feature("d1").set("linsolver","mumps");\n'
    elif solver=='iterative':
        sol+='''    m.sol("sol1").feature("s2").feature("fc1").set("linsolver","i1");
    m.sol("sol1").feature("s2").feature("i1").set("maxlinit",800);
    m.sol("sol1").feature("s2").set("stol",1e-6);
'''
    elif solver!='direct':raise ValueError(solver)
    sol+='    m.sol("sol1").runAll();'
    code=code.replace('    m.study("std").run();',sol)
    p=ROOT/'data/raw'/name/'PilotRecursive.java';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(code)
    info=c.diagnostics();info.update(mapping=mapping,current_cycles=groups,mesh_radius_divisor=mesh_div,zscale=zscale,
        air_radius_mm=air_radius_mm,solver=solver,electrical=electrical,mesh_kind=mesh_kind,profile_frame=profile_frame,
        axial_layers=axial_layers if mesh_kind=='sweep' else None,profile_curvature=profile_curvature.tolist(),
        formal_baseline=False,validation_status='pending')
    (p.parent/'inputs.json').write_text(json.dumps(info,indent=2))
    return p


if __name__=='__main__':run(build())
