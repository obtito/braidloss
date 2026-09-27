"""Export a local, closed copper-surface mesh for visualization, not FEM/CAD."""
from pathlib import Path
import json, struct
import numpy as np

out=Path(__file__).resolve().parent
D=json.loads((out/'render_data.json').read_text(encoding='utf-8'))
c=np.array(D['coeff484']);p=D['parameters'];q=p['station_advance_mm']
N=484; longitudinal=64; radial=16; radius=p['diameter_um']/2000
positions=[];normals=[];colors=[];faces=[]
for j in range(N):
    z=np.linspace(0,4*q,longitudinal+1)
    u=j+z/q; idx=np.floor(u).astype(int)%N;t=u-np.floor(u);a=c[idx]
    xy=((a[:,3]*t[:,None]+a[:,2])*t[:,None]+a[:,1])*t[:,None]+a[:,0]
    v=((3*a[:,3]*t[:,None]+2*a[:,2])*t[:,None]+a[:,1])/q
    tangent=np.column_stack([v,np.ones(len(v))]);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    n1=np.column_stack([np.ones(len(v)),np.zeros(len(v)),-v[:,0]]);n1/=np.linalg.norm(n1,axis=1)[:,None]
    n2=np.cross(tangent,n1)
    angles=np.arange(radial)*2*np.pi/radial
    nr=n1[:,None,:]*np.cos(angles)[None,:,None]+n2[:,None,:]*np.sin(angles)[None,:,None]
    ctr=np.column_stack([xy,z]);vs=(ctr[:,None,:]+radius*nr).reshape(-1,3)
    vs=np.vstack([vs,ctr[0],ctr[-1]])/1000 # glTF uses metres
    ns=np.vstack([nr.reshape(-1,3),-tangent[0],tangent[-1]])
    base=sum(len(x) for x in positions)
    ff=[]
    for k in range(longitudinal):
        for h in range(radial):
            aa=k*radial+h;bb=k*radial+(h+1)%radial;cc=aa+radial;dd=bb+radial
            ff.extend([[aa,bb,cc],[cc,bb,dd]])
    for h in range(radial):
        ff.append([(longitudinal+1)*radial,(h+1)%radial,h])
        ff.append([(longitudinal+1)*radial+1,longitudinal*radial+h,longitudinal*radial+(h+1)%radial])
    rgb=[75,136,188,255] if j==p['closeup']['loop_indices'][0] else [227,140,58,255] if j==p['closeup']['loop_indices'][1] else [148,155,160,255]
    positions.append(vs);normals.append(ns);colors.append(np.tile(rgb,(len(vs),1)));faces.append(np.array(ff)+base)
v=np.vstack(positions).astype('<f4');nn=np.vstack(normals).astype('<f4');col=np.vstack(colors).astype('u1');f=np.vstack(faces).astype('<u4')
# Check connectivity on one strand; every strand uses the same closed indexing.
f0=faces[0];edge=np.sort(np.vstack([f0[:,[0,1]],f0[:,[1,2]],f0[:,[2,0]]]),axis=1)
assert np.all(np.unique(edge,axis=0,return_counts=True)[1]==2)
buffers=[];views=[];off=0
for array,target in [(v,34962),(nn,34962),(col,34962),(f,34963)]:
    raw=array.tobytes();raw+=b'\0'*((-len(raw))%4)
    views.append({'buffer':0,'byteOffset':off,'byteLength':len(raw),'target':target});buffers.append(raw);off+=len(raw)
gltf={'asset':{'version':'2.0','generator':'Exact existing Litz spline -> finite copper tube visualization'},
 'scene':0,'scenes':[{'nodes':[0]}],'nodes':[{'mesh':0,'name':'484 copper strands, 4.132 mm segment; visualization mesh only'}],
 'meshes':[{'primitives':[{'attributes':{'POSITION':0,'NORMAL':1,'COLOR_0':2},'indices':3,'material':0}]}],
 'materials':[{'name':'Copper display colors','pbrMetallicRoughness':{'metallicFactor':.15,'roughnessFactor':.5}}],
 'buffers':[{'byteLength':off}],'bufferViews':views,
 'accessors':[{'bufferView':0,'componentType':5126,'count':len(v),'type':'VEC3','min':v.min(axis=0).tolist(),'max':v.max(axis=0).tolist()},
              {'bufferView':1,'componentType':5126,'count':len(nn),'type':'VEC3'},
              {'bufferView':2,'componentType':5121,'count':len(col),'type':'VEC4','normalized':True},
              {'bufferView':3,'componentType':5125,'count':f.size,'type':'SCALAR'}],
 'extras':{'copper_diameter_um':p['diameter_um'],'full_period_mm':500,'segment_length_mm':4*q,'insulation_included':False,
           'purpose':'Visualization; polygonal tubes are not a COMSOL swept-solid validation.','radial_facets':radial}}
j=json.dumps(gltf,separators=(',',':')).encode();j+=b' '*((-len(j))%4);binary=b''.join(buffers)
filename=out/'484-copper-tubes-local.glb'
filename.write_bytes(struct.pack('<4sII',b'glTF',2,12+8+len(j)+8+len(binary))+struct.pack('<I4s',len(j),b'JSON')+j+struct.pack('<I4s',len(binary),b'BIN\0')+binary)
print(json.dumps({'model':str(filename),'strands':N,'vertices':len(v),'triangles':len(f),'bytes':filename.stat().st_size,'closed_connectivity':True}))
