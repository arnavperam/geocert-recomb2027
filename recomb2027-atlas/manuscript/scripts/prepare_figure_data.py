"""Extract plotting inputs from the published, restored results. No experiments run.

The LDHB singleton thresholds are evaluated from the archived geometry and checked
against every saved budget row. All other values are copied from result records.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, csv, hashlib, json, sys
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data-root',type=Path,required=True)
p.add_argument('--public-root',type=Path,required=True)
a=p.parse_args();R=a.data_root;P=a.public_root
O=Path(__file__).resolve().parents[1]/'figure_data';O.mkdir(exist_ok=True)
sys.path.insert(0,str(R))
from atlas.model import compose
from atlas.scalar_reference import prepare,query
from atlas.transport import chart_to_pixel
inputs={}
def read(path):
 inputs[str(path.relative_to(R) if path.is_relative_to(R) else path.relative_to(P))]=hashlib.sha256(path.read_bytes()).hexdigest()
 return json.loads(path.read_text())
def save(name,records):
 (O/(name+'.json')).write_text(json.dumps(records,indent=2)+'\n')
def csvsave(name,rows):
 with (O/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

D=R/'external';atlas=read(D/'data/atlas_breast_01.json')
spectra=read(D/'results/cell_spectra.json')
s=next(q for q in spectra if q['roi']=='breast_01' and q['cell']==81 and q['gene']=='LDHB')
local={(q['patch'],q['gene']):q for q in read(D/'results/patch_spectra.json') if q['roi']=='breast_01'}
adds={(q['patch'],q['gene']):q for q in read(D/'results/additive_local_spectra.json') if q['roi']=='breast_01'}
pm={q['id']:q for q in atlas['patches']};pixel=float(F(atlas['pixel_size_um']))
maskpath=D/'data/source_crops/breast_01/masks.npz';inputs[str(maskpath.relative_to(R))]=hashlib.sha256(maskpath.read_bytes()).hexdigest()
m=np.load(maskpath)['cell_mask'];yy,xx=np.nonzero(m==81)
x0,x1=max(0,xx.min()-15),min(m.shape[1],xx.max()+16);y0,y1=max(0,yy.min()-15),min(m.shape[0],yy.max()+16)
np.savez_compressed(O/'ldhb_mask.npz',mask=m[y0:y1,x0:x1],extent=np.array([x0,x1,y1,y0])*pixel)
geometry=[];single=[];parts=[]
for pid in s['patches']:
 patch=pm[pid];pre=prepare(patch)
 for i,r in enumerate(patch['points']):
  if r['gene']!='LDHB':continue
  answer=query(patch,{i:1-r['source']},pre)
  if answer is not None:single.append({'cost':str(answer[1]),'out':patch['cell_labels'][r['source']]==81})
 ad=adds[pid,'LDHB'];sign=1 if patch['cell_labels'][0]==81 else -1
 parts.append((pid,{str((int(k)-ad['source_left'])*sign):v for k,v in ad['additive_minima'].items()}))
 delta=s['witness_traces']['17'][pid]
 if delta==0:continue
 sp=local[pid,'LDHB'];left=sp['source_left']+delta*sign;w=sp['witnesses'][str(left)]
 seg=[list(map(F,v)) for v in w['boundaries'][0]];ts={u for v in seg for u in v[:2]}
 sx,sy=patch['sx'],patch['sy'];org=sx*patch['source_path_pixels'][0][0]+sy*patch['source_path_pixels'][0][1]
 ts|={F(sx*x+sy*y-org) for x,y in patch['source_path_pixels']}
 def at(t):
  for u,v,x,y in seg:
   if u<=t<=v:return x+(y-x)*(t-u)/(v-u)
  raise ValueError(t)
 alt=[[float(v)*pixel for v in chart_to_pixel(patch,t,at(t))] for t in sorted(ts)]
 moved=[]
 for r,owner in zip(patch['points'],w['owners']):
  if owner!=r['source']:
   moved.append({'gene':r['gene'],'id':r['transcript_id'],'xy':[float(v)*pixel for v in chart_to_pixel(patch,F(r['t']),F(r['x']))]})
 geometry.append({'patch':pid,'delta':delta,'area_um2':w['physical_area_um2'],'reference':(np.asarray(patch['source_path_pixels'])*pixel).tolist(),'attaining':alt,'moved':moved})
save('ldhb_geometry',geometry)
exact={int(k):F(v) for k,v in s['minima'].items()};add={k:v for k,(v,tr) in compose(parts,22).items()}
def endpoints(b):
 e=[k for k,v in exact.items() if v<=b];ad=[k for k,v in add.items() if v<=b]
 return dict(budget_um2=float(b),budget_exact=str(b),exact_lower=min(e),exact_upper=max(e),independent_lower=22-sum(r['out'] and F(r['cost'])<=b for r in single),independent_upper=22+sum(not r['out'] and F(r['cost'])<=b for r in single),additive_lower=min(ad),additive_upper=max(ad))
thresholds={F(0),F('0.65'),F('.25'),F('.5')}|set(exact.values())|set(add.values())|{F(r['cost']) for r in single}
rows=[endpoints(b) for b in sorted(thresholds) if b<=F('.65')]
reference=D/'results/count_baseline_comparison.csv';inputs[str(reference.relative_to(R))]=hashlib.sha256(reference.read_bytes()).hexdigest()
checked=0
for r in csv.DictReader(reference.open()):
 if r['roi']=='breast_01' and r['cell']=='81' and r['gene']=='LDHB':
  got=endpoints(F(r['budget_um2']))
  for source,target in [('lower','exact_lower'),('upper','exact_upper'),('independent_lower','independent_lower'),('independent_upper','independent_upper'),('additive_inner_lower','additive_lower'),('additive_inner_upper','additive_upper')]:assert int(r[source])==got[target]
  checked+=1
assert checked==8
assert sum(x['gene']=='LDHB' for g in geometry for x in g['moved'])==5
assert sum(x['gene']=='KRT8' for g in geometry for x in g['moved'])==1
csvsave('ldhb_intervals',rows)
save('ldhb_verification',{'saved_budget_rows_matched':checked,'source_count':22,'attaining_count':17,'minimum_area_um2':str(exact[17]),'source_spectrum':s,'singleton_thresholds':single})

for name,path in [('width_summary',R/'results/reach/summary.json'),('width_budget',R/'results/reach/budget_summary.json'),('chain_scaling',R/'results/chain_benchmark_summary.json'),('donors',P/'extensions/external_validation/results/donor_summary.json'),('regions',P/'extensions/external_validation/results/region_summary.json'),('donor_budgets',P/'extensions/external_validation/results/budget_summary.json'),('proseg_context',P/'extensions/external_validation/results/proseg_summary.json'),('union_reach',P/'extensions/atlas_union/results/reference_summary.json')]:
 obj=read(path);save(name,obj)
 if isinstance(obj,list) and obj and name!='proseg_context':csvsave(name,obj)

# A native lung-mask excerpt, retaining its recorded tissue coordinates.
atlas=read(R/'data/atlas_roi_01.json');pixel=float(F(atlas['pixel_size_um']))
mp=R/'data/source_crops/roi_01/masks.npz';inputs[str(mp.relative_to(R))]=hashlib.sha256(mp.read_bytes()).hexdigest();m=np.load(mp)['cell_mask']
x0,x1=int(45/pixel),int(125/pixel)+1;y0,y1=0,int(44/pixel)+1
np.savez_compressed(O/'lung_excerpt.npz',mask=m[y0:y1,x0:x1],extent=np.array([x0,x1,y1,y0])*pixel)
ox,oy=atlas['bbox_pixels'][:2]
save('lung_excerpt',{'paths':[(np.asarray(q['source_path_pixels'])*pixel).tolist() for q in atlas['patches']], 'rna':[[float(F(r['x_um']))-ox*pixel,float(F(r['y_um']))-oy*pixel] for r in atlas['records']]})
save('provenance',{'source_repository':'https://github.com/arnavperam/geocert-recomb2027','distribution_version':'1.1.1','inputs_sha256':inputs,'operation':'Extract saved results and render; no new experimental runs','ldhb_saved_budget_rows_matched':checked})
print('Plotting inputs extracted; LDHB endpoints match all eight saved budgets.')
