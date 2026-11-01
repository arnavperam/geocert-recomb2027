"""Render manuscript plots from the supplied, verified figure_data directory.

White backgrounds, open axes, explicit symbols and vector text. CSV files beside
the inputs can be imported into GraphPad Prism; no Prism software is required.
"""
from pathlib import Path
from fractions import Fraction as F
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap

ROOT=Path(__file__).resolve().parents[1];D=ROOT/'figure_data';O=ROOT/'figures';O.mkdir(exist_ok=True)
TEAL='#237E76';BRICK='#B65A40';BLUE='#426D94';BERRY='#A64F67';INK='#242424';GRAY='#737373';LIGHT='#D3D3D3'
plt.rcParams.update({'font.family':'Nimbus Sans','font.size':10.5,'axes.labelsize':10.5,'axes.titlesize':10.5,'axes.titleweight':'normal',
 'axes.linewidth':.85,'axes.edgecolor':INK,'axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,'xtick.direction':'out','ytick.direction':'out',
 'xtick.major.size':3.5,'ytick.major.size':3.5,'xtick.major.width':.85,'ytick.major.width':.85,'xtick.labelsize':10,'ytick.labelsize':10,
 'legend.frameon':False,'legend.fontsize':10,'legend.handlelength':2.2,'legend.borderaxespad':.3,'lines.linewidth':1.5,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','mathtext.fontset':'dejavusans'})
def read(n):return json.loads((D/(n+'.json')).read_text())
def csvread(n):return list(csv.DictReader((D/(n+'.csv')).open()))
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=220)
 plt.close(fig);print(name,flush=True)
def label(fig,ax,letter,title):
 x=ax.get_position().x0;y=ax.get_position().y1+.035
 fig.text(x-.045,y,letter,fontweight='bold',fontsize=12,va='bottom')
 fig.text(x+.006,y,title,fontsize=10.5,va='bottom')
def mask_display(ax,name,target=None):
 z=np.load(D/(name+'.npz'));m=z['mask'];extent=z['extent'];shade=np.full(m.shape,.965)
 if target is not None:shade[m==target]=.89
 edge=np.zeros_like(m,dtype=bool);edge[:-1]|=m[:-1]!=m[1:];edge[:,:-1]|=m[:,:-1]!=m[:,1:];shade[edge]=.75
 ax.imshow(shade,cmap='gray',vmin=0,vmax=1,extent=extent,interpolation='nearest',origin='upper',rasterized=True)
 ax.set_aspect('equal');return extent
def tissue(ax,legend=True):
 mask_display(ax,'ldhb_mask',81)
 for g in read('ldhb_geometry'):
  p=np.array(g['reference']);q=np.array(g['attaining'])
  ax.plot(p[:,0],p[:,1],color=INK,ls=(0,(4,2)),lw=1.1,zorder=3)
  ax.plot(q[:,0],q[:,1],color=TEAL,lw=1.9,zorder=4)
  for r in g['moved']:
   ax.plot(*r['xy'],marker='x' if r['gene']=='LDHB' else '+',color=BRICK if r['gene']=='LDHB' else INK,ms=5.2,mew=1.3,zorder=5)
 ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_xticks([35,40,45,50]);ax.set_yticks([150,155,160])
 if legend:
  handles=[Line2D([],[],color=INK,ls='--',label='Reference'),Line2D([],[],color=TEAL,label='Attaining'),Line2D([],[],color=BRICK,marker='x',ls='none',label='LDHB'),Line2D([],[],color=INK,marker='+',ls='none',label='KRT8')]
  ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.10,-.27),ncol=2,columnspacing=1.0,handlelength=1.5,fontsize=10)
def count_curves(ax,legend=True):
 rows=csvread('ldhb_intervals');x=np.array([float(r['budget_um2']) for r in rows])
 for side in ['lower','upper']:
  for prefix,color,ls,width in [('independent',GRAY,(0,(5,2)),1.5),('additive',BRICK,(0,(2,2)),1.6),('exact',TEAL,'-',1.8)]:
   ax.step(x,[int(r[prefix+'_'+side]) for r in rows],where='post',color=color,ls=ls,lw=width,zorder=4 if prefix=='exact' else 3)
 ax.fill_between(x,[int(r['exact_lower']) for r in rows],[int(r['exact_upper']) for r in rows],step='post',color=TEAL,alpha=.08,zorder=0)
 ax.set(xlim=(0,.65),ylim=(16.5,23.5),xticks=[0,.25,.5],xticklabels=['0','0.25','0.50'],yticks=[17,19,21,23],xlabel='Area budget (µm²)',ylabel='LDHB count')
 if legend:
  handles=[Line2D([],[],color=TEAL,label='Exact'),Line2D([],[],color=GRAY,ls=(0,(5,2)),label='Independent'),Line2D([],[],color=BRICK,ls=(0,(2,2)),label='Additive')]
  ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.10,-.27),ncol=2,columnspacing=.9,handlelength=1.8,fontsize=10)

fig=plt.figure(figsize=(6.5,3.15));a=fig.add_axes([.075,.33,.34,.53]);b=fig.add_axes([.56,.33,.42,.53]);tissue(a);count_curves(b)
label(fig,a,'a','An attaining boundary');label(fig,b,'b','LDHB count intervals');save(fig,'ldhb_composite')
fig,ax=plt.subplots(figsize=(6.5,3.3));fig.subplots_adjust(left=.12,right=.98,bottom=.27,top=.96);count_curves(ax,False)
fig.legend(handles=[Line2D([],[],color=TEAL,label='Exact'),Line2D([],[],color=GRAY,ls=(0,(5,2)),label='Independent outer bound'),Line2D([],[],color=BRICK,ls=(0,(2,2)),label='Additive inner approximation')],loc='lower center',ncol=2,bbox_to_anchor=(.52,.005));save(fig,'count_curve')

donors={r['donor']:r for r in read('donors')};regions=read('regions');budgets=read('donor_budgets');colors={'IDC':BLUE,'ILC':BERRY}
fig=plt.figure(figsize=(6.5,2.8));a=fig.add_axes([.10,.25,.30,.60]);b=fig.add_axes([.57,.25,.40,.60])
for j,donor in enumerate(['IDC','ILC']):
 rs=[r for r in regions if r['donor']==donor];col=colors[donor]
 a.scatter(j+np.array([-.12,0,.12]),[100*r['tightened']/r['eligible'] for r in rs],s=25,facecolor='white',edgecolor=col,linewidth=1.1,zorder=3)
 pooled=100*donors[donor]['tightened']/donors[donor]['eligible'];a.scatter([j],[pooled],marker='D',s=35,color=col,zorder=4)
 a.text(j,18.7,f'{pooled:.1f}%',color=col,ha='center',fontsize=10.5)
 rs=sorted([r for r in budgets if r['donor']==donor and F(r['budget_um2'])>0],key=lambda r:F(r['budget_um2']))
 b.plot([float(F(r['budget_um2'])) for r in rs],[100*r['tightened']/donors[donor]['eligible'] for r in rs],marker='o' if j==0 else 's',ms=4,mfc='white',mew=1.1,color=col,label='Ductal' if j==0 else 'Lobular')
a.set(xlim=(-.45,1.45),ylim=(0,21),yticks=[0,5,10,15,20],xticks=[0,1],xticklabels=['Ductal','Lobular'],ylabel='Eligible pairs tightened (%)')
b.set(xscale='log',xlim=(.04,5),ylim=(-.3,8.5),yticks=[0,2,4,6,8],xticks=[.05,.1,.25,.5,1,2,4],xticklabels=['.05','.1','.25','.5','1','2','4'],xlabel='Area budget (µm²)',ylabel='Eligible pairs tightened (%)');b.minorticks_off();b.legend(loc='upper right',fontsize=10)
label(fig,a,'a','Across the budget grid');label(fig,b,'b','At each tested budget')
fig.text(.08,.035,'Open circles: individual regions     Diamonds: pooled donor fractions',fontsize=10);save(fig,'donor_replication')

width=read('width_summary');wd={(r['dataset'],r['setting']):r for r in width};order=['baseline','.5','1','2','3'];order=['baseline','0.5','1','2','3']
series=[('lung','proseg','Lung / Proseg',BLUE,'o'),('lung','nuclear','Lung / nuclear',TEAL,'s'),('breast','nuclear','Breast / nuclear',BRICK,'^')]
union={(r['dataset'],r['reference']):r for r in read('union_reach') if r['budget_um2']=='unbounded'}
fig=plt.figure(figsize=(6.5,2.9));a=fig.add_axes([.10,.30,.36,.55]);b=fig.add_axes([.64,.30,.33,.55])
for ds,ref,name,color,mark in series:
 rows=[wd[ds,k] for k in order];xs=[100*r['edge_coverage'] for r in rows];ys=[100*r[ref+'_asymptotic_reach'] for r in rows]
 a.plot(xs,ys,color=color,marker=mark,ms=4,mfc='white',mew=1.1,lw=1.1,label=name)
 # Label support widths once, along the breast series and the Proseg series.
 if ds=='breast' or ref=='proseg':
  for x,y,k in zip(xs,ys,order):a.annotate('A' if k=='baseline' else k,(x,y),xytext=(0,5 if ds=='breast' else -12),textcoords='offset points',ha='center',fontsize=10)
 r=union[ds,ref];b.plot([0,1],[100*r['baseline_reach'],100*r['union_reach']],color=color,marker=mark,ms=5,mfc='white',mew=1.2,lw=1.1)
 b.annotate(f"{100*r['union_reach']:.1f}%",(1,100*r['union_reach']),xytext=(6,4 if ds=='lung' and ref=='nuclear' else -8 if ds=='breast' else 0),textcoords='offset points',color=color,fontsize=10)
a.set(xlim=(0,78),ylim=(0,24),xticks=[0,20,40,60],yticks=[0,10,20],xlabel='Editable boundary (%)',ylabel='Reference reach (%)')
b.set(xlim=(-.2,1.7),ylim=(0,32),xticks=[0,1],xticklabels=['Primary','Union'],yticks=[0,10,20,30],ylabel='Reference reach (%)')
label(fig,a,'a','Fixed-width atlases');label(fig,b,'b','Alternative complete atlases')
fig.legend(*a.get_legend_handles_labels(),loc='lower center',ncol=3,bbox_to_anchor=(.52,.015),columnspacing=1,handletextpad=.5);save(fig,'reach_composite')

wb=read('width_budget');palette=[INK,'#B8B8B8',BLUE,TEAL,BRICK];marks=['D','o','s','^','v']
fig,axs=plt.subplots(1,2,figsize=(6.5,3.0));fig.subplots_adjust(left=.10,right=.98,bottom=.27,top=.86,wspace=.34)
for ax,ds in zip(axs,['lung','breast']):
 for k,col,mark in zip(order,palette,marks):
  rows=sorted([r for r in wb if r['dataset']==ds and r['setting']==k and F(r['budget_um2'])>0],key=lambda r:F(r['budget_um2']))
  ax.plot([float(F(r['budget_um2'])) for r in rows],[100*r['non_singleton']/r['questions'] for r in rows],color=col,marker=mark,ms=3,mfc='white',mew=.8,lw=1.1,label='Primary' if k=='baseline' else k+' µm')
 ax.set(xscale='log',xlim=(.004,20),ylim=(0,30),yticks=[0,10,20,30],xticks=[.01,.1,1,10],xticklabels=['.01','.1','1','10'],xlabel='Area budget (µm²)',ylabel='Non-singleton pairs (%)');ax.minorticks_off()
label(fig,axs[0],'a','Lung');label(fig,axs[1],'b','Breast format-test')
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=5,bbox_to_anchor=(.53,.015),columnspacing=.7,handlelength=1.1,handletextpad=.4);save(fig,'width_budget')

chain=read('chain_scaling');fig,ax=plt.subplots(figsize=(6.5,2.9));fig.subplots_adjust(left=.12,right=.97,bottom=.23,top=.95)
for n in [16,32,64,128]:
 vals=[r['total_seconds'] for r in chain['records'] if r['n']==n];assert len(vals)==3
 ax.scatter(n*np.array([.97,1,1.03]),vals,s=23,facecolor='white',edgecolor=GRAY,lw=1,zorder=3)
ax.plot([r['n'] for r in chain['table']],[r['median_total_seconds'] for r in chain['table']],color=TEAL,marker='D',ms=4,label='Median of three runs')
ax.set(xscale='log',yscale='log',xticks=[16,32,64,128],xticklabels=['16','32','64','128'],yticks=[.5,1,5,10,50],yticklabels=['0.5','1','5','10','50'],xlabel='Explicit point transcripts',ylabel='Runtime with witnesses (s)');ax.minorticks_off();ax.legend(loc='upper left');save(fig,'chain_scaling')

# Source tissue overview, without arbitrary colors assigned to cell identities.
fig,ax=plt.subplots(figsize=(6.5,3.4));fig.subplots_adjust(left=.10,right=.98,bottom=.18,top=.90);mask_display(ax,'lung_excerpt')
r=read('lung_excerpt');ax.add_collection(LineCollection(r['paths'],colors=TEAL,linewidths=1.05));pts=np.array(r['rna']);inside=(pts[:,0]>=45)&(pts[:,0]<=125)&(pts[:,1]>=0)&(pts[:,1]<=44)
ax.scatter(*pts[inside].T,s=4,c=INK,alpha=.60,linewidths=0)
ax.plot([52,62],[39,39],color=INK,lw=2);ax.text(57,36,'10 µm',ha='center',fontsize=10)
ax.set(xlim=(45,125),ylim=(44,0),xlabel='x (µm)',ylabel='y (µm)')
fig.legend(handles=[Line2D([],[],color=TEAL,label='Editable source arcs'),Line2D([],[],color=INK,marker='o',ms=3,ls='none',label='Observed RNA')],ncol=2,loc='upper center',bbox_to_anchor=(.54,1.0));save(fig,'atlas_reference')

# Each line connects the same core identifiers, not independent replicates.
context=read('proseg_context');fig,axs=plt.subplots(1,2,figsize=(6.5,3.2));fig.subplots_adjust(left=.10,right=.98,bottom=.22,top=.84,wspace=.40)
for ax,dslist,title in zip(axs,[['LUNG'],['IDC','ILC']],['Lung cores','Breast donor cores']):
 for ds in dslist:
  regions=sorted({r['region'] for r in context if r['donor']==ds})
  for i,reg in enumerate(regions):
   rs=sorted([r for r in context if r['donor']==ds and r['region']==reg and r['factorization']=='no_factorization'],key=lambda r:r['input_um']);assert len(rs)==2 and all(r['returncode']==0 for r in rs)
   ys=[r['background_percent'] for r in rs];col=[BLUE,BERRY][i] if ds=='LUNG' else colors[ds]
   ax.plot([0,1],ys,color=col,marker='o' if ds!='ILC' else 's',mfc='white',ms=4,mew=1,lw=1.1)
   if ds=='LUNG':
    for j,y in enumerate(ys):ax.annotate(f'{y:.1f}%',(j,y),xytext=(0,7),textcoords='offset points',ha='center',color=col,fontsize=10)
 ax.set(xlim=(-.25,1.25),xticks=[0,1],xticklabels=['200 µm','800 µm'],xlabel='Input window side',ylabel='Core RNA in background (%)',ylim=(0,105 if dslist==['LUNG'] else 35))
label(fig,axs[0],'a','Lung cores');label(fig,axs[1],'b','Breast donor cores')
axs[1].legend(handles=[Line2D([],[],color=BLUE,marker='o',mfc='white',label='Ductal'),Line2D([],[],color=BERRY,marker='s',mfc='white',label='Lobular')],loc='upper right',fontsize=10);save(fig,'proseg_context')

def export_csv(name, rows):
 with (D/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

# Flat plotting tables for Prism or another editor; denominators are explicit.
export_csv('donor_regions_plot', [dict(donor=r['donor'],region=r['region'],eligible_pairs=r['eligible'],tightened_pairs=r['tightened'],tightened_percent=100*r['tightened']/r['eligible']) for r in read('regions')])
export_csv('donor_budget_plot', [dict(donor=r['donor'],budget_um2=float(F(r['budget_um2'])),eligible_pairs_full_grid=donors[r['donor']]['eligible'],tightened_pairs=r['tightened'],tightened_percent=100*r['tightened']/donors[r['donor']]['eligible']) for r in budgets])
export_csv('width_budget_plot', [dict(dataset=r['dataset'],support_halfwidth_um=r['setting'],budget_um2=float(F(r['budget_um2'])),cell_gene_pairs=r['questions'],non_singleton_pairs=r['non_singleton'],non_singleton_percent=100*r['non_singleton']/r['questions']) for r in wb])
export_csv('chain_scaling_plot', [dict(explicit_RNA=r['n'],run=i%3+1,total_seconds=r['total_seconds']) for i,r in enumerate(chain['records'])])
export_csv('proseg_context_plot', [dict(donor=r['donor'],region=r['region'],input_window_um=r['input_um'],common_core_RNA=r['common_retained_core_RNA'],background_core_RNA=r['background_core_RNA'],background_percent=r['background_percent']) for r in context if r['factorization']=='no_factorization' and r['returncode']==0])
print('All quantitative figures rendered from recorded observations and results.')
