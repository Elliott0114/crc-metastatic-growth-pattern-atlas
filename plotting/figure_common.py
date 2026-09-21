"""Fixed-size views of frozen results. No fitting or result-dependent selection."""
from pathlib import Path
import hashlib
import json
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import QuadMesh
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parents[1]
Q = Path('analysis_results/rec_q2_integration_2026-09-19')
N = Path('analysis_results/rec_manuscript_reinforcement_2026-09-17')
B = Path('analysis_results/rec_evidence_bridge_2026-09-16')
C = Path('analysis_results/rec_junction_context_2026-09-16')
P = Path('manuscript/rec_discovery_reinforced_2026-09-17/source_data/preserved_figure_inputs')
COL = dict(dHGP='#3B6FB6', rHGP='#CB7144', Claudin='#218C84', Polarity='#8064A2',
           F11R='#B38C38', HRC='#336B60', Junction='#555555')
REGION = {'Other retained spots':'#C7CDD1', 'Liver-like':'#61A89E',
          'Near tumour-side':'#D89958', 'Deep tumour-side':'#8064A2'}
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'axes.labelsize':8,
    'axes.titlesize':9, 'xtick.labelsize':8, 'ytick.labelsize':8, 'legend.fontsize':8,
    'axes.spines.top':False, 'axes.spines.right':False, 'axes.linewidth':.65,
    'xtick.major.width':.6, 'ytick.major.width':.6, 'lines.linewidth':.8,
    'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none', 'savefig.facecolor':'white'})
for folder in ['figures','editable_figures','review','source_data/panels','source_data/frozen','revision']:
    (OUT/folder).mkdir(parents=True,exist_ok=True)
MANIFEST = OUT/'panel_manifest.json'
RECORDS = json.loads(MANIFEST.read_text(encoding='utf-8')) if MANIFEST.exists() else []
ACTIVE = None
SOURCES = {}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def begin(name):
    global ACTIVE, SOURCES, RECORDS
    ACTIVE=name; SOURCES={}; RECORDS=[r for r in RECORDS if r['figure_id']!=name]

def source(path):
    path=Path(path)
    relative=path.relative_to(ROOT) if path.is_absolute() else path
    actual=ROOT/relative
    if not actual.is_file(): raise FileNotFoundError(actual)
    digest=sha(actual)
    SOURCES[str(relative)]=digest
    return actual

def read(path):
    return pd.read_csv(source(path),sep='\t')

def panel_data(letter,frame,**meta):
    path=OUT/'source_data/panels'/f'{ACTIVE}_{letter}.tsv'
    frame.to_csv(path,sep='\t',index=False,encoding='utf-8',na_rep='NA')
    base={'figure_id':ACTIVE,'panel_id':letter,'old_panel':'new encoding of frozen results',
          'source_family':'see source study/patient columns','patient_set':'see source data',
          'score_version':'frozen, as labelled','gene_set':'frozen, as labelled',
          'normalization':'as labelled','model':'as labelled','contrast':'as labelled',
          'estimand':'as labelled','unit_of_inference':'patient',
          'ci_method':'frozen nominal 95% intervals where applicable',
          'multiplicity_family':'unchanged; see result q columns and legend',
          'source_result_paths':json.dumps(list(SOURCES),ensure_ascii=False),
          'source_hashes':json.dumps(SOURCES,sort_keys=True),
          'panel_source_data':str(path.relative_to(OUT)),'rows':len(frame),
          'output_path':f'figures/{ACTIVE}.pdf','notes':''}
    base.update(meta);RECORDS.append(base)

def canvas(height=175):
    fig=plt.figure(figsize=(170/25.4,height/25.4));fig._height_mm=height
    return fig

def ax(fig,x,top,width,height):
    return fig.add_axes([x/170,1-(top+height)/fig._height_mm,width/170,height/fig._height_mm])

def text(fig,x,top,value,**kwargs):
    return fig.text(x/170,1-top/fig._height_mm,value,va=kwargs.pop('va','top'),fontsize=kwargs.pop('fontsize',8),**kwargs)

def heading(fig,letter,title,x,top):
    text(fig,x,top,letter,fontweight='bold',fontsize=11)
    text(fig,x+6,top+.2,title,fontweight='bold',fontsize=9)

def zero(axis):
    axis.axvline(0,color='.55',ls=(0,(3,3)),lw=.8,zorder=0)

def patient_offsets(identifiers,width=.14):
    ids=[str(x) for x in identifiers]
    ordered=sorted(set(ids))
    mapping=dict(zip(ordered,np.linspace(-width,width,len(ordered))))
    return np.array([mapping[x] for x in ids])

def hgp_handles():
    return [Line2D([],[],ls='',marker='o' if h=='dHGP' else 's',color=COL[h],
                   markersize=3.7,label=h) for h in ['dHGP','rHGP']]

def forest(axis,frame,labels,estimate='effect',low='low',high='high',color='.25',marker='o'):
    for i,(_,r) in enumerate(frame.iterrows()):
        if pd.isna(r[estimate]):
            axis.text(.97,i,'NA',transform=axis.get_yaxis_transform(),ha='right',va='center',color='.45')
            continue
        co=color[i] if isinstance(color,list) else color
        if pd.notna(r[low]) and pd.notna(r[high]):axis.plot([r[low],r[high]],[i,i],c=co,lw=1)
        axis.plot(r[estimate],i,marker,c=co,ms=3.7)
    axis.set(yticks=range(len(frame)),yticklabels=labels,ylim=(len(frame)-.5,-.5))
    zero(axis);axis.xaxis.set_major_locator(MaxNLocator(4))

def groups(axis,frame,column,ylabel,percent=False):
    for i,h in enumerate(['dHGP','rHGP']):
        v=frame.loc[frame.hgp.eq(h),column].to_numpy()*(100 if percent else 1)
        axis.scatter(i+np.linspace(-.12,.12,len(v)),v,s=15,c=COL[h],marker='o' if h=='dHGP' else 's',zorder=3)
        axis.plot([i-.21,i+.21],[v.mean()]*2,c='.2',lw=1.2)
    axis.set(xticks=[0,1],xticklabels=[f'dHGP ({sum(frame.hgp.eq("dHGP"))})',f'rHGP ({sum(frame.hgp.eq("rHGP"))})'],xlim=(-.45,1.45),ylabel=ylabel)
    axis.yaxis.set_major_locator(MaxNLocator(4))

def members():
    d=read(C/'definitions.tsv')
    return sorted(d.loc[d.component.eq('Claudin'),'gene'],key=lambda g:int(g[4:]))+['CRB3','PALS1','PARD3','PARD6A','PARD6B','PARD6G','PATJ','PRKCI','F11R']

def member_color(g):
    return COL['Claudin' if g.startswith('CLDN') else 'F11R' if g=='F11R' else 'Polarity']

def row_lines(axis):
    for y in [20.5,28.5]:axis.axhline(y,c='.65',lw=.65)

def heat(axis,matrix,vmin=None,vmax=None,cmap='PuOr_r',missing='NA'):
    data=np.asarray(matrix,dtype=float)
    if vmax is None:vmax=float(np.nanmax(np.abs(data))) or 1
    if vmin is None:vmin=-vmax
    cm=plt.colormaps[cmap].copy();cm.set_bad('#DFDFDF')
    ny,nx=data.shape
    im=axis.pcolormesh(np.arange(nx+1)-.5,np.arange(ny+1)-.5,np.ma.masked_invalid(data),
                      cmap=cm,vmin=vmin,vmax=vmax,shading='flat',edgecolors='none',rasterized=False)
    axis.set(xlim=(-.5,nx-.5),ylim=(ny-.5,-.5))
    axis.tick_params(length=0)
    # Masked QuadMesh cells are transparent; draw the original missing-value grey.
    for y,x in zip(*np.where(~np.isfinite(data))):
        axis.add_patch(Rectangle((x-.5,y-.5),1,1,facecolor='#DFDFDF',edgecolor='none',zorder=-1))
    for y,x in zip(*np.where(~np.isfinite(data))):axis.text(x,y,missing,ha='center',va='center',fontsize=8,c='.4')
    return im

def cbar(fig,im,x,top,width,label,ticks=None):
    cb=fig.colorbar(im,cax=ax(fig,x,top,width,2.4),orientation='horizontal',ticks=ticks)
    cb.ax.tick_params(labelsize=8,length=2,pad=1)
    cb.set_label(label,labelpad=1,fontsize=8)
    return cb

def model_legend(fig,x,top):
    handles=[Line2D([],[],ls='',marker='o',mfc='white',mec='.5',label='Baseline'),
             Line2D([],[],ls='',marker='s',color='.15',label='Source-adjusted')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(x/170,1-top/fig._height_mm),
               ncol=2,frameon=False,borderaxespad=0,handletextpad=.4,columnspacing=1)

def regional_forest(axis,frame,genes,labels=True):
    for j,(model,co,mark) in enumerate([('unadjusted','.5','o'),('source_adjusted','.15','s')]):
        d=frame[frame.model.eq(model)].set_index('gene').reindex(genes)
        for i,(_,r) in enumerate(d.iterrows()):
            yy=i+(-.15 if j==0 else .15)
            if pd.isna(r.logFC):continue
            axis.plot([r.approximate_QL_Wald_low,r.approximate_QL_Wald_high],[yy,yy],color=co,lw=.85)
            axis.plot(r.logFC,yy,mark,ms=3.6,mec=co,mfc='white' if j==0 else co)
    axis.set(yticks=range(len(genes)),yticklabels=genes if labels else [],ylim=(len(genes)-.5,-.5))
    zero(axis);axis.xaxis.set_major_locator(MaxNLocator(4))

def save(fig,name):
    for axis in fig.axes:
        for collection in axis.collections:
            if isinstance(collection,QuadMesh):collection.set_rasterized(False)
    fig.canvas.draw()
    # Matplotlib retains Text objects for ticks beyond the view; those are not drawn.
    hidden_ticks=set()
    for axis in fig.axes:
        for dimension,limits in [(axis.xaxis,axis.get_xlim()),(axis.yaxis,axis.get_ylim())]:
            lo,hi=sorted(limits)
            for tick in dimension.get_major_ticks()+dimension.get_minor_ticks():
                if not lo-1e-9<=tick.get_loc()<=hi+1e-9:hidden_ticks.update([id(tick.label1),id(tick.label2)])
    labels=[]
    for obj in fig.findobj(matplotlib.text.Text):
        if not obj.get_visible() or not obj.get_text() or id(obj) in hidden_ticks:continue
        box=obj.get_window_extent(fig.canvas.get_renderer())
        labels.append({'text':obj.get_text(),'size_pt':obj.get_fontsize(),'x0':box.x0,'y0':box.y0,'x1':box.x1,'y1':box.y1})
    for folder,extension in [('figures','pdf'),('editable_figures','svg'),('review','png')]:
        fig.savefig(OUT/folder/f'{name}.{extension}',dpi=300)
    (OUT/'review'/f'{name}_text_bounds.json').write_text(json.dumps({'width_px':fig.bbox.width,'height_px':fig.bbox.height,'labels':labels},ensure_ascii=False,indent=2),encoding='utf-8')
    axes=[{'x_mm':a.get_position().x0*170,'y_from_bottom_mm':a.get_position().y0*fig._height_mm,
           'width_mm':a.get_position().width*170,'height_mm':a.get_position().height*fig._height_mm,
           'xlabel':a.get_xlabel(),'ylabel':a.get_ylabel(),
           'xlim':list(a.get_xlim()),'ylim':list(a.get_ylim())} for a in fig.axes]
    (OUT/'review'/f'{name}_axes.json').write_text(json.dumps(axes,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    from build_review import caption_page
    caption_page(fig,name)
    MANIFEST.write_text(json.dumps(RECORDS,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    pd.DataFrame(RECORDS).to_csv(OUT/'panel_manifest.tsv',sep='\t',index=False,encoding='utf-8')
    plt.close(fig);print('Rendered',name,flush=True)
