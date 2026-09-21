"""Five redesigned main figures; only frozen estimates and original image data."""
from figure_common import *

STUDIES=['Che_2021_Cell_Discov','Liu_2024_Cancer_Res','Wang_2023_Sci_Adv','Giguelay_2022_Theranostics','Sathe_2023_Clin_Cancer_Res']

def paired_components(axis,d):
    offsets=patient_offsets(d.study.astype(str)+'|'+d.patient.astype(str))
    for offset,(_,r) in zip(offsets,d.iterrows()):axis.plot(np.array([0,1])+offset,[r.rho_claudin,r.rho_polarity],c='.79',lw=.65,zorder=1)
    for x,k in enumerate(['Claudin','Polarity']):axis.scatter(x+offsets,d['rho_'+k.lower()],s=15,c=COL[k],zorder=3)
    axis.axhline(0,c='.7',lw=.65,ls='--')
    axis.set(xticks=[0,1],xticklabels=['Claudin','Polarity'],xlim=(-.35,1.35),ylabel='Partial Spearman ρ')

def figure2():
    name='Figure_2';begin(name)
    pts=read(Q/'external/patient_component_effects.tsv').query("scoring=='fixed_background' and model=='common'")
    study=read(Q/'external/study_component_effects.tsv').query("scoring=='fixed_background' and model=='common' and endpoint=='Claudin_minus_Polarity'")
    meta=read(Q/'external/random_effects_meta.tsv').query("scoring=='fixed_background' and model=='common' and endpoint=='Claudin_minus_Polarity'")
    matching=read(N/'cells/comparison_summary.tsv')
    fig=canvas(170)
    a=ax(fig,22,21,53,43);b=ax(fig,111,21,53,43)
    og=pts[pts.study.eq('Ogden_2025')];external=pts[pts.study.isin(STUDIES)]
    paired_components(a,og);paired_components(b,external)
    lim=(pts[['rho_claudin','rho_polarity']].min().min()-.04,pts[['rho_claudin','rho_polarity']].max().max()+.04)
    a.set_ylim(lim);b.set_ylim(lim)
    heading(fig,'A','Discovery: common model',5,5);heading(fig,'B','Independent malignant cells',92,5)
    text(fig,24,14,'14 patients');text(fig,112,14,'28 patients / 5 studies')
    panel_data('A',og,old_panel='Main 2a',source_family='Ogden 2025',patient_set='14 patients',score_version='harmonized fixed background',model='common',estimand='partial Spearman rho')
    panel_data('B',external,old_panel='Main 2b',source_family='five original studies; Atlas container not counted',patient_set='28 patients',score_version='harmonized fixed background',model='common',estimand='partial Spearman rho')
    c=ax(fig,25,96,66,56)
    d=study.set_index('study').loc[STUDIES].reset_index()
    pooled=meta.iloc[0];joined=pd.concat([d,meta.assign(study='Random effects')],ignore_index=True)
    labels=[f'{s.split("_")[0]} ({int(n)})' for s,n in zip(d.study,d.n_patients)]+['Pooled (28)']
    forest(c,joined,labels,estimate='estimate',color=[COL['Claudin']]*5+['black'])
    c.lines[-2].set_marker('D')
    c.set_ylim(7.25,-.5)
    c.set_xlim(joined.low.min()-.02,joined.high.max()+.02)
    c.set_xticks([-.1,0,.1,.2])
    c.set_xlabel('Claudin − Polarity (Δ Fisher z)')
    heading(fig,'C','Study and pooled differences',5,83)
    c.text(.02,5.78,f'Δz = {pooled.estimate:.3f}; 95% CI {pooled.low:.3f}–{pooled.high:.3f}',transform=c.get_yaxis_transform(),va='center',bbox={'facecolor':'white','edgecolor':'none','pad':.5})
    c.text(.02,6.48,f'P = {pooled.p:.3f}; I² = {pooled.I2:.1f}%',transform=c.get_yaxis_transform(),va='center',bbox={'facecolor':'white','edgecolor':'none','pad':.5})
    panel_data('C',joined,old_panel='Main 2c',model='common; study-wise random effects',estimand='Claudin minus Polarity Fisher-z difference',ci_method='study patient bootstrap; pooled Hartung–Knapp',patient_set='28 patients / 5 studies')
    d_ax=ax(fig,124,108,40,37)
    for i,scope in enumerate(['full_reference','size_matched','measurement_matched']):
        for offset,mode,mark,co in [(-.13,'fixed_background','o',COL['Claudin']),(.13,'target_mean','s','.45')]:
            r=matching[(matching.scope==scope)&(matching.scoring==mode)].iloc[0]
            d_ax.plot([r.low,r.high],[i+offset]*2,c=co,lw=1)
            d_ax.plot(r.effect,i+offset,mark,c=co,ms=3.5)
    zero(d_ax);d_ax.set(yticks=range(3),yticklabels=['Full 20:8','Size 8:8','Matched 7:7'],ylim=(2.6,-.6),xlabel='Δ Fisher z',xlim=(-.005,.16),xticks=[0,.08,.16])
    heading(fig,'D','Measurement checks',101,83)
    text(fig,107,92,'Original scoring /\ndiscovery model')
    fig.legend(handles=[Line2D([],[],marker='o',c=COL['Claudin'],ls='',label='Fixed background'),Line2D([],[],marker='s',c='.45',ls='',label='Target mean')],loc='upper left',bbox_to_anchor=(105/170,1-156/170),frameon=False,borderaxespad=0,handletextpad=.4)
    panel_data('D',matching,old_panel='S13b',source_family='Ogden 2025',patient_set='14 patients',score_version='original discovery; HRC/covariates fixed',model='original state-adjusted conditioning',gene_set='full 20:8; 500 size-matched 8:8 subsets; 7 matched pairs',estimand='paired Fisher-z difference',multiplicity_family='four secondary matching tests; unchanged')
    save(fig,name)

def figure3():
    name='Figure_3';begin(name);genes=members()
    detection=read(C/'cells/compartment_detection_summary.tsv')
    studies=read(Q/'external/study_member_effects.tsv')
    og=studies.query("study=='Ogden_2025' and scoring=='fixed_background' and model=='state_enhanced'").set_index('gene').reindex(genes)
    ext=read(Q/'external/random_effects_meta.tsv').query("scoring=='fixed_background' and model=='common'").set_index('endpoint').reindex(genes)
    per=studies.query("scoring=='fixed_background' and model=='common'")
    per=per[per.study.isin(STUDIES)].copy();per['rho']=np.tanh(per.estimate)
    for frame in [og,ext]:
        for source_column,display_column in [('estimate','display_rho'),('low','display_rho_low'),('high','display_rho_high')]:frame[display_column]=np.tanh(frame[source_column])
    assert len(per)==150 and per.estimate.notna().sum()==125
    fig=canvas(170)
    a=ax(fig,23,26,22,108);b=ax(fig,51,26,27,108);c=ax(fig,88,26,34,108);d=ax(fig,132,26,33,108)
    comp=['Epithelial','Hepatocyte','Endothelial','Stromal']
    mat=detection.pivot(index='gene',columns='compartment',values='mean_detection').reindex(index=genes,columns=comp)
    im=heat(a,mat,vmin=0,vmax=1,cmap='Blues',missing='·')
    a.set(yticks=range(30),yticklabels=genes,xticks=range(4),xticklabels=['Epi.','Hep.','Endo.','Stroma'])
    a.tick_params(axis='x',labelrotation=65)
    for tick in a.get_yticklabels():
        if tick.get_text() in ['CLDN1','CLDN2','CLDN3','CLDN4','CLDN5','CLDN7','CRB3','F11R']:tick.set_fontweight('bold')
    for i,g in enumerate(genes):
        fig.patches.append(Rectangle((8/170,1-(26+(i+1)*108/30)/170),.65/170,(108/30)/170,
                                    transform=fig.transFigure,facecolor=member_color(g),edgecolor='none'))
    text(fig,3.5,26+10.5*108/30,'Claudins',rotation=90,va='center',ha='center')
    text(fig,3.5,26+25*108/30,'Polarity',rotation=90,va='center',ha='center')
    text(fig,1,26+29.5*108/30,'F11R',va='center')
    cbar(fig,im,23,155,22,'Detection',[0,.5,1])
    for axis,frame,kind in [(b,og,'n'),(c,ext,'k')]:
        for i,(gene,r) in enumerate(frame.iterrows()):
            co=member_color(gene)
            if pd.notna(r.estimate):
                e,lo,hi=np.tanh([r.estimate,r.low,r.high])
                if np.isfinite(lo) and np.isfinite(hi):axis.plot([lo,hi],[i,i],c=co,lw=.85)
                filled=axis is b or (pd.notna(r.q_member_family) and r.q_member_family<.05)
                axis.plot(e,i,'o',mec=co,mfc=co if filled else 'white',ms=3.3)
            else:axis.text(.97,i,'NA',transform=axis.get_yaxis_transform(),ha='right',va='center',c='.45')
            count=r.n_patients if kind=='n' else r.k
            axis.text(1.13 if axis is b else 1.10,i,str(int(count)) if pd.notna(count) else 'NA',transform=axis.get_yaxis_transform(),ha='center',va='center')
        axis.set(yticks=range(30),yticklabels=[],ylim=(29.5,-.5),xlim=(-.30,.40),xticks=[-.2,0,.2,.4],xlabel='Partial ρ')
        axis.tick_params(axis='y',length=0);zero(axis)
        text(fig,81.5 if axis is b else 125.4,21,kind,ha='center')
    values=per.pivot(index='gene',columns='study',values='rho').reindex(index=genes,columns=STUDIES)
    sm=heat(d,values,missing='·')
    d.set(yticks=range(30),yticklabels=[],xticks=range(5),xticklabels=['Che','Liu','Wang','Giguelay','Sathe'])
    d.tick_params(axis='x',labelrotation=65)
    cbar(fig,sm,132,155,33,'Study ρ')
    for axis in [a,b,c,d]:row_lines(axis)
    heading(fig,'A','Sources',5,5);heading(fig,'B','Ogden',49,5);heading(fig,'C','External',87,5);heading(fig,'D','Studies',131,5)
    text(fig,20.5,19,'n',ha='right')
    for j,count in enumerate([15,6,8,10]):text(fig,23+(j+.5)*22/4,19,str(count),ha='center')
    text(fig,51,13,'State-adjusted');text(fig,51,18,'Bootstrap 95% CI')
    text(fig,88,13,'Random effects / HK');text(fig,132,13,'5 original studies')
    text(fig,88,148,'External BH q,\n30 members')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',mfc='.35',mec='.35',label='< 0.05'),Line2D([],[],ls='',marker='o',mfc='white',mec='.35',label='≥ 0.05')],loc='upper left',bbox_to_anchor=(88/170,1-158/170),frameon=False,borderaxespad=0,handletextpad=.3,ncol=2,columnspacing=.5,handlelength=.7)
    panel_data('A',detection[detection.gene.isin(genes)],old_panel='Main 3a',source_family='Ogden broad compartments',patient_set='Epithelial 15; hepatocyte 6; endothelial 8; stromal 10',estimand='equal-patient mean detection fraction',ci_method='descriptive')
    panel_data('B',og.reset_index(),old_panel='Main 3b',source_family='Ogden',patient_set='14 eligible patients; gene-specific availability',model='harmonized state_enhanced',estimand='Fisher-z mean back-transformed to rho',ci_method='patient bootstrap, back-transformed')
    panel_data('C',ext.reset_index(names='gene'),old_panel='Main 3c',source_family='five external original studies',model='common / random effects',estimand='Fisher-z pooled mean back-transformed to rho',ci_method='Hartung–Knapp, back-transformed',multiplicity_family='30-member BH family')
    panel_data('D',per,source_family='five external original studies',model='common; equal-patient study means',estimand='study Fisher-z mean back-transformed to rho',ci_method='heatmap descriptive; study intervals supplied',notes='All 150 study–member combinations: 125 finite, 25 unavailable; no new heterogeneity tests.')
    save(fig,name)

def figure5():
    name='Figure_5';begin(name);defs=read(C/'definitions.tsv')
    scores=read(Q/'regional/patient_region_display_scores.tsv').query("dataset=='all_pairs' and normalization=='TMM'")
    proxies=read(Q/'regional/patient_region_proxy_scores.tsv').query("dataset=='all_pairs' and normalization=='TMM'")
    allgenes=read(Q/'regional/genome_results.tsv.gz').query("dataset=='all_pairs' and normalization=='TMM'")
    tests=read(Q/'regional/programme_tests.tsv').query("dataset=='all_pairs' and normalization=='TMM' and component=='HRC'")
    effects=read(Q/'regional/all30_member_results.tsv').query("dataset=='all_pairs' and normalization=='TMM'")
    fig=canvas(175);a=ax(fig,23,21,52,39);b=ax(fig,111,21,53,39)
    w=scores.pivot(index='patient',columns='region',values='HRC_score')
    for offset,(_,r) in zip(patient_offsets(w.index),w.iterrows()):
        a.plot(np.array([0,1])+offset,[r.micro_tumour,r.macro_tumour],c='.78',lw=.65)
        a.plot(offset,r.micro_tumour,'o',mec='.45',mfc='white',ms=3.5);a.plot(1+offset,r.macro_tumour,'s',c='.2',ms=3.3)
    a.set(xticks=[0,1],xticklabels=['Micro','Macro'],xlim=(-.35,1.35),ylabel='HRC programme score')
    a.axhline(0,c='.75',ls='--',lw=.6)
    heading(fig,'A','Paired regional HRC',5,5);text(fig,24,14,'11 patients; 88 eligible genes')
    diff=[]
    for i,comp in enumerate(['Epithelial','Liver','Endothelial']):
        p=proxies.pivot(index='patient',columns='region',values=comp);v=p.macro_tumour-p.micro_tumour
        b.scatter(i+np.linspace(-.14,.14,len(v)),v,s=15,color=['#547F8C','#B89255','#8E799B'][i])
        b.plot(i,v.mean(),'D',c='black',ms=4)
        diff.extend({'patient':p,'proxy':comp,'macro_minus_micro':value} for p,value in v.items())
    b.axhline(0,c='.65',ls='--',lw=.7)
    b.set(xticks=range(3),xticklabels=['Epithelial','Liver','Endothelial'],ylabel='Proxy change: macro − micro',xlim=(-.45,2.45))
    heading(fig,'B','Source-associated expression',92,5)
    c=ax(fig,24,112,51,51)
    hrc=set(defs.loc[defs.component.eq('HRC'),'gene']);h=allgenes[allgenes.gene.isin(hrc)]
    hc=h.pivot(index='gene',columns='model',values='logFC');assert hc.shape==(88,2) and hc.notna().all().all()
    lim=(-1.25,2.10)
    c.plot(lim,lim,c='.55',ls='--',lw=.7);c.axhline(0,c='.8',lw=.6);zero(c)
    c.scatter(hc.unadjusted,hc.source_adjusted,s=13,color=COL['HRC'],alpha=.78,edgecolors='white',linewidths=.25)
    c.set(xlim=lim,ylim=lim,xlabel='Baseline log2FC',ylabel='Source-adjusted log2FC');c.set_aspect('equal')
    heading(fig,'C','HRC gene coefficients',5,82)
    qt=tests.set_index('model').q_direction
    text(fig,50,92,'Baseline',ha='center');text(fig,88,92,'Source-adjusted',ha='right')
    for y,label,base,adjusted in [(98,'Median coefficient',f'{hc.unadjusted.median():.3f}',f'{hc.source_adjusted.median():.3f}'),(104,'Directional q',f'{qt.unadjusted:.1e}',f'{qt.source_adjusted:.3f}')]:
        text(fig,5,y,label);text(fig,50,y,base,ha='center');text(fig,77,y,adjusted,ha='center')
    d=ax(fig,112,109,52,54);genes=['CLDN3','CLDN4','CLDN7','CLDN2','CLDN5','CRB3']
    regional_forest(d,effects,genes);d.set_xlabel('Macro − micro log2FC')
    heading(fig,'D','Junction-member coefficients',92,82);model_legend(fig,96,94)
    panel_data('A',scores,old_panel='Main 5a',patient_set='11 pairs',gene_set='88 count-eligible HRC genes',normalization='TMM',model='display only; inference in frozen paired count tests',estimand='mean gene-standardized logCPM')
    panel_data('B',pd.DataFrame(diff),old_panel='Main 5b',patient_set='11 pairs',gene_set='disjoint epithelial/liver/endothelial proxies',normalization='TMM',estimand='within-patient source-proxy score change')
    panel_data('C',hc.reset_index(),patient_set='11 pairs',gene_set='88 count-eligible HRC genes',normalization='TMM',model='patient + region, with/without three source proxies',estimand='paired gene coefficients, descriptive',ci_method='not displayed; no test across genes',notes='Programme tests retained in separate panel source file; genes are not independent biological replicates.')
    tests.to_csv(OUT/'source_data/panels/Figure_5_C_programme_tests.tsv',sep='\t',index=False,encoding='utf-8')
    panel_data('D',effects[effects.gene.isin(genes)],old_panel='Main 5c',patient_set='11 pairs',normalization='TMM',model='baseline and source_adjusted',estimand='gene macro-minus-micro log2FC',ci_method='approximate QL-Wald; QL test P values in source table')
    save(fig,name)

def figure1():
    name='Figure_1';begin(name)
    maps=read(P/'Figure_1/figure1a.tsv');nc=read(P/'Figure_1/figure1c.tsv');liver=read(P/'Figure_1/figure1d.tsv')
    topo=read(P/'Figure_1/figure_phase2_topology_residual_summary_source_data.tsv');mapping=read(P/'Figure_2/figure2a.tsv')
    crop=json.loads(source('manuscript/rec_figure1_matched_frames_2026-09-10/review/crop_definition.json').read_text(encoding='utf-8'))
    fig=canvas(175);heading(fig,'A','Study framework',5,3)
    boxes=[(5,11,44,11,'Discovery: Ogden\n14 patients'),(61,11,104,11,'Independent validation\n5 original studies · 28 patients'),
           (5,25,44,11,'ISS spatial context\n17 patients'),(53,25,56,11,'HGP tissues\nBulk 15 · spatial 6'),(113,25,52,11,'Paired regions\n11 patients')]
    for x,y,w,h,label in boxes:
        fig.patches.append(Rectangle((x/170,1-(y+h)/175),w/170,h/175,transform=fig.transFigure,facecolor='#F2F4F4',edgecolor='#BCC6C9',lw=.65))
        text(fig,x+w/2,y+2.2,label,ha='center',fontsize=8)
    arrow=ax(fig,0,0,170,175);arrow.set_axis_off();arrow.annotate('',xy=(.35,.906),xytext=(.29,.906),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.8,'color':'.35'})
    framework=pd.DataFrame([{'role':'discovery','source':'Ogden','n_patients':14},{'role':'external','source':'five original studies via Atlas','n_patients':28},{'role':'bulk HGP','source':'GSE151165','n_patients':15},{'role':'spatial HGP','source':'E-MTAB-12043','n_patients':6},{'role':'regional','source':'GSE294385','n_patients':11},{'role':'same-section subset','source':'GSE294385','n_patients':10}])
    panel_data('A',framework,gene_set='NC3 / REC / HRC / REC50 roles distinct',ci_method='schematic',notes='E-MTAB-12022 and E-MTAB-12043 share four patients; no causal arrows.')
    heading(fig,'B','Original ISS coordinate fields',5,41)
    left=maps[maps.Sample.eq('ENC-P02-S03')];right=maps[maps.Sample.eq('REP-P18-S24')]
    dy=left.y.max()-left.y.min();dx=max(left.x.max()-left.x.min(),right.x.max()-right.x.min())
    # Same source-coordinate display scale and the existing central P18 window.
    palettes={'Other tissue':'#C7CDD1','Liver epithelium':'#D8AB66','Other tumour cells':'#97AFBA','NC3':'#218C84'}
    displayed=[]
    for j,(data,label,x0) in enumerate([(left,'P02 · dHGP',5),(right,'P18 · rHGP detail',86)]):
        axis=ax(fig,x0,53,65,32)
        lo,hi=(left.y.min(),left.y.max()) if j==0 else (crop['y_min'],crop['y_max'])
        mid=(data.x.min()+data.x.max())/2
        for category,co in palettes.items():
            d=data[data.display_class_new.eq(category)];axis.scatter(d.x,d.y,c=co,s=.70,linewidths=0,alpha=1,rasterized=True)
        axis.set(xlim=(mid-dx/2,mid+dx/2),ylim=(hi,lo));axis.set_aspect('equal');axis.set_axis_off()
        text(fig,x0+4,48,label)
        z=data.copy();z['in_display_window']=z.y.between(lo,hi);displayed.append(z)
        if j==1:
            overview=ax(fig,153,54,11,31);text(fig,158,49,'Overview',ha='center')
            for category,co in palettes.items():
                d=data[data.display_class_new.eq(category)];overview.scatter(d.x,d.y,c=co,s=.22,linewidths=0,alpha=1,rasterized=True)
            overview.add_patch(Rectangle((data.x.min(),lo),data.x.max()-data.x.min(),hi-lo,fill=False,edgecolor='.15',lw=.6))
            overview.invert_yaxis();overview.set_aspect('equal');overview.set_axis_off()
            from matplotlib.patches import ConnectionPatch
            for yy,edge in [(lo,0),(hi,1)]:
                fig.add_artist(ConnectionPatch(xyA=(data.x.min(),yy),coordsA=overview.transData,
                    xyB=(1,1-edge),coordsB=axis.transAxes,color='.45',lw=.5,clip_on=False))
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=co,label=lab) for lab,co in palettes.items()],loc='upper left',bbox_to_anchor=(.03,1-86/175),ncol=4,frameon=False,borderaxespad=0,handletextpad=.3,columnspacing=.8)
    panel_data('B',pd.concat(displayed),old_panel='S17a',source_family='ISS original cell coordinates',patient_set='representative P02/P18; inference uses complete fields',estimand='coordinate display',ci_method='none',notes=json.dumps(crop))
    c=ax(fig,22,101,51,23);d=ax(fig,111,101,51,23)
    groups(c,nc,'nc3_fraction_neoplastic_equal_roi','NC3 (%)',True);groups(d,liver,'any_target_neighbour_rate_equal_roi','Liver neighbour (%)',True)
    heading(fig,'C','NC3 abundance',5,93);heading(fig,'D','Tumour–liver proximity',92,93)
    panel_data('C',nc,old_panel='Main 1a',patient_set='7 dHGP / 8 rHGP',estimand='NC3 fraction, equal ROI weighting')
    panel_data('D',liver,old_panel='Main 1b',patient_set='7 dHGP / 8 rHGP',estimand='ten-nearest-neighbour liver contact frequency')
    e=ax(fig,32,143,35,19);rows=[];labels=[]
    for metric,label in [('homotypic_any_neighbour','NC3–NC3'),('dhc_dual_interface','NC3 + DHC')]:
        for h,display in [('EHGP','dHGP'),('RHGP','rHGP')]:
            r=topo[(topo.metric==metric)&(topo.hgp==h)].iloc[0]
            rows.append({'effect':100*r.effect,'low':100*r.bootstrap_ci_lower,'high':100*r.bootstrap_ci_upper});labels.append(display)
    forest(e,pd.DataFrame(rows),labels,color=[COL['dHGP'],COL['rHGP']]*2)
    for i,line in enumerate([x for x in e.lines if x.get_marker()=='o']):line.set_marker('o' if i%2==0 else 's')
    e.set_xlabel('Conditional residual (pp)',labelpad=1)
    text(fig,3,147.75,'NC3–NC3',va='center');text(fig,1,157.25,'NC3 + DHC',va='center')
    heading(fig,'E','Conditional adjacency',5,133)
    panel_data('E',topo[topo.metric.isin(['homotypic_any_neighbour','dhc_dual_interface'])],old_panel='Main 1c',model='conditional null, k=10',estimand='observed minus conditional expectation',ci_method='patient bootstrap')
    f=ax(fig,80,143,84,19)
    order=['Goblet','Hypoxia','iREC','REC','UPR','Colonocyte','TA1','Stem','Intermediate','Stem NOTUM']
    assert set(order)==set(mapping.Ogden_state)
    # Retain the fourth source state's original NCγ label.
    states=['NC1','NC2','NC3']+[x for x in mapping.ISS_state.unique() if x not in ['NC1','NC2','NC3']]
    matrix=mapping.pivot(index='ISS_state',columns='Ogden_state',values='pearson_correlation').reindex(index=states,columns=order)
    short={'Colonocyte':'Colono.','Intermediate':'Interm.','Stem NOTUM':'Stem-N.'}
    im=heat(f,matrix);f.set(yticks=range(4),yticklabels=states,xticks=range(10),xticklabels=[short.get(v,v) for v in order])
    f.tick_params(axis='x',labelrotation=45,pad=1);f.tick_params(axis='y',pad=1)
    f.add_patch(Rectangle((order.index('REC')-.5,states.index('NC3')-.5),1,1,fill=False,edgecolor='black',lw=1.2))
    f.text(order.index('REC'),states.index('NC3'),f'{matrix.loc["NC3","REC"]:.3f}',ha='center',va='center',fontsize=8,color='black')
    heading(fig,'F','Reference mapping',77,133)
    # Compact labelled scale above the matrix, independent of any detection scale.
    cb=fig.colorbar(im,cax=ax(fig,138,138,26,1.5),orientation='horizontal',ticks=[-.5,0,.5]);cb.ax.xaxis.set_ticks_position('top');cb.ax.tick_params(labelsize=8,length=1,pad=0)
    text(fig,116,138,'Pearson r',fontsize=8)
    panel_data('F',mapping,old_panel='S9a',estimand='centroid Pearson correlation across 94 shared genes',ci_method='descriptive; full mapping sensitivity in S1',notes='Fixed source/reference ordering; no clustering or probability interpretation.')
    save(fig,name)

def figure4():
    name='Figure_4';begin(name)
    spots=read(B/'spatial/spot_scores.tsv.gz');meta=read(B/'spatial/sample_qc.tsv')
    images=read(P/'supplementary_figure_3/representative_spatial_summary.tsv')
    sc=read(B/'spatial/patient_correlations.tsv');matched=read(B/'integration/matched_gene_patient_HGP_scores.tsv')
    bulk=read(N/'bulk/patient_scores.tsv');effects=read(B/'spatial/region_effects.tsv')
    fig=canvas(175);heading(fig,'A','H&E and expression-defined regions',5,3)
    samples=[]
    for j,patient in enumerate(['PT55','PT54']):
        r=images[images.patient.eq(patient)].iloc[0];s=meta[meta.patient.eq(patient)].iloc[0]
        im=plt.imread(source(P/'supplementary_figure_3'/Path(r.image_path).name))
        d=spots[spots['sample'].eq(s['sample'])].copy()
        original=read(P/'supplementary_figure_3/representative_spatial_spots.tsv.gz')
        o=original[original.patient.eq(patient)].copy()
        # Use the archived, metadata-scaled coordinates joined by barcode.
        d=d.merge(o[['barcode','x_hires','y_hires']],on='barcode',how='left',validate='one_to_one')
        assert d[['x_hires','y_hires']].notna().all().all()
        category=np.full(len(d),'Other retained spots',dtype=object)
        category[d.liver_side]='Liver-like'
        category[d.tumour_side&d.distance.between(1,5)]='Near tumour-side'
        category[d.tumour_side&(d.distance>5)&np.isfinite(d.distance)]='Deep tumour-side'
        d['display_region']=category;samples.append(d)
        for k in range(2):
            x=5+j*84+k*39;axis=ax(fig,x,23,37,37);axis.imshow(im,origin='upper')
            if k:
                for label,co in REGION.items():
                    v=d[d.display_region.eq(label)];axis.scatter(v.x_hires,v.y_hires,s=1.1,c=co,linewidths=0,alpha=.95,rasterized=True)
            length=500*float(r.pixels_per_micrometre);x0=im.shape[1]*.06;y0=im.shape[0]*.93
            axis.plot([x0,x0+length],[y0,y0],c='black',lw=1.4)
            axis.text(x0,im.shape[0]*.84,'500 μm',fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':.5})
            axis.set(xlim=(0,im.shape[1]),ylim=(im.shape[0],0));axis.set_axis_off()
            text(fig,x+18.5,17,'H&E' if k==0 else 'Expression-defined\nregions',ha='center',va='center')
        text(fig,5+j*84,10,f'{patient} · {r.hgp}')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=co,label=lab) for lab,co in REGION.items()],loc='upper left',bbox_to_anchor=(.03,1-62/175),ncol=2,frameon=False,borderaxespad=0,handletextpad=.4,columnspacing=1)
    panel_data('A',pd.concat(samples),old_panel='S1e/f with current frozen masks',source_family='E-MTAB-12043',patient_set='PT55/PT54; all six shown in S4',estimand='original H&E plus expression-defined masks',ci_method='none')
    b=ax(fig,22,89,52,27);c=ax(fig,111,89,53,27)
    w=sc.query("restriction=='all_tumour' and model=='adjusted'");w=w[w.component.isin(['Claudin','Polarity'])]
    bw=w.pivot(index=['patient','hgp'],columns='component',values='rho')
    for offset,((patient,hgp),row) in zip(patient_offsets(bw.index.get_level_values('patient'),.10),bw.iterrows()):
        b.plot(np.array([0,1])+offset,row[['Claudin','Polarity']],c=COL[hgp],alpha=.45,lw=.7)
        b.plot(np.array([0,1])+offset,row[['Claudin','Polarity']],ls='',marker='o' if hgp=='dHGP' else 's',c=COL[hgp],ms=3.3)
    b.axhline(0,c='.7',ls='--',lw=.6);b.set(xticks=[0,1],xticklabels=['Claudin','Polarity'],xlim=(-.25,1.25),ylabel='Spatial partial ρ')
    heading(fig,'B','Within-tissue co-expression',5,78)
    near=matched.query("dataset=='Spatial_rHGP_minus_dHGP_near' and normalization=='TMM' and component=='Claudin' and score=='mean_logCPM'")
    groups(c,near,'value','Mean log2CPM (13 genes)')
    for x,hgp in enumerate(['dHGP','rHGP']):
        block=near[near.hgp.eq(hgp)]
        for offset,r in zip(np.linspace(-.12,.12,len(block)),block.itertuples()):
            move={'PT55':(0,13),'PT61':(11,3),'PT68':(8,-4),'PT36':(-13,2),'PT44':(8,2),'PT54':(7,-2)}[r.patient]
            c.annotate(r.patient,(x+offset,r.value),xytext=move,textcoords='offset points',fontsize=8,
                       ha='right' if move[0]<0 else 'center' if move[0]==0 else 'left',va='center',arrowprops={'arrowstyle':'-','color':'.65','lw':.4})
    b.legend(handles=hgp_handles(),loc='lower left',bbox_to_anchor=(0,1.02),ncol=2,frameon=False,
             borderaxespad=0,handlelength=.8,handletextpad=.3,columnspacing=.8)
    heading(fig,'C','Near-region Claudin · TMM',92,78);text(fig,115,83,'13/21 genes; exact P = 0.10')
    panel_data('B',w,old_panel='Main 4a',patient_set='six patients, three per HGP',model='technical and source-score adjusted',estimand='spatial tumour-side HRC partial rho')
    panel_data('C',near,old_panel='Main 4b',patient_set='three dHGP / three rHGP',gene_set='13 count-eligible claudins',normalization='TMM',estimand='mean log2CPM',ci_method='group points/means; exact patient allocation P=0.10')
    d=ax(fig,22,140,52,20);e=ax(fig,115,140,49,20)
    for j,comp in enumerate(['HRC','Junction']):
        for i,hgp in enumerate(['dHGP','rHGP']):
            subset=bulk[bulk.hgp.eq(hgp)];v=subset[comp];x=j*3+i
            d.scatter(x+patient_offsets(subset.sample_id,.18),v,c=COL[hgp],s=11,marker='o' if hgp=='dHGP' else 's')
            d.plot([x-.2,x+.2],[v.mean()]*2,c='.25',lw=1)
    d.axhline(0,c='.7',ls='--',lw=.6);d.set(xticks=[.5,3.5],xticklabels=['HRC (97)\n9 dHGP / 6 rHGP','Junction (28)\n9 dHGP / 6 rHGP'],ylabel='Gene-z score')
    heading(fig,'D','Independent bulk · 15 patients',5,133)
    loc=effects.query("hops==5 and restriction=='all_tumour' and score=='mean_logCPM' and endpoint=='HGP_localization'")
    loc=loc.set_index('component').loc[['HRC','Claudin','Polarity','Junction']].reset_index()
    forest(e,loc,loc.component.tolist());e.set_xlabel('rHGP − dHGP in near − deep\n(mean log₂CPM)',labelpad=1)
    heading(fig,'E','Localization contrast',92,133)
    panel_data('D',bulk,old_panel='Main 4d',source_family='GSE151165',patient_set='9 dHGP / 6 rHGP',gene_set='full HRC 97; junction 28',normalization='frozen full-programme bulk projection',estimand='full-programme gene-z score')
    panel_data('E',loc,old_panel='Main 4f',source_family='E-MTAB-12043',patient_set='six paired near/deep profiles, three per HGP',normalization='frozen library-scaled logCPM, paired bands',model='five-hop all tumour-side localization contrast',estimand='rHGP minus dHGP in patient near-minus-deep expression',ci_method='patient bootstrap; exact allocation tests')
    save(fig,name)

if __name__=='__main__':
    import sys
    choices=sys.argv[1:] or ['2','3','5','1','4']
    for number in choices:globals()['figure'+number]()
