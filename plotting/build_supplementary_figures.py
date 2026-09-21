"""Six single-page supplements drawn from archived numerical outputs."""
from figure_common import *
from build_main_figures import STUDIES

def supplementary1():
    name='Supplementary_Figure_1';begin(name)
    composition=read(P/'Figure_1/figure1b.tsv');adj=read(P/'supplementary_figure_1/figureS1a.tsv')
    loo=read('analysis_results/cross_platform_state_anchor/leave_one_patient_out_mapping.tsv')
    topo=read('analysis_results/deep_biology_upgrade_2026-08-31/phase2_mechanistic_specificity/iss_cohesive_dual_interface/iss_cohesive_dual_interface_group_summary.tsv')
    fig=canvas(165);a=ax(fig,20,26,59,33);b=ax(fig,111,19,52,40)
    w=composition.pivot(index=['hgp','canonical_patient'],columns='state',values='fraction').fillna(0).reset_index().sort_values(['hgp','NC3'])
    states=['NC1','NC2','NC3']+[s for s in composition.state.unique() if s not in ['NC1','NC2','NC3']]
    bottom=np.zeros(len(w));positions=np.arange(len(w),dtype=float);positions[7:]+=.8
    for state,co in zip(states,['#607D8B','#A3B8C2',COL['Claudin'],'#BEABD0']):
        a.bar(positions,w[state]*100,bottom=bottom,color=co,width=.86,label=state);bottom+=w[state]*100
    a.set(xticks=positions,xticklabels=w.canonical_patient,ylim=(0,100),ylabel='Neoplastic cells (%)');a.tick_params(axis='x',labelrotation=90)
    a.text(positions[:7].mean(),1.04,'dHGP',transform=a.get_xaxis_transform(),ha='center',va='bottom')
    a.text(positions[7:].mean(),1.04,'rHGP',transform=a.get_xaxis_transform(),ha='center',va='bottom')
    a.legend(loc='upper center',bbox_to_anchor=(.5,1.40),ncol=4,frameon=False,handlelength=.8,columnspacing=.65)
    heading(fig,'A','Patient state composition',5,4)
    d=adj.rename(columns={'mean_difference_RHGP_minus_EHGP':'effect','mean_difference_bootstrap_CI_low':'low','mean_difference_bootstrap_CI_high':'high'}).copy()
    d[['effect','low','high']]*=100
    forest(b,d,d.focal_cluster.tolist());b.set_xlabel('rHGP − dHGP (pp)')
    heading(fig,'B','Liver-adjacency differences',92,4)
    panel_data('A',composition,old_panel='S17b',patient_set='15 pure-HGP patients',estimand='state proportions with equal ROI weighting',ci_method='descriptive')
    panel_data('B',adj,old_panel='S1a',patient_set='7 dHGP / 8 rHGP',estimand='state-specific liver-adjacency difference',ci_method='patient bootstrap')
    c=ax(fig,20,91,59,51)
    for i,dataset in enumerate(['ISS','Ogden']):
        ds=loo[loo.omitted_dataset.str.contains(dataset,case=False)].copy()
        valid=ds.top_Ogden_state.eq('REC')
        v=ds.loc[valid,'top_minus_second_margin']
        c.scatter(i+np.linspace(-.16,.16,len(v)),v,c=COL['Claudin'] if i==0 else COL['Polarity'],s=15)
    c.axhline(0,c='.65',ls='--',lw=.7);c.set(xticks=[0,1],xticklabels=['ISS omissions','Ogden omissions'],xlim=(-.4,1.4),ylabel='REC − next-reference Pearson r')
    heading(fig,'C','Mapping after patient omission',5,78)
    text(fig,5,151,'17 ISS / 10 evaluable Ogden omissions;\n5 further Ogden omissions changed coverage.')
    panel_data('C',loo,old_panel='S9b',estimand='NC3 REC-leading Pearson similarity margin',ci_method='omission sensitivity',notes='Coverage-changing omissions retained in source; displayed evaluable REC-leading omissions.')
    d_ax=ax(fig,115,91,48,51)
    t=topo[(topo.state=='NC3')&(topo.value_type=='spatial_null_residual')&(topo.test=='within_hgp_mean_vs_zero')&topo.metric.isin(['homotypic_any_neighbour','dhc_dual_interface'])]
    labels=[]
    for i,(metric,h) in enumerate([(m,h) for m in ['homotypic_any_neighbour','dhc_dual_interface'] for h in ['EHGP','RHGP']]):
        labels.append(('NC3–NC3' if metric.startswith('homotypic') else 'NC3 + DHC')+(' / dHGP' if h=='EHGP' else ' / rHGP'))
        for k,offset,style in [(5,-.23,'-'),(10,0,'--'),(20,.23,':')]:
            row=t[(t.metric==metric)&(t.hgp==h)&(t.k==k)].iloc[0]
            co=COL['dHGP' if h=='EHGP' else 'rHGP']
            d_ax.plot(100*np.array([row.bootstrap_ci_lower,row.bootstrap_ci_upper]),[i+offset]*2,c=co,lw=1,ls=style)
            d_ax.plot(100*row.effect,i+offset,'o' if h=='EHGP' else 's',c=co,ms=3.4)
    d_ax.set(yticks=range(4),yticklabels=labels,ylim=(3.6,-.6),xlabel='Conditional residual (pp)');zero(d_ax)
    heading(fig,'D','Neighbourhood scales',92,78)
    fig.legend(handles=[Line2D([],[],ls=style,color='.3',label=f'k = {k}') for k,style in [(5,'-'),(10,'--'),(20,':')]],loc='upper left',bbox_to_anchor=(.58,.09),frameon=False,ncol=3,handlelength=1.4,handletextpad=.3,columnspacing=.5)
    panel_data('D',t,source_family='ISS',patient_set='7 dHGP / 8 rHGP',model='conditional null at k=5/10/20',estimand='observed minus conditional expectation',ci_method='frozen patient bootstrap')
    save(fig,name)

def supplementary2():
    name='Supplementary_Figure_2';begin(name)
    g=read('manuscript/rec_reviewer_revision_2026-09-07/source_data/Figure_3/table_phase2_direct_program_gsea.tsv')
    summary=read('analysis_results/rec_public_upgrade_2026-09-07/within_cell/conditional_summary.tsv')
    pats=read('analysis_results/rec_public_upgrade_2026-09-07/within_cell/patient_conditional_correlations.tsv')
    paired=read('analysis_results/rec_hgp_focused_revision_2026-09-07/paired_states/paired_state_contrasts.tsv')
    subsets=read(N/'cells/subset_mean_differences.tsv');features=read(N/'cells/matching_features.tsv');pairs=read(N/'cells/measurement_matching_pairs.tsv')
    component=read(C/'cells/component_summary.tsv')
    fig=canvas(175);a=ax(fig,44,18,30,37);b=ax(fig,123,18,40,37)
    programs=['Core HRC','Regenerative epithelium','Partial EMT','E2F targets','G2M checkpoint','Homotypic cell–cell adhesion','Tight-junction interactions','Keratinization'];states=['Hypoxia','UPR','iREC']
    mat=g.pivot(index='program_label',columns='comparator_state',values='nes').loc[programs,states]
    q=g.pivot(index='program_label',columns='comparator_state',values='fdr_q').loc[programs,states]
    im=heat(a,mat);a.set(yticks=range(8),yticklabels=['HRC','Regenerative','Partial EMT','E2F','G2M','Homotypic adhesion','Tight junction','Keratinization'],xticks=range(3),xticklabels=['Hyp.','UPR','iREC'])
    for i in range(8):
        for j in range(3):
            if q.iloc[i,j]<.05:a.plot(j+.27,i-.27,'o',ms=1.5,mec='none',mfc='white' if abs(mat.iloc[i,j])>2 else '.15')
    cbar(fig,im,44,63,30,'NES')
    heading(fig,'A','Paired state enrichment',5,4)
    selected=summary.query("method=='matched' and state=='REC' and minimum_cells==20")
    selected=selected.set_index(['rival','adjustment']).loc[[(r,a) for r in ['WHITE_RSC','ICMS3_MINUS_ICMS2'] for a in ['technical_only','technical_plus_rival']]].reset_index()
    forest(b,selected,['Technical','+ RSC','Technical','+ iCMS'],estimate='rho',low='CI_low',high='CI_high')
    for y,label in [(.5,'RSC'),(2.5,'iCMS')]:b.text(-.68,y,label,transform=b.get_yaxis_transform(),ha='right',va='center',fontweight='bold')
    b.axhline(1.5,color='.84',lw=.5)
    text(fig,5,63,'• FDR < 0.05\nwithin each GSEA run')
    for i,r in selected.iterrows():
        v=pats[(pats.rival==r.rival)&(pats.adjustment==r.adjustment)&pats.method.eq('matched')&pats.state.eq('REC')&pats.n_cells.ge(20)]
        b.scatter(v.rho,i+np.linspace(-.15,.15,len(v)),s=10,facecolors='none',edgecolors='.65',linewidths=.55)
    b.set_xlabel('HRC–junction partial ρ');heading(fig,'B','Within-REC conditioning',92,4)
    panel_data('A',g[g.program_label.isin(programs)],old_panel='S12a',patient_set='paired REC vs Hypoxia/UPR/iREC',estimand='GSEA NES',multiplicity_family='original within-run GSEA FDR')
    panel_data('B',selected,old_panel='S12d',patient_set='nine REC patients',score_version='original disjoint matched scores',estimand='conditional partial correlation',ci_method='patient bootstrap')
    pats.to_csv(OUT/'source_data/panels/Supplementary_Figure_2_B_patient_values.tsv',sep='\t',index=False,encoding='utf-8')
    c=ax(fig,37,89,38,30);d=ax(fig,111,89,52,30)
    pp=paired[paired.analysis_role.eq('primary_extension')].copy()
    labels=[('RSC' if r.rival=='WHITE_RSC' else 'iCMS')+' / '+r.comparator.replace('Hypoxia','Hyp.')+f' ({r.n_pairs})' for r in pp.itertuples()]
    forest(c,pp,labels,estimate='fisher_z_difference',low='CI_low',high='CI_high');c.set_xlabel('REC − state (Δ Fisher z)')
    heading(fig,'C','Between-state differences',5,78)
    for mode,co,label in [('fixed_background',COL['Claudin'],'Fixed background'),('target_mean','.45','Target mean')]:
        v=subsets.loc[subsets.scoring.eq(mode),'difference_z'];d.hist(v,bins=np.linspace(subsets.difference_z.min(),subsets.difference_z.max(),29),histtype='step',color=co,lw=1.1,linestyle='-' if mode=='fixed_background' else '--',label=label)
    zero(d);d.set(xlabel='Mean-patient Δ Fisher z',ylabel='Subsets');d.legend(frameon=False,fontsize=8,loc='upper left',handlelength=1)
    heading(fig,'D','500 membership subsets',92,78)
    panel_data('C',pp,old_panel='S2a',patient_set='nine/six/five paired patients',estimand='paired REC-minus-state Fisher-z difference',multiplicity_family='six paired comparisons')
    panel_data('D',subsets,old_panel='S13c',patient_set='14 patients',score_version='original discovery',estimand='subset mean-patient Fisher-z difference',ci_method='membership perturbation distribution, not patient CI')
    e=ax(fig,23,143,51,20);f=ax(fig,124,143,39,20)
    ff=features.set_index('gene')
    for r in pairs.itertuples():
        v=ff.loc[[r.polarity_gene,r.claudin_gene]];e.plot(v.z_log_mean,v.z_logit_detection,c='.7',lw=.6)
    for k,mark in [('Claudin','o'),('Polarity','s')]:
        for matched in [False,True]:
            v=features[(features.component==k)&features.matched.eq(matched)];e.scatter(v.z_log_mean,v.z_logit_detection,s=15,marker=mark,facecolors=COL[k] if matched else 'none',edgecolors=COL[k],linewidths=.7)
    e.set(xlabel='Standardized log expression',ylabel='Logit detection (z)');e.xaxis.set_major_locator(MaxNLocator(3));e.yaxis.set_major_locator(MaxNLocator(3))
    heading(fig,'E','Seven matched gene pairs',5,135)
    rr=[];labels=[]
    for scope,scoring in [('All_epithelial','matched'),('All_epithelial','raw'),('REC','matched'),('REC','raw')]:
        r=component[(component.scope==scope)&(component.scoring==scoring)&component.model.eq('conditional')&component.component.str.startswith('Polarity_minus')].iloc[0].copy()
        r['effect'],r['low'],r['high']=-r.effect,-r.high,-r.low;rr.append(r)
        labels.append(('All epithelial' if scope=='All_epithelial' else 'REC')+' / '+('matched' if scoring=='matched' else 'all raw'))
    result=pd.DataFrame(rr);forest(f,result,labels);f.set_xlabel('Claudin − Polarity (Δz)');heading(fig,'F','Original scoring sensitivity',92,135)
    panel_data('E',features,old_panel='S13d',source_family='Ogden',estimand='expression/detection matching features',ci_method='none; fixed seven-pair selection')
    pairs.to_csv(OUT/'source_data/panels/Supplementary_Figure_2_E_pairs.tsv',sep='\t',index=False,encoding='utf-8')
    panel_data('F',result,old_panel='S10b',score_version='original matched versus all-raw scoring',estimand='Claudin minus Polarity Fisher-z difference',notes='All-raw scoring changes HRC and covariates; differs from Figure 2D target-only sensitivity.')
    save(fig,name)

def supplementary3():
    name='Supplementary_Figure_3';begin(name);genes=members()
    p=read(Q/'external/patient_component_effects.tsv').query("scoring=='fixed_background' and model=='common'");p=p[p.study.isin(STUDIES)].copy()
    loo=read(Q/'external/leave_one_study.tsv').query("scoring=='fixed_background' and endpoint=='Claudin_minus_Polarity'")
    m=read(Q/'external/random_effects_meta.tsv');comp=m[m.endpoint.eq('Claudin_minus_Polarity')].copy();member=m.query("scoring=='fixed_background' and model=='common'").set_index('endpoint').reindex(genes)
    fig=canvas(175);a=ax(fig,22,20,54,94);offset=0;positions=[];labels=[];p=p.set_index('study').loc[STUDIES].reset_index()
    for s in STUDIES:
        g=p[p.study.eq(s)]
        for r in g.itertuples():
            a.plot([0,r.delta_z],[offset]*2,c='.82',lw=.6);a.plot(r.delta_z,offset,'o',c=COL['Claudin'],ms=3.4)
            positions.append(offset);labels.append(str(r.patient).split('.')[-1]);offset+=1
        a.axhline(offset-.5,c='.83',lw=.6);text(fig,76,20+(offset-len(g)/2)/33*94,s.split('_')[0],rotation=90,ha='left',va='center');offset+=1
    a.set(yticks=positions,yticklabels=labels,ylim=(offset-1,-1),xlabel='Patient Δ Fisher z');zero(a)
    heading(fig,'A','All external patient differences',5,5)
    b=ax(fig,120,20,43,34);forest(b,loo,['Omit '+x.split('_')[0] for x in loo.excluded_study],estimate='estimate');b.set_xlabel('Pooled Δ Fisher z')
    heading(fig,'B','Study omissions',92,5)
    c=ax(fig,120,88,43,22);forest(c,comp,['Fixed\nbackground','All target\nmeans'],estimate='estimate');c.set_xlabel('Pooled Δ Fisher z')
    limits=pd.concat([loo,comp])[['low','high']]
    lo,hi=min(0,limits.low.min()),limits.high.max();pad=(hi-lo)*.06
    for axis in [b,c]:
        axis.set_xlim(min(-.05,lo-pad),max(.20,hi+pad))
        axis.set_xticks([-.05,0,.10,.20]);axis.tick_params(axis='x',pad=2)
    heading(fig,'C','External scoring sensitivity',92,74)
    d=ax(fig,18,143,146,6)
    data=np.array([member.I2.to_numpy()]);im=heat(d,data,vmin=0,vmax=100,cmap='Blues',missing='NA')
    d.set(yticks=[0],yticklabels=['I²'],xticks=[])
    counts=ax(fig,18,150,146,4)
    counts.set(xlim=(-.5,29.5),ylim=(-.5,.5),yticks=[0],yticklabels=['k'],xticks=range(30),xticklabels=genes)
    counts.tick_params(length=0);counts.tick_params(axis='x',labelrotation=90,pad=1)
    for spine in counts.spines.values():spine.set_visible(False)
    for i,r in enumerate(member.itertuples()):counts.text(i,0,str(int(r.k)) if pd.notna(r.k) else 'NA',ha='center',va='center',color='.15')
    heading(fig,'D','Member heterogeneity',5,130)
    cb=fig.colorbar(im,cax=ax(fig,122,134,40,2.4),orientation='horizontal',ticks=[0,50,100]);cb.ax.xaxis.set_ticks_position('top');cb.ax.tick_params(labelsize=8,length=2,pad=1)
    text(fig,105,134,'I² (%)')
    panel_data('A',p,source_family='five external original studies',patient_set='all 28 eligible patients',model='harmonized common',estimand='patient Claudin-minus-Polarity Fisher-z difference',ci_method='individual patient estimates')
    panel_data('B',loo,old_panel='Main 2d',model='random effects with study omission',estimand='pooled component difference',ci_method='Hartung–Knapp')
    panel_data('C',comp,score_version='fixed background versus all target means',model='common / random effects',estimand='pooled component difference',ci_method='Hartung–Knapp')
    panel_data('D',member.reset_index(names='gene'),gene_set='all 30 junction members',estimand='I-squared with evaluable study count',ci_method='descriptive; no heterogeneity-derived classification')
    save(fig,name)

def supplementary4():
    name='Supplementary_Figure_4';begin(name)
    spots=read(B/'spatial/spot_scores.tsv.gz');meta=read(B/'spatial/sample_qc.tsv')
    fig=canvas(165);order=['PT55','PT61','PT68','PT36','PT44','PT54'];counts=[]
    for i,patient in enumerate(order):
        row=meta[meta.patient.eq(patient)].iloc[0];d=spots[spots['sample'].eq(row['sample'])].copy()
        region=np.full(len(d),'Other retained spots',dtype=object);region[d.liver_side]='Liver-like';region[d.tumour_side&d.distance.between(1,5)]='Near tumour-side';region[d.tumour_side&(d.distance>5)&np.isfinite(d.distance)]='Deep tumour-side'
        d['display_region']=region
        x=5+(i%3)*55;top=16+(i//3)*75;axis=ax(fig,x,top,50,50)
        for label,co in REGION.items():
            v=d[d.display_region.eq(label)];axis.scatter(v.pixel_col,v.pixel_row,s=1.5,c=co,linewidths=0,alpha=1,rasterized=True)
            counts.append({'patient':patient,'hgp':row.hgp,'region':label,'spots':len(v)})
        axis.set_aspect('equal');axis.invert_yaxis();axis.set_axis_off();heading(fig,chr(65+i),f'{patient} · {row.hgp}',x,top-10)
        count=d.display_region.value_counts()
        text(fig,x+25,top+51,f'{len(d):,} spots\nNear {count.get("Near tumour-side",0):,} · Deep {count.get("Deep tumour-side",0):,}',ha='center')
        panel_data(chr(65+i),d,old_panel='S8 '+patient,source_family='E-MTAB-12043',patient_set=patient,model='frozen epithelial-minus-liver expression mask, five hops',estimand='spot coordinates and region labels',ci_method='descriptive',notes='Source orientation and aspect preserved; no mask smoothing or manual pathology boundary.')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=co,label=label) for label,co in REGION.items()],loc='upper center',bbox_to_anchor=(.5,1-153/165),frameon=False,ncol=2,handletextpad=.4,columnspacing=1.5,borderaxespad=0)
    pd.DataFrame(counts).to_csv(OUT/'source_data/panels/Supplementary_Figure_4_region_counts.tsv',sep='\t',index=False,encoding='utf-8')
    save(fig,name)

def supplementary5():
    name='Supplementary_Figure_5';begin(name);genes=members()
    members0=read(N/'integration/all_junction_members.tsv').set_index('gene').reindex(genes)
    coverage=read(B/'integration/coverage_only_claudin_HGP.tsv');bulk=read(N/'bulk/score_effects.tsv');models=read(N/'bulk/composition_models.tsv')
    scores=read(B/'spatial/region_scores.tsv');diag=read(N/'diagnostics/patient_normalization_diagnostics.tsv');status=read(N/'integration/member_measurement_status.tsv')
    fig=canvas(175);a=ax(fig,20,20,28,94)
    columns=['Bulk_rHGP_minus_dHGP_logFC','Spatial_rHGP_minus_dHGP_tumour_logFC','Spatial_rHGP_minus_dHGP_near_logFC']
    mat=members0[columns];im=heat(a,mat,missing='');a.set(yticks=range(30),yticklabels=genes,xticks=range(3),xticklabels=['Bulk','Tumour','Near']);a.tick_params(axis='x',labelrotation=55);row_lines(a)
    reasons=status.pivot(index='gene',columns='column',values='status')
    for i,gene in enumerate(genes):
        for j,column in enumerate(columns):
            if pd.isna(mat.loc[gene,column]):
                glyph={'not_measured':'×','abundance_ineligible':'/','insufficient_detection_or_variation':'·'}[reasons.loc[gene,column]]
                a.text(j,i,glyph,ha='center',va='center',c='.4',fontsize=8)
    cb=fig.colorbar(im,cax=ax(fig,53,31,2.3,67),orientation='vertical')
    cb.ax.tick_params(labelsize=8,length=2,pad=1);cb.set_label('TMM log2FC',fontsize=8,labelpad=3)
    heading(fig,'A','HGP members',5,5);text(fig,20,14,'rHGP − dHGP')
    b=ax(fig,108,18,55,19)
    d=coverage[coverage.dataset.eq('Spatial_rHGP_minus_dHGP_near')];rows=[];labels=[]
    for membership,norm,label in [('all_assayed','TMM','21 / TMM'),('all_assayed','library_only','21 / library'),('count_eligible','TMM','13 / TMM'),('count_eligible','library_only','13 / library')]:
        rows.append(d[(d.membership==membership)&(d.normalization==norm)].iloc[0]);labels.append(label)
    forest(b,pd.DataFrame(rows),labels);b.set_xlabel('Mean log2CPM difference')
    heading(fig,'B','Coverage and normalization',68,5)
    c=ax(fig,113,60,50,24);rows=[];labels=[]
    for k in ['HRC','Junction']:
        for version,label in [(k,'full'),(k+'_disjoint','disjoint')]:rows.append(bulk[bulk.component.eq(version)].iloc[0]);labels.append(k+' / '+label)
        for cv in ['Epithelial','Liver']:rows.append(models[(models.component==k)&(models.covariate==cv)].iloc[0]);labels.append(k+' / +'+('epi.' if cv=='Epithelial' else 'liver'))
    bb=pd.DataFrame(rows);forest(c,bb,labels);c.set_xlabel('Programme-score coefficient')
    heading(fig,'C','Bulk definition and adjustment',68,49)
    d_ax=ax(fig,104,107,59,16);ss=scores.query("hops==5 and restriction=='all_tumour' and score=='mean_logCPM' and basis=='paired_bands'");records=[]
    for i,comp in enumerate(['Claudin','Polarity','HRC']):
        wide=ss[ss.component.eq(comp)].pivot(index=['patient','hgp'],columns='region',values='value')
        for j,((p,h),r) in enumerate(wide.iterrows()):
            value=r.near-r.deep;d_ax.scatter(i+(j-2.5)*.045,value,c=COL[h],s=12,marker='o' if h=='dHGP' else 's');records.append({'patient':p,'hgp':h,'component':comp,'near_minus_deep':value})
    d_ax.axhline(0,c='.7',ls='--',lw=.6);d_ax.set(xticks=range(3),xticklabels=['Claudin','Polarity','HRC'],ylabel='Near − deep\n(mean log2CPM)')
    heading(fig,'D','Patient localization differences',68,97)
    contexts=[('Bulk_rHGP_minus_dHGP','Bulk'),('Spatial_rHGP_minus_dHGP_tumour','Tumour'),('Spatial_rHGP_minus_dHGP_near','Near')]
    heading(fig,'E','RNA composition and normalization',5,132)
    for i,(ctx,label) in enumerate(contexts):
        axis=ax(fig,22+53*i,141,32,25);data=diag[diag.dataset.eq(ctx)]
        for h in ['dHGP','rHGP']:
            v=data[data.hgp.eq(h)];axis.scatter(v.top10_count_fraction*100,v.log2_norm_factor,c=COL[h],s=13,marker='o' if h=='dHGP' else 's')
        axis.axhline(0,c='.7',lw=.6);axis.xaxis.set_major_locator(MaxNLocator(3));axis.yaxis.set_major_locator(MaxNLocator(2))
        text(fig,54+53*i,138,label,ha='right')
        if i==0:axis.set_ylabel('log2 factor')
    text(fig,87,172,'Top-ten transcript fraction (%)',ha='center')
    panel_data('A',members0[columns].reset_index(),old_panel='S14c first three columns',normalization='TMM',estimand='gene rHGP-minus-dHGP log2FC',notes='Full range shown; count-ineligible entries distinguished from valid zero.')
    status.to_csv(OUT/'source_data/panels/Supplementary_Figure_5_A_measurement_status.tsv',sep='\t',index=False,encoding='utf-8')
    panel_data('B',d,old_panel='Main 4c',normalization='TMM/library only',gene_set='21 assayed or 13 eligible claudins',estimand='near-region mean-logCPM HGP difference',ci_method='within-HGP patient bootstrap')
    panel_data('C',bb,old_panel='S3b/S7d',source_family='GSE151165',estimand='full/disjoint scores and separate source-adjusted HGP coefficients',ci_method='patient bootstrap')
    panel_data('D',pd.DataFrame(records),old_panel='Main 4e',estimand='patient near-minus-deep mean logCPM',normalization='frozen paired-band scores')
    panel_data('E',diag,old_panel='S11b/d/f',estimand='top-ten transcript fraction and TMM factor',ci_method='descriptive; no normalization selection')
    save(fig,name)

def supplementary6():
    name='Supplementary_Figure_6';begin(name);genes=members()
    d=read(Q/'regional/all30_member_results.tsv');tests=read(Q/'regional/programme_tests.tsv');vif=read(Q/'regional/within_patient_VIF.tsv');diag=read(Q/'regional/model_diagnostics.tsv')
    fig=canvas(175);a=ax(fig,20,22,43,88);b1=ax(fig,88,22,30,88);b2=ax(fig,134,22,30,88)
    specs=[(dataset,norm,model) for dataset in ['all_pairs','same_section'] for norm in ['TMM','library_only'] for model in ['unadjusted','source_adjusted']]
    matrix=[];q=[]
    for dataset,norm,model in specs:
        dd=d[(d.dataset==dataset)&(d.normalization==norm)&(d.model==model)].set_index('gene').reindex(genes)
        matrix.append(dd.logFC.to_numpy());q.append(dd.q_junction.to_numpy())
    mat=np.array(matrix).T;im=heat(a,mat,missing='/');a.set(yticks=range(30),yticklabels=genes,xticks=range(8),xticklabels=['B','A']*4)
    for i in range(30):
        for j in range(8):
            if np.isfinite(q[j][i]) and q[j][i]<.05:a.plot(j+.27,i-.27,'o',ms=1.5,mec='none',mfc='white' if abs(mat[i,j])>.6*np.nanmax(abs(mat)) else 'black')
    for j,label in enumerate(['11\nTMM','11\nLib.','10\nTMM','10\nLib.']):text(fig,20+(j+.5)*43/4,115,label,ha='center')
    heading(fig,'A','All specifications',5,5);text(fig,20,14,'B: baseline; A: adjusted')
    row_lines(a);cbar(fig,im,22,127,39,'Macro − micro log2FC')
    adjust=d.query("normalization=='TMM' and model=='source_adjusted'")
    low=adjust.approximate_QL_Wald_low.min();high=adjust.approximate_QL_Wald_high.max();padding=(high-low)*.06
    for axis,dataset,label,x in [(b1,'all_pairs','11 pairs',88),(b2,'same_section','10 same-section',130)]:
        dd=adjust[adjust.dataset.eq(dataset)].set_index('gene').reindex(genes).reset_index()
        forest(axis,dd,['']*30,estimate='logFC',low='approximate_QL_Wald_low',high='approximate_QL_Wald_high',color=[member_color(g) for g in genes])
        axis.set_xlim(low-padding,high+padding);axis.set_xlabel('Adjusted log2FC');axis.tick_params(axis='y',length=0);row_lines(axis);text(fig,x,14,label)
    heading(fig,'B','Adjusted estimates and intervals · TMM',73,5)
    t=tests.query("normalization=='TMM'");components=['HRC','Claudin','Polarity','Junction'];data=[];labels=[]
    for metric,short in [('q_direction','Dir.'),('q_mixed','Mixed')]:
        for comp in components:
            data.append([t[(t.component==comp)&(t.dataset==ds)&(t.model==model)][metric].iloc[0] for ds,model in [('all_pairs','unadjusted'),('all_pairs','source_adjusted'),('same_section','unadjusted'),('same_section','source_adjusted')]])
            labels.append(comp+' / '+short)
    data=np.concatenate([np.array(data[:4]),np.array(data[4:])],axis=1)
    heading(fig,'C','Programme q values · TMM',5,141)
    text(fig,38.5,147,'Directional',ha='center');text(fig,77.5,147,'Mixed',ha='center')
    centers=19+(np.arange(8)+.5)*78/8
    for x,label in zip(centers,['11 B','11 A','10 B','10 A']*2):text(fig,x,153,label,ha='center')
    for i in range(4):
        y=159+i*4
        text(fig,17,y,components[i],ha='right')
        for j in range(8):
            value=data[i,j]
            label=f'{value:.0e}'.replace('e-0','e−') if value<.001 else f'{value:.3f}'
            text(fig,centers[j],y,label,ha='center',color='black',fontweight='bold' if value<.05 else 'normal')
    # A single separator distinguishes the two unchanged test families.
    fig.lines.append(Line2D([58/170,58/170],[1-157/175,1-174/175],transform=fig.transFigure,c='.7',lw=.5))
    table=[];spec4=[('all_pairs','TMM'),('all_pairs','library_only'),('same_section','TMM'),('same_section','library_only')]
    for variable in ['macro','Epithelial','Liver','Endothelial']:
        table.append([f'{vif[(vif.dataset==ds)&(vif.normalization==norm)&(vif.variable==variable)].within_patient_VIF.iloc[0]:.1f}' for ds,norm in spec4])
    table.append([str(int(diag[(diag.dataset==ds)&(diag.normalization==norm)&diag.model.eq('source_adjusted')].residual_df.iloc[0])) for ds,norm in spec4])
    heading(fig,'D','VIF and df',112,141)
    text(fig,132,147,'11 pairs',ha='center');text(fig,154,147,'10 same-section',ha='center')
    for x,label in zip([126.5,137.5,148.5,159.5],['TMM','Library','TMM','Library']):text(fig,x,153,label,ha='center')
    for i,(label,row) in enumerate(zip(['Region','Epithelial','Liver','Endothelial','Residual df'],table)):
        y=159+i*3.2
        text(fig,120,y,label,ha='right')
        for x,value in zip([126.5,137.5,148.5,159.5],row):text(fig,x,y,value,ha='center')
    panel_data('A',d,old_panel='S15a–d',patient_set='11 pairs / 10 same-section subset',normalization='TMM and library_only',model='baseline and source_adjusted',estimand='all 30 member log2FC across eight specifications',multiplicity_family='frozen measured-member BH; q_junction',notes='Complete intervals retained in Table 8; no clipping of heatmap limits.')
    panel_data('B',adjust,patient_set='11 pairs / 10 same-section subset',normalization='TMM',model='source_adjusted',estimand='all 30 member log2FC',ci_method='approximate QL-Wald; identical full-range axes')
    panel_data('C',t,patient_set='11 pairs / 10 same-section subset',normalization='TMM',model='baseline and source_adjusted',estimand='directional and mixed-direction programme tests',multiplicity_family='separate four-programme families for each test/specification')
    panel_data('D',vif,model='source_adjusted',estimand='within-patient VIF and residual model df',ci_method='descriptive diagnostics')
    diag.to_csv(OUT/'source_data/panels/Supplementary_Figure_6_D_model_diagnostics.tsv',sep='\t',index=False,encoding='utf-8')
    save(fig,name)

if __name__=='__main__':
    import sys
    for number in sys.argv[1:] or ['1','2','3','4','5','6']:globals()['supplementary'+number]()
