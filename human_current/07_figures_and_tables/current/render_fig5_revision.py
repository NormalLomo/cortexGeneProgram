#!/usr/bin/env python3
"""Fig5 retained evidence plate with unified panel scale.

The page combines the complete 54-program same-bin matrix with the retained
cross-subclass outline subset of the original 331 pairs in panel a. All marks are drawn from approved retained tables,
the approved program-program median/IQR NPZ and the original Fig5 vector
source. It recovers the original unstandardized regional column order
with the original clustering settings, and applies display-only within-pair
Z scores to f/g/h. It performs no new inference, score fitting, or grouping.
"""
from pathlib import Path
import csv, io, sys, re, textwrap, subprocess
import numpy as np
import pyarrow.parquet as pq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as PlotPath
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.font_manager import FontProperties, findfont, fontManager
fontManager.addfont('/usr/share/fonts/truetype/liberation/LiberationSansNarrow-Regular.ttf')
fontManager.addfont('/usr/share/fonts/truetype/liberation/LiberationSansNarrow-Bold.ttf')
import fitz

ROOT = Path((__import__("os").environ["NMF_WORK_ROOT"] + ""))
CROSS = Path((__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1"))
OUT = ROOT / 'analysis/fig5_revision'
W, H = 720, 867
FINAL_WIDTH_MM = 170.0
FS, SMALL, TITLE, TAG = 10.5, 7.5, 11.4, 13.5
PAIRS = [('P9','P51'), ('P36','P39'), ('P8','P47')]
CHIPS = [('B01012B2','DLPFC'), ('D00865B3','V1'), ('A01186A4','M1')]
REGIONS = ['DLPFC','SMG','M1','V1','AG','VLPFC','SPL','FPPFC','S1']
DONORS = ['DonorB','DonorF','DonorG','DonorH','DonorI']
DC = ['#d58f00','#0072b2','#009e73','#cc79a7','#d55e00']
AREA_COLOURS = ['#7b9cc7','#d5a06e','#87ad8a','#c98289','#a99bc4','#a78471',
                '#c595b5','#939ba3','#b4b77c','#7eafb5','#8290b7','#c9b58e',
                '#83a9a0','#b49d9d']
AREA_STYLES = ['-',(0,(4,2)),(0,(1,1.4)),(0,(5,1.5,1,1.5))]
GLOBAL_COLOUR = '#173e5c'
PREF = ['L2-L3 IT LINC00507','L3-L4 IT RORB','L4-L5 IT RORB','L6 IT','L6 CAR3',
        'L6 CT','L6B','ET','NP','PVALB','CHANDELIER','SST','VIP','NDNF','LAMP5',
        'PAX6','AST','OLIGO','OPC','MICRO','ENDO','VLMC']
PREF_SHORT = {'L2-L3 IT LINC00507':'L2-L3 IT','L3-L4 IT RORB':'L3-L4 IT',
              'L4-L5 IT RORB':'L4-L5 IT','CHANDELIER':'Chand.'}
COL = {**{p:'#4C78A8' for p in PREF[:9]}, **{p:'#8F63B8' for p in PREF[9:16]},
       'AST':'#2E9D58','OLIGO':'#C65353','OPC':'#B7791F',
       'MICRO':'#E88945','ENDO':'#795548','VLMC':'#8C8C8C'}
plt.rcParams.update({'font.family':'Liberation Sans Narrow','font.size':FS,
                     'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,
                     'axes.titlesize':TITLE,'axes.linewidth':.65,
                     'xtick.major.width':.65,'ytick.major.width':.65,
                     'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'none',
                     'figure.facecolor':'none','axes.facecolor':'none'})

def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def key(a,b):
    return frozenset((a,b))

def star(p,ann):
    return p + ('*' if ann[p]['confidence'] == 'Lower confidence' else '')

SHORT = {'P1': 'Alt. mRNA splicing reg.', 'P2': 'Axoneme assembly', 'P3': 'Glu. receptor sig.', 'P4': 'Mito. resp. chain asm.', 'P5': 'Neg. lipid biosynth. reg.', 'P6': 'Pos. secretion reg.', 'P7': 'Cation channel (interneur.)', 'P8': 'Neurofil. cytosk. (pan-neur.)', 'P9': 'ECM organiz. reg.', 'P10': 'SR Ca2+ release reg.', 'P11': 'Postsyn. Glu. sig.', 'P12': 'Oligodend. develop.', 'P13': 'Pos. cation-chan. reg.', 'P14': 'Astrocyte Glu. transp.', 'P15': 'Axon guid./outgrowth', 'P16': 'Axon guid./syn. adhes.', 'P17': 'Glu.-receptor sig.', 'P18': 'Trans-syn. sig. reg.', 'P19': 'Axon guid./neural crest', 'P20': 'Heteroph. cell adhes.', 'P21': 'Syn.-ves. exocyt. (inh.)', 'P22': 'Neuropeptide sig. (exc.)', 'P23': 'Myelination reg.', 'P24': 'Mixed/low-specificity', 'P25': 'Voltage-gated ion chan.', 'P26': 'Postsyn. Glu. receptor a', 'P27': 'Body-fluid secretion', 'P28': 'Axon guid./cell adhes. (L6b)', 'P29': 'Ionotropic Glu. receptor', 'P30': 'Acetylcholine rec. sig.', 'P31': 'Activity-dep. IEG', 'P32': 'Oligodend. (diffuse)', 'P33': 'Oligodend./myelin', 'P34': 'Oligodend. diff.', 'P35': 'Syn. adhes./assembly', 'P36': 'Microglial immune act.', 'P37': 'Ca2+ homeostasis', 'P38': 'Cell secretion reg.', 'P39': 'ECM/OPC proteoglycan', 'P40': 'L-Glu. transport', 'P41': 'Myelination', 'P42': 'NDNF interneur. syn. sig.', 'P43': 'Postsyn. Glu. receptor b', 'P44': 'Sprouting angiogen.', 'P45': 'Microglial immune adhes.', 'P46': 'ECM organization', 'P47': 'Reactive astrocyte/vasc.', 'P48': 'Astrocyte metal homeost.', 'P49': 'Microglial compl./MHC', 'P50': 'Endoth. chemotaxis reg.', 'P51': 'Blood-ves. morphog.', 'P52': 'Neuropeptide sig. (inh.)', 'P53': 'Blood-ves. morphog. (vasc.)', 'P54': 'Microglial chemokine sig.'}

def full(p,ann):
    return SHORT.get(p,ann[p]['functional_name']) + ' (' + star(p,ann) + ')'

def axat(fig,x,y,w,h):
    return fig.add_axes([x/W, 1-(y+h)/H, w/W, h/H])

def tx(fig,x,y,s,size=FS,weight='normal',ha='left',va='baseline',
       color='#202124',rotation=0,**kw):
    return fig.text(x/W, 1-y/H, s, fontsize=size, weight=weight,
                    ha=ha, va=va, color=color, rotation=rotation, **kw)

def tag(fig,x,y,letter):
    tx(fig,x,y,letter,TAG,'bold',color='#111111')

def axes_style(ax, signed=False):
    ax.spines[['top','right']].set_visible(False)
    ax.set_facecolor('none')
    ax.tick_params(length=1.8,pad=1.5,labelsize=SMALL)
    ax.grid(axis='y',color='#e8ebee',lw=.4)
    ax.set_axisbelow(True)

def groups(programs,ann):
    out=[]
    for pref in PREF:
        ix=[i for i,p in enumerate(programs)
            if ann[p]['dominant_subclass']==pref]
        if ix:
            out.append((pref,ix[0],ix[-1]+1))
    return out

def retained_edges(ann,effects,sim):
    # Same old-component upper triangle as the original Fig5d producer.
    raw=sorted(ann,key=lambda p:int(ann[p]['cnmf_component']))
    return [(a,b,effects[(a,b)]) for i,a in enumerate(raw) for b in raw[i+1:]
            if float(sim[key(a,b)]['gene_cosine'])<.25 and effects[(a,b)]>.32]

def heatmap(fig,ann,effects,sim):
    order=sorted(ann,key=lambda p:(PREF.index(ann[p]['dominant_subclass']),int(p[1:])))
    n=len(order);mx,my,side=105.,50.,300.
    gs=groups(order,ann);gap=.3;cell=(side-gap*(len(gs)-1))/n
    starts=[];widths=[];original_index=[];cursor=0.
    group_ranges=[]
    for gi,(pref,a,b) in enumerate(gs):
        left=cursor
        for i in range(a,b):
            starts.append(cursor);widths.append(cell);original_index.append(i);cursor+=cell
        group_ranges.append((pref,a,b,left,cursor))
        if gi<len(gs)-1:
            starts.append(cursor);widths.append(gap);original_index.append(None);cursor+=gap
    bounds=np.r_[0.,np.cumsum(widths)]
    lookup={v:i for i,v in enumerate(original_index) if v is not None}
    centers={i:starts[j]+cell/2 for i,j in lookup.items()}
    values=np.ma.masked_all((len(widths),len(widths)))
    for i,a in enumerate(order):
        for j,b in enumerate(order):values[lookup[i],lookup[j]]=effects[(a,b)]
    norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('white')
    ax=axat(fig,mx,my,side,side)
    ax.pcolormesh(bounds,bounds,values,cmap=cmap,norm=norm,shading='flat',edgecolors='none',rasterized=False)
    ax.set_xlim(0,side);ax.set_ylim(side,0);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_linewidth(.5)
    for gi,(pref,a,b,lo,hi) in enumerate(group_ranges):
        # Gaps occur at identical row and column positions; every data cell remains.
        ax.add_patch(Rectangle((lo,lo),hi-lo,hi-lo,fill=False,edgecolor='#8b9ca5',lw=.4))
        strip=axat(fig,mx+lo,my-4,hi-lo,4)
        strip.set_facecolor('none')
        strip.add_patch(Rectangle((0,0),1,1,transform=strip.transAxes,facecolor=COL[pref],edgecolor='none'))
        strip.set_xticks([]);strip.set_yticks([])
        for sp in strip.spines.values():sp.set_visible(False)
        label_x=mx+(gi+.5)*side/len(gs)
        tx(fig,label_x,my-16,PREF_SHORT.get(pref,pref),SMALL,
           ha='center',va='bottom',rotation=90)
        fig.add_artist(Line2D([label_x/W,(mx+(lo+hi)/2)/W],
            [1-(my-15)/H,1-(my-5)/H],transform=fig.transFigure,
            color='#9ba7af',lw=.4,zorder=2))
    for i,p in enumerate(order):
        tx(fig,mx-5 if i%2==0 else mx+side+24,my+centers[i],full(p,ann),SMALL,
           ha='right' if i%2==0 else 'left',va='center')
    index={p:i for i,p in enumerate(order)};degree={p:0 for p in order}
    edge_set=retained_edges(ann,effects,sim)
    # Display emphasis only: compare known canonical S3 subclass keys.
    display_edges=[(a,b,w) for a,b,w in edge_set
                   if ann[a].get('dominant_subclass') in PREF
                   and ann[b].get('dominant_subclass') in PREF
                   and ann[a]['dominant_subclass'] != ann[b]['dominant_subclass']]
    for a,b,w in display_edges:
        degree[a]+=1;degree[b]+=1
        i,j=sorted((index[a],index[b]))
        # Cross-subclass subset, marked once in the display upper triangle.
        inset=.35
        ax.add_patch(Rectangle((centers[j]-cell/2+inset,centers[i]-cell/2+inset),
            cell-2*inset,cell-2*inset,fill=False,edgecolor='black',lw=.5,zorder=4))
    da=axat(fig,mx+side+4,my,18,side)
    da.barh([centers[i] for i in range(n)],[degree[p] for p in order],height=cell*.78,
            color='#657985',edgecolor='none')
    da.set_ylim(side,0);da.set_xlim(0,max(degree.values())*1.05)
    da.set_yticks([]);da.set_xticks([0,40]);da.tick_params(labelsize=SMALL,length=2,pad=2)
    da.spines[['top','left','right']].set_visible(False)
    tx(fig,mx+side+12,my-7,'Partners',SMALL,ha='center')
    tag(fig,4,14,'a');tx(fig,26,14,'Program matrix',TITLE,'bold')
    tx(fig,10,31,'54 programs · 1,431 pairs',SMALL)
    cax=axat(fig,110,363,120,4.8)
    sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap);sm.set_array([])
    cb=fig.colorbar(sm,cax=cax,orientation='horizontal',ticks=[-2,0,2],extend='both')
    cb.ax.tick_params(labelsize=SMALL,length=1.5,pad=1)
    tx(fig,8,367,'Same-bin log2 median g',SMALL)
    tx(fig,247,369,f'Black boxes: {len(display_edges)} cross-subclass low-loading-overlap positive pairs',SMALL)
    tx(fig,247,381,'Bars: cross-subclass partners; * lower functional-annotation confidence',SMALL)
    print('Retained compact layout: a matrix 300 x 300 native pt; alternating-side short names;',
          'original selected unordered pairs',len(edge_set),
          'displayed cross-subclass black boxes',len(display_edges),flush=True)

def scatter(fig,ann,effects,sim,source):
    tag(fig,542,18,'b');tx(fig,562,18,'Loading similarity',TITLE,'bold')
    ax=axat(fig,555,30,152,94)
    order=sorted(ann,key=lambda p:int(ann[p]['cnmf_component']))
    xy=np.array([(float(sim[key(a,b)]['gene_cosine']),effects[(a,b)])
                 for i,a in enumerate(order) for b in order[i+1:]])
    mask=(xy[:,0]<.25)&(xy[:,1]>.32)
    ax.scatter(xy[:,0],xy[:,1],s=1.0,color='#b8bec3',alpha=.55,linewidths=0,rasterized=True)
    ax.scatter(xy[mask,0],xy[mask,1],s=2.2,color='#b2182b',alpha=.72,linewidths=0,rasterized=True)
    # Recover the unchanged fitted line using the original labelled grid coordinates.
    for drawing in source[0].get_drawings():
        r=drawing['rect'];c=drawing['color']
        if c and c[0]<.2 and c[2]>.3 and r.x0>320 and r.x1<528 and r.y0>20 and r.y1<194 and r.width>150:
            for item in drawing['items']:
                if item[0]=='l':
                    pts=item[1:3]
                    xs=[.25+(q.x-362.9364929)/(2*(470.4082031-362.9364929)) for q in pts]
                    ys=[(160.2611084-q.y)/((160.2611084-60.3684692)/2) for q in pts]
                    ax.plot(xs,ys,color='#1b6ca8',lw=.8)
    ax.axvline(.25,color='#8c969d',lw=.65,ls='--');ax.axhline(.32,color='#8c969d',lw=.65,ls='--')
    ax.set_xlim(.05,1.0);ax.set_ylim(-.55,2.7);ax.set_xticks([.25,.5,.75]);ax.set_yticks([0,1,2])
    ax.set_xlabel('Loading cosine',labelpad=1.5,fontsize=8.5);ax.set_ylabel('Same-bin log2 g',labelpad=2,fontsize=8.5)
    axes_style(ax);ax.tick_params(labelsize=SMALL)
    ax.text(.06,.92,'r = 0.78',transform=ax.transAxes,fontsize=8.5,color='#1b6ca8')

def load_spatial(ann):
    meta=pq.read_table(CROSS/'spatial_bin50_meta.parquet',
                       columns=['bin','chip','x','y','region']).to_pandas()
    meta=meta.loc[meta.chip.isin([c for c,_ in CHIPS])];bins=set(meta.bin)
    mask=pq.read_table(CROSS/'spatial_bin50_rctd_weights.parquet',
                       columns=['bin','rctd_pass_mask']).to_pandas()
    mask=mask.loc[mask.bin.isin(bins)]
    pids=sorted({p for pair in PAIRS for p in pair},key=lambda p:int(p[1:]))
    score=pq.read_table(CROSS/'spatial_bin50_program_score_SCT.parquet',
                        columns=['bin','bin_total_umi']+
                        ['program_'+ann[p]['cnmf_component'] for p in pids]).to_pandas()
    score=score.loc[score.bin.isin(bins)]
    data=meta.merge(mask,on='bin',how='left').merge(score,on='bin',how='left')
    data['valid']=data.rctd_pass_mask.fillna(False).astype(bool)&(data.bin_total_umi.fillna(0)>=200)
    geometry={}
    for chip,_ in CHIPS:
        d=data.loc[data.chip.eq(chip)]
        geometry[chip]=(float(d.x.min()),float(d.x.max()),float(d.y.min()),float(d.y.max()))
    span=max(max(b-a,d-c) for a,b,c,d in geometry.values())*1.04
    return data,geometry,span

def tissue(fig,x,y,chip,p,ann,data,geometry,span,side=76.8):
    ax=axat(fig,x,y,side,side)
    d=data.loc[data.chip.eq(chip)]
    v=d['program_'+ann[p]['cnmf_component']].to_numpy(float)
    good=d.valid.to_numpy()&np.isfinite(v)
    lo,hi=np.percentile(v[good],[1,99]);lo=lo if lo<0 else 0.
    if hi<=lo:hi=lo+1e-6
    size=(side*50/span)**2
    ax.scatter(d.loc[~d.valid,'x'],d.loc[~d.valid,'y'],s=size,marker='s',
               color='#dcdcdc',linewidths=0,rasterized=True)
    ax.scatter(d.loc[good,'x'],d.loc[good,'y'],c=v[good],cmap='viridis',
               vmin=lo,vmax=hi,s=size,marker='s',linewidths=0,rasterized=True)
    a,b,c,e=geometry[chip];cx=(a+b)/2;cy=(c+e)/2
    ax.set_xlim(cx-span/2,cx+span/2);ax.set_ylim(cy+span/2,cy-span/2)
    ax.set_aspect('equal');ax.set_axis_off();ax.set_facecolor('none')
    # Source coordinates remain DNB pixels: 2000 px = 1000 µm = 1 mm.
    ruler_fraction=2000.0/span
    left=.05
    ax.plot([left,left+ruler_fraction],[-.035,-.035],transform=ax.transAxes,
            color='#27323a',lw=.8,clip_on=False,solid_capstyle='butt')
    ax.text(left+ruler_fraction/2,-.085,'1 mm',transform=ax.transAxes,
            ha='center',va='top',fontsize=SMALL,color='#27323a',clip_on=False)

def distance(fig,x,y,pair,profiles,edges,w=64,h=64):
    ax=axat(fig,x,y,w,h)
    (v,lo,hi),regional,area_edges=profiles[pair]
    ax.fill_between(edges[1:],lo,hi,color=GLOBAL_COLOUR,alpha=.055,linewidth=0,zorder=1)
    for ri,values in enumerate(regional):
        if not np.isfinite(values).any():
            raise ValueError(f'Missing regional curve for {pair}, area index {ri}')
        ax.plot(area_edges,values,color=AREA_COLOURS[ri],
                linestyle=AREA_STYLES[ri%len(AREA_STYLES)],lw=.65,alpha=.65,zorder=2)
    ax.plot(edges[1:],v,'o-',color=GLOBAL_COLOUR,ms=2.7,lw=1.6,mew=0,zorder=4)
    ax.axhline(0,color='#8d969d',lw=.7,ls='--',zorder=0)
    # Automatic bounds include every original regional extreme and the global IQR.
    ax.set_xlim(0,515);ax.set_xticks([25,250,500])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4,steps=[1,2,5,10]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'.replace('-','−')))
    ax.set_xlabel('Outer edge (µm)',labelpad=1.5,fontsize=8.5)
    ax.set_ylabel('log2 median g',labelpad=2,fontsize=8.5)
    axes_style(ax);ax.tick_params(labelsize=SMALL)
    pos=ax.get_position()
    print('Retained distance axes',pair,':',pos.width*FINAL_WIDTH_MM,'x',
          pos.height*H/W*FINAL_WIDTH_MM,'mm;',len(regional),'regional lines + all-44 median/IQR',flush=True)
    return ax

def area_key(fig,names,counts):
    handles=[Line2D([],[],color=AREA_COLOURS[i],linestyle=AREA_STYLES[i%len(AREA_STYLES)],
                    lw=.95,label=f'{name} (n={int(counts[i])})') for i,name in enumerate(names)]
    handles.append(Line2D([],[],color=GLOBAL_COLOUR,lw=1.6,label='All 44 (IQR)'))
    tx(fig,8,537,'c–e: regions',SMALL)
    tx(fig,8,548,'(n sections)',SMALL)
    fig.legend(handles=handles,loc='upper left',
               bbox_to_anchor=(100/W,1-551/H,614/W,22/H),mode='expand',ncol=8,
               frameon=False,fontsize=SMALL,handlelength=1.5,handletextpad=.35,
               columnspacing=.8,labelspacing=.35,borderaxespad=0,borderpad=0)
    print('Retained regional lines and source section counts:',list(zip(names,counts)),flush=True)

def map_key(fig,x,y,w=130):
    bar=axat(fig,x,y,w,4)
    bar.imshow(np.linspace(0,1,256)[None,:],cmap='viridis',aspect='auto');bar.set_axis_off()
    tx(fig,x,y-3,'low',SMALL);tx(fig,x+w,y-3,'high',SMALL,ha='right')
    tx(fig,x+w/2-16.5,y+16,'Relative program score',FS,ha='center')

TYPE_BY_PROGRAM={'P9':'Exc','P51':'Endo','P8':'Exc','P47':'Ast','P36':'Micro','P39':'OPC'}
TYPE_COLOUR={'Exc':'#4C78A8','Inh':'#8F63B8','Ast':'#2E9D58','Endo':'#795548','Micro':'#E88945','OPC':'#B7791F'}
# Text hues match the approved Fig4 example palette; Endo uses its corresponding dark tone.
TYPE_TEXT={'Exc':'#2F4A68','Inh':'#593D72','Ast':'#1D6137','Endo':'#4B352D','Micro':'#90552B','OPC':'#714B13'}

def type_strip(fig,x,y,height,pair):
    for j,p in enumerate(pair):
        label=TYPE_BY_PROGRAM[p];yy=y+j*height/2
        fig.add_artist(Rectangle((x/W,1-(yy+height/2)/H),8/W,height/2/H,
            transform=fig.transFigure,facecolor=TYPE_COLOUR[label],edgecolor='none',zorder=3))
        tx(fig,x+4,yy+height/4,label,SMALL,'bold',ha='center',va='center',rotation=90,
           color='#202124' if label in ['Micro','OPC'] else 'white')

def example_header(fig,x,y,p,ann):
    colour=TYPE_TEXT[TYPE_BY_PROGRAM[p]]
    tx(fig,x,y,SHORT[p],8.0,'bold',ha='center',color=colour)
    tx(fig,x,y+11,star(p,ann)+' · '+ann[p]['dominant_subclass'],8.5,'bold',ha='center',color=colour)

def primary_c(fig,ann,data,geom,span,profiles,edges):
    tag(fig,542,159,'c');tx(fig,562,159,'ECM / vessel programs',TITLE,'bold')
    example_header(fig,596.4,176,'P9',ann);example_header(fig,679.4,176,'P51',ann)
    type_strip(fig,548,206,289,('P9','P51'))
    # Only the M1 display row is omitted; loading and shared physical span stay unchanged.
    for row,(chip,area) in enumerate(CHIPS[:2]):
        yy=206+row*94
        tx(fig,596.4,yy-6,area,SMALL,ha='center')
        for col,p in enumerate(['P9','P51']):
            tissue(fig,558+col*83,yy,chip,p,ann,data,geom,span,side=76.8)
    tx(fig,646,404,'P9–P51',FS,'bold',ha='center')
    distance(fig,590,415,('P9','P51'),profiles,edges,w=112,h=80)

def pair_panel(fig,x,y,letter,title,pair,ann,data,geom,span,profiles,edges):
    tag(fig,x+6,y+14,letter);tx(fig,x+28,y+14,title,TITLE,'bold')
    type_strip(fig,x+4,y+46,76.8,pair)
    for j,p in enumerate(pair):
        xx=x+19+j*83;yy=y+46
        example_header(fig,xx+38.4,y+29,p,ann)
        tissue(fig,xx,yy,'D00865B3',p,ann,data,geom,span,side=76.8)
    tx(fig,x+236,y+35,'–'.join(star(p,ann) for p in pair),8.5,'bold',ha='center')
    # The axes use the transparent gutter beside c, not an obsolete panel box.
    distance(fig,x+204,y+48,pair,profiles,edges,w=64,h=64)
    if letter=='d':
        tx(fig,x+236,y+14,'All profiles: 44 sections',SMALL,ha='center')
        tx(fig,x+236,y+25,'median and Q1–Q3',SMALL,ha='center')
    else:
        map_key(fig,x+204,y+6,w=64)

def within_pair_z(vals,axis):
    # Display-only standardization of the exact matrix shown, matching R scale().
    means=np.nanmean(vals,axis=axis,keepdims=True)
    sd=np.nanstd(vals,axis=axis,ddof=1,keepdims=True)
    z=np.divide(vals-means,sd,out=np.full_like(vals,np.nan),where=sd>0)
    finite=z[np.isfinite(z)]
    norm=TwoSlopeNorm(vmin=float(finite.min()),vcenter=0,vmax=float(finite.max()))
    return np.ma.masked_invalid(z),norm

def donor_panel(fig,x,y,ann,per,loo,source):
    tag(fig,x+6,y+4,'f');tx(fig,x+28,y+4,'Donor and section profiles',TITLE,'bold')
    vals=np.array([[float(per[key(*pair)][d]['donor_log2_median_g']) for d in DONORS] for pair in PAIRS])
    z,norm=within_pair_z(vals,axis=1)
    ax=axat(fig,x+29,y+36,156,54)
    im=ax.imshow(z,cmap='RdBu_r',norm=norm,aspect='auto',interpolation='none')
    ax.set_xticks(range(5),['B (1)','F (10)','G (13)','H (16)','I (4)'])
    ax.xaxis.tick_top();ax.set_yticks(range(3),['P9–P51','P36–P39','P8–P47'])
    ax.tick_params(length=0,pad=2,labelsize=SMALL)
    ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True)
    ax.grid(which='minor',color='white',lw=.4);ax.tick_params(which='minor',length=0)
    for j in range(3):
        for i in range(5):
            c=im.cmap(norm(z[j,i]));ink='white' if .2126*c[0]+.7152*c[1]+.0722*c[2]<.5 else '#202124'
            ax.text(i,j,f'{vals[j,i]:.2f}',ha='center',va='center',fontsize=SMALL,color=ink)
    cbax=axat(fig,x+29,y+99,156,4.8)
    cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=[norm.vmin,0,norm.vmax],format='%.2f')
    cb.ax.tick_params(labelsize=SMALL,length=1.5,pad=1)
    bx=axat(fig,x+195,y+36,65,54)
    for j,pair in enumerate(PAIRS):
        r=next(iter(loo[key(*pair)].values()))
        bx.plot([float(r['loo_min_log2_median_g']),float(r['loo_max_log2_median_g'])],[j,j],color='#222222',lw=.7)
        bx.vlines([float(r['loo_min_log2_median_g']),float(r['loo_max_log2_median_g'])],
                  j-.12,j+.12,color='#222222',lw=.7)
        bx.scatter(float(r['all_sections_log2_median_g']),j,marker='D',s=9,color='#222222',zorder=3)
    bx.set_ylim(2.5,-.5);bx.set_xlim(-.5,1.5);bx.set_yticks([]);bx.set_xticks([-.5,0,.5,1,1.5])
    bx.axvline(0,color='#9ba7af',lw=.5,ls='--',zorder=0)
    bx.tick_params(labelsize=SMALL,length=1.5,pad=1);bx.spines[['top','right','left']].set_visible(False)
    bx.set_xlabel('log2 median g',fontsize=SMALL,labelpad=1)
    tx(fig,x+227.5,y+26,'Leave-one-donor-out',SMALL,ha='center')
    tx(fig,x+227.5,y+116,'Point: all 44',SMALL,ha='center')
    tx(fig,x+227.5,y+126,'Range: LOO min–max',SMALL,ha='center')
    tx(fig,x+313,y+14,'Selected-pair distribution',SMALL,ha='center')
    tx(fig,x+29,y+17,'5 donors · within-pair Z; log2 median g values',SMALL)

def region_panel(fig,x,y,ann,regional):
    tag(fig,x+6,y+4,'g');tx(fig,x+28,y+4,'Regional support',TITLE,'bold')
    vals=np.array([[float(regional[key(*pair)][r]['log2_median_g']) for r in REGIONS] for pair in PAIRS])
    z,norm=within_pair_z(vals,axis=1)
    ax=axat(fig,x+23,y+36,326,54)
    im=ax.imshow(z,cmap='RdBu_r',norm=norm,aspect='auto',interpolation='none')
    ax.set_xticks(range(9),REGIONS);ax.xaxis.tick_top();ax.set_yticks(range(3),['–'.join(star(p,ann) for p in pair) for pair in PAIRS])
    ax.tick_params(length=0,pad=2,labelsize=SMALL)
    ax.set_xticks(np.arange(-.5,9,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True)
    ax.grid(which='minor',color='white',lw=.4);ax.tick_params(which='minor',length=0)
    for j in range(3):
        for i in range(9):
            c=im.cmap(norm(z[j,i]));ink='white' if .2126*c[0]+.7152*c[1]+.0722*c[2]<.5 else '#202124'
            ax.text(i,j,f'{vals[j,i]:.2f}',ha='center',va='center',fontsize=SMALL,color=ink)
    cbax=axat(fig,x+23,y+102,320,4.8)
    cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=[norm.vmin,0,norm.vmax],format='%.2f')
    cb.ax.tick_params(labelsize=SMALL,length=1.5,pad=1)
    tx(fig,x+23,y+18,'Within-pair Z-score; numbers: log2 median g',SMALL)

def original_regional_overview(fig,ann,effects,source):
    # Recover the original selected 40 pairs in R upper-triangle traversal order.
    # This is the adopted old producer's selection, not a new example selection.
    excluded={9,18,19,35,52,57}
    raw=sorted((p for p in ann if int(ann[p]['cnmf_component']) not in excluded),
               key=lambda p:int(ann[p]['cnmf_component']))
    upper=[(raw[i],raw[j]) for j in range(len(raw)) for i in range(j)]
    selected=sorted(upper,key=lambda pair:-abs(effects[pair]))[:40]
    with np.load(CROSS/'markcorr_v2/final/progprog_byarea_median_iqr.npz',allow_pickle=True) as d:
        areas=[str(v) for v in d['area_names']]
        stable=[a for a,u in zip(areas,d['unstable']) if int(u)==0]
        aa=[str(v) for v in d['A_names']];bb=[str(v) for v in d['B_names']]
        raw_values=np.array([[float(d['log2_median_g'][areas.index(region),
            aa.index('program_'+ann[a]['cnmf_component']),
            bb.index('program_'+ann[b]['cnmf_component']),0])
            for a,b in selected] for region in stable])
    # The original ComplexHeatmap defaults on the UNSTANDARDIZED full-precision
    # matrix: Euclidean, complete linkage, then mean-branch reordering weighted
    # by negative column means. No clustering is ever applied to Z scores.
    rcode="""options(future.globals.maxSize = 50 * 1024^3)
    m <- as.matrix(read.table(file('stdin'), header=FALSE, sep='\t'))
    tree <- as.dendrogram(stats::hclust(stats::dist(t(m), method='euclidean'), method='complete'))
    tree <- stats::reorder(tree, wts=-colMeans(m), agglo.FUN=mean)
    cat(order.dendrogram(tree), sep=' ')
    """
    payload='\n'.join('\t'.join(format(float(v),'.17g') for v in row) for row in raw_values)+'\n'
    result=subprocess.run(['/usr/local/bin/Rscript','--vanilla','-e',rcode],
                          input=payload,text=True,capture_output=True,check=True)
    order=[int(v)-1 for v in result.stdout.split()]
    rows=['S1','FPPFC','DLPFC','VLPFC','SPL','AG','M1','V1','SMG']
    vals=raw_values[[stable.index(r) for r in rows],:][:,order]
    z,norm=within_pair_z(vals,axis=0)
    tag(fig,10,706,'h');tx(fig,32,706,'40-pair regional overview',TITLE,'bold')
    ax=axat(fig,63,732,398.8,89.73)
    im=ax.pcolormesh(np.arange(41),np.arange(10),z,cmap='RdBu_r',norm=norm,
                     shading='flat',edgecolors='white',linewidth=.35,antialiased=False,rasterized=False)
    ax.set_xlim(0,40);ax.set_ylim(9,0);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)
    # The two existing raw-effect trees are traced from their original vectors,
    # with native final-size stroke weight. They are not recalculated from Z.
    for target,clip in [((63,710,398.8,22),fitz.Rect(65.2,646,269,670.3)),
                        ((3,732,35,89.73),fitz.Rect(9,670.3,38,772.05))]:
        tree_ax=axat(fig,*target);tree_ax.set_xlim(clip.x0,clip.x1);tree_ax.set_ylim(clip.y1,clip.y0)
        tree_ax.set_axis_off()
        for drawing in source[0].get_drawings():
            stroke=drawing['color'];r=drawing['rect']
            if stroke is None or max(stroke)>.35:continue
            # PDF axis-aligned lines have zero-area rectangles; inclusive interval
            # overlap retains those original horizontal and vertical segments.
            if r.x1<clip.x0 or r.x0>clip.x1 or r.y1<clip.y0 or r.y0>clip.y1:continue
            for item in drawing['items']:
                if item[0]=='l':tree_ax.plot([item[1].x,item[2].x],[item[1].y,item[2].y],color=stroke,lw=.5)
    for j,name in enumerate(rows):
        tx(fig,58,732+(j+.5)*89.73/9,name,SMALL,ha='right',va='center',
           color='#E08A00' if name in ['S1','M1','V1'] else '#1B6CA8')
    cbax=axat(fig,134,832.13,327.8,4.8)
    cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=[norm.vmin,0,norm.vmax],format='%.2f')
    cb.ax.tick_params(labelsize=SMALL,length=1.5,pad=1)
    tx(fig,8,838.13,'Within-pair Z-score',SMALL)
    tx(fig,8,856.13,'Trees: original raw effects',SMALL)
    print('Panel h: square cells',398.8/40*FINAL_WIDTH_MM/W,'x',89.73/9*FINAL_WIDTH_MM/W,'mm;',
          'retained matrix',398.8*FINAL_WIDTH_MM/W,'x',89.73*FINAL_WIDTH_MM/W,'mm',flush=True)
    print('Panel h: original 40-pair selection and raw-effect column order recovered;',
          '9 regions retained; within-pair Z-score range',norm.vmin,norm.vmax,flush=True)

def original_ranges(fig,source):
    # Retain the original 22 min/max pairs directly from their vector positions.
    ax=axat(fig,514.2,722,166.4,114)
    ax.set_xlim(331,528);ax.set_ylim(769,641);ax.set_axis_off()
    ax.axvline(333.4,color='#b8bdc2',lw=.7,zorder=0)
    kr=fitz.Rect(331,641,528,769)
    for drawing in source[0].get_drawings():
        r=drawing['rect']
        if not kr.contains(r):continue
        fill=drawing['fill'];stroke=drawing['color']
        if fill and r.width<8 and r.height<8 and max(fill)-min(fill)>.2:
            ax.scatter((r.x0+r.x1)/2,(r.y0+r.y1)/2,s=5,color=fill,zorder=3)
        elif stroke and max(stroke)-min(stroke)<.06:
            for item in drawing['items']:
                if item[0]=='l':ax.plot([item[1].x,item[2].x],[item[1].y,item[2].y],color='#93999e',lw=.5,zorder=1)
    labels=[]
    for b in source[0].get_text('dict')['blocks']:
        for line in b.get('lines',[]):
            for span in line.get('spans',[]):
                r=fitz.Rect(span['bbox']);value=span['text'].strip()
                if 303<r.x0<331 and 638<r.y0<770 and re.match(r'^P\d',value):labels.append((r.y0,value))
    labels.sort()
    for i,(yy,value) in enumerate(labels):
        ypos=722+(yy+3.2-641)/128*114
        tx(fig,507.2 if i%2==0 else 688,ypos,value,SMALL,ha='right' if i%2==0 else 'left',va='center')
    tag(fig,472.2,706,'i');tx(fig,494.2,706,'22-pair regional ranges',TITLE,'bold')
    tx(fig,516.2,842,'0',SMALL,ha='center')
    tx(fig,597.4,851,'Regional log2 g: blue min, red max',8.5,ha='center')

def original_direction(fig,source):
    # Native re-emission of the original visible histogram geometry, not new bins.
    # The legacy pool first excludes the six technical components and intersects
    # the between-chip headline keys (q, effect and Z-direction selection).
    # Its plotted NPZ FSS is raw-g sign agreement, not that Z-direction statistic,
    # and is descriptive after selection rather than independent validation.
    clip=fitz.Rect(140,526,188.5,604)
    drawings=source[0].get_drawings()
    # Earlier obscured page content must remain obscured when the final panel's
    # decorative white backing is omitted. Start after its last covering backdrop.
    back=[i for i,d in enumerate(drawings) if d['fill']==(1.,1.,1.) and d['rect'].contains(clip)]
    start=back[-1]+1 if back else 0
    ax=axat(fig,291,595,52.8,64);ax.set_xlim(clip.x0,clip.x1);ax.set_ylim(clip.y1,clip.y0)
    ax.set_axis_off()
    for drawing in drawings[start:]:
        r=drawing['rect'];fill=drawing['fill'];stroke=drawing['color']
        if r.x1<clip.x0 or r.x0>clip.x1 or r.y1<clip.y0 or r.y0>clip.y1:continue
        if fill==(1.,1.,1.) and r.contains(clip):continue
        verts=[];codes=[]
        def move(q):
            if not verts or tuple(verts[-1])!=(q.x,q.y):verts.append((q.x,q.y));codes.append(PlotPath.MOVETO)
        for item in drawing['items']:
            if item[0]=='re':
                q=item[1]
                verts.extend([(q.x0,q.y0),(q.x1,q.y0),(q.x1,q.y1),(q.x0,q.y1),(q.x0,q.y0)])
                codes.extend([PlotPath.MOVETO,PlotPath.LINETO,PlotPath.LINETO,PlotPath.LINETO,PlotPath.CLOSEPOLY])
            elif item[0]=='l':
                move(item[1]);verts.append((item[2].x,item[2].y));codes.append(PlotPath.LINETO)
            elif item[0]=='c':
                move(item[1])
                for q in item[2:]:verts.append((q.x,q.y));codes.append(PlotPath.CURVE4)
        if verts:
            if drawing.get('closePath') and codes[-1]!=PlotPath.CLOSEPOLY:
                verts.append(verts[0]);codes.append(PlotPath.CLOSEPOLY)
            patch=PathPatch(PlotPath(verts,codes),facecolor=fill if fill is not None else 'none',
                            edgecolor=stroke if stroke is not None else 'none',
                            linewidth=.18 if fill is not None else .4)
            dash=re.search(r'\[([^]]+)\]',drawing.get('dashes',''))
            if dash:
                vals=[float(v) for v in dash.group(1).split()]
                if vals:patch.set_linestyle((0,vals))
            ax.add_patch(patch)
    for val,sy in [('0',603.7),('300',582.7),('600',561.7),('900',540.7)]:
        tx(fig,286,595+(sy-526)/78*64+2,val,SMALL,ha='right')
    for val,sx in [('0.85',143.8),('0.90',158.7),('0.95',173.6),('1.00',188.5)]:
        tx(fig,291+(sx-140)/48.5*52.8,668,val,SMALL,ha='center')
    tx(fig,268,627,'Selected pairs',SMALL,ha='center',va='center',rotation=90)
    tx(fig,317.4,588,'Median 1.00',SMALL,ha='center')
    tx(fig,317.4,679,'Same-direction fraction',SMALL,ha='center')
    tx(fig,317.4,689,'44 sections',SMALL,ha='center')

def main():
    source=fitz.open(stream=sys.stdin.buffer.read(),filetype='pdf')
    ann={r['new_P']:r for r in read(ROOT/'tables/TableS3_program_annotation.tsv')}
    old={r['cnmf_component']:r['new_P'] for r in ann.values()}
    sim={}
    for r in read(CROSS/'markcorr_v2/similarity/program_similarity_pairs.tsv'):
        a,b=r['A'].removeprefix('program_'),r['B'].removeprefix('program_')
        if a in old and b in old:sim[key(old[a],old[b])]=r
    with np.load(CROSS/'markcorr_v2/final/progprog_median_iqr.npz',allow_pickle=True) as d:
        aa=list(d['A_names']);bb=list(d['B_names'])
        ia={p:aa.index('program_'+ann[p]['cnmf_component']) for p in ann}
        ib={p:bb.index('program_'+ann[p]['cnmf_component']) for p in ann}
        effects={(a,b):float(d['log2_median_g'][ia[a],ib[b],0]) for a in ann for b in ann}
        edges=d['ring_edges_um'].copy()
        profiles={pair:tuple(d[f][ia[pair[0]],ib[pair[1]],:].copy()
                             for f in ['log2_median_g','log2_q1','log2_q3'])
                  for pair in PAIRS}
    # Read the retained program-program area curves, never the Fig4 cell-program product.
    with np.load(CROSS/'markcorr_v2/final/progprog_byarea_median_iqr.npz',allow_pickle=True) as d:
        area_names=[str(v) for v in d['area_names']]
        area_counts=d['n_chips_area'].tolist()
        area_edges=d['ring_edges_um'][1:].copy()
        aa=list(d['A_names']);bb=list(d['B_names'])
        for pair in PAIRS:
            ia=aa.index('program_'+ann[pair[0]]['cnmf_component'])
            ib=bb.index('program_'+ann[pair[1]]['cnmf_component'])
            profiles[pair]=(profiles[pair],d['log2_median_g'][:,ia,ib,:].copy(),area_edges)
    per,loo={},{}
    for fn,target,field in [('same_bin_per_donor_effects.tsv',per,'donor'),
                            ('same_bin_true_donor_loo.tsv',loo,'omitted_donor')]:
        for r in read(ROOT/'analysis/02_true_donor_spatial_support'/fn):
            if r['mode']=='progprog' and key(r['A_label'],r['B_label']) in [key(*p) for p in PAIRS]:
                target.setdefault(key(r['A_label'],r['B_label']),{})[r[field]]=r
    regional={}
    with (CROSS/'markcorr_v2/final/progprog_byarea_median_iqr.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            a,b=r['A'].removeprefix('program_'),r['B'].removeprefix('program_')
            if a in old and b in old and int(a)<int(b) and float(r['ring_um'])==25:
                k=key(old[a],old[b])
                if k in [key(*p) for p in PAIRS]:regional.setdefault(k,{})[r['area']]=r
    # The legacy area TSV contains only effect-threshold-selected pairs.
    # The approved modest-negative example uses the same retained estimator
    # directly from the already loaded complete area NPZ, without reaggregation.
    regional[key('P8','P47')]={area:{'log2_median_g':float(values[0])}
        for area,values in zip(area_names,profiles[('P8','P47')][1])}
    data,geom,span=load_spatial(ann)
    fig=plt.figure(figsize=(W/72,H/72),facecolor='none')
    fig.patch.set_alpha(0)
    heatmap(fig,ann,effects,sim);scatter(fig,ann,effects,sim,source)
    primary_c(fig,ann,data,geom,span,profiles,edges)
    pair_panel(fig,4,384,'d','V1 · Modest avoidance',('P8','P47'),ann,data,geom,span,profiles,edges)
    pair_panel(fig,274,384,'e','V1 · Immune / OPC-ECM',('P36','P39'),ann,data,geom,span,profiles,edges)
    area_key(fig,area_names,area_counts)
    donor_panel(fig,4,563,ann,per,loo,source)
    region_panel(fig,366,563,ann,regional)
    original_regional_overview(fig,ann,effects,source)
    original_ranges(fig,source)
    original_direction(fig,source)
    source_overlays_doc=io.BytesIO()
    fig.savefig(source_overlays_doc,format='pdf',dpi=300,transparent=True,facecolor='none')
    plt.close(fig)
    doc=fitz.open(stream=source_overlays_doc.getvalue(),filetype='pdf')
    # Coordinate-unit conversion after native reflow: physical output width is 170 mm.
    physical=fitz.open()
    page=physical.new_page(width=FINAL_WIDTH_MM/25.4*72,height=H/W*FINAL_WIDTH_MM/25.4*72)
    page.show_pdf_page(page.rect,doc,0)
    doc.close();doc=physical
    OUT.mkdir(parents=True,exist_ok=True)
    pdf=OUT/'Fig5_revised.pdf';png=OUT/'Fig5_revised.png'
    doc.save(pdf,garbage=3,deflate=True)
    pix=doc[0].get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=True)
    pix.save(png)
    print('Produced',pdf,flush=True);print('Produced',png,flush=True)
    print('Retained PDF page:',doc[0].rect.width,'x',doc[0].rect.height,'pt; PNG:',pix.width,'x',pix.height,'px',flush=True)
    print('Retained single-page figure:',FINAL_WIDTH_MM,'x',H/W*FINAL_WIDTH_MM,'mm;',
          'minimum final font:',SMALL*FINAL_WIDTH_MM/(W/72*25.4),'pt;',
          'PDF bytes:',pdf.stat().st_size,'PNG bytes:',png.stat().st_size,
          'PNG alpha:',pix.alpha,'channels:',pix.n,flush=True)
    doc.close();source.close()

if __name__=='__main__':
    main()
