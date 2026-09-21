#!/usr/bin/env python3
"""Approved six-main-figure composition, preserving source vector marks.

Latest targeted mapping: Fig1 a/b/c retained, old e -> d; rank plot remains S1. Fig2 a-c retained with condensed display. Fig3 unchanged. Fig4 old a/b/c/d retained, old f -> e; old e -> S14d. Fig5 old b/c/e/f -> a/b/c/d; old a -> S15e; old d -> S16d. Fig6 old a/b retained; old c/d -> S17g/h. Current supplements end at S19; the duplicate former S19 is withdrawn and the former S20 Neurosynth figure is S19.

Prior mapping (for source identification):
Fig1 a=new human workflow; b=old1b text repair; c=old1c; d=S1a; e=S2c.
Fig2 current native-cleaning and local repair actions only; no retired subset panels.
Fig3 a=old3a; b=old3b; c=old3f; d=old3g; e=old3h (inclusive reference).
Fig4 a=old4a text repair; b=old4b; c=old4d; d=old4c;
     e=old4e label repair; f=old4f V1 row + old4g overlays.
Fig5 a=old5a; b=old5b; c=old5g; d=old5h left section histogram ONLY;
     e=A P16 same-bin donor/LOO; f=A P16 0-500 donor/LOO.
     Old5h right pseudo-donor LOO is superseded, not relabelled or retained as evidence.
Fig6 a=old6a; b=old6f; c=old6e; d=old6j (current6 is historical8).
Current figures and supplemental panels are produced by their dedicated current
producers; this file contains only explicit native-PDF repair actions.
"""
from pathlib import Path
import fitz
ROOT=Path((__import__("os").environ["NMF_WORK_ROOT"] + ""))
SRC=ROOT/'inputs/current_six_figures';OUT=ROOT/'figures/human_revision';PAN=OUT/'panels'
OVERALL_BOX=ROOT/'analysis/external_single_cell_program_scores/cell_nearest_reference_cosine_overall_boxplot.pdf'
F_HALF=ROOT/'figures/human_revision/panels/Fig1_f_halfwidth.pdf'
CACHE={}
CROPS={
 '1c':(350,172,547,319),'1d':(352,318,547,541),'1e':(0,540,202,736),'1f':(199,540,547,736),
 '2a':(0,0,547,211),'2b':(0,210,96,360),'2c':(94,210,222,360),'2d':(220,210,385,360),'2e':(382,210,547,360),'2f':(0,359,547,506),'2g':(0,505,166,661),'2h':(163,508,369,661),'2i':(369,508,547,661),
 '3a':(0,0,503.643,138),'3b':(0,138,122,267),'3c':(121,138,253,267),'3d':(250,138,367,267),'3e':(365,138,503.643,267),'3f':(0,266,310,533),'3g':(0,531,322,669.818),'3h':(306,273,503.643,492),'3i':(322,491,503.643,669.818),
 '4b':(0,283,112,408),'4c':(112,284,207,408),'4d':(207,284,339,408),'4f':(0,407,224,600),'4fv1':(0,477,224,532),'4fscale':(38,586,224,600),'4g':(312,430,536,531),'4h':(212,526,354,601),'4i':(354,531,539,606),'4j':(0,588,133,721),'4k':(132,594,281,721),'4l':(280,601,539,721),
 '5a':(0,0,287,246),'5b':(287,0,532,219),'5c':(0,244,276,502),'5d':(276,219,532,361),'5e':(274,360,428,503),'5f':(427,354,532,502),'5g':(0,504,130,628),'5hist':(130,522,194,627),'5i':(290,502,532,627),'5j':(0,621,290,791.136),'5k':(290,627,532,791.136),
 '6a':(0,0,253,247),'6b':(250,0,510.236,207),'6c':(0,246,250,328),'6d':(250,207,341,325),'6e':(340,207,510.236,326),'6f':(0,325,290,505),'6g':(289,328,510.236,507),'6h':(0,501,168,629.291),'6i':(165,509,386,629.291),'6j':(383,507,510.236,629.291),
 'S1a':(0,0,350,224),'S2c':(0,218,277,378)
}

def source(key):
    if key.startswith('S'):n='Fig'+key[:2]+'.pdf'
    else:n='Fig'+key[0]+'.pdf'
    return SRC/n,fitz.Rect(CROPS[key])

def getdoc(path):
    path=Path(path)
    if path not in CACHE:
        d=fitz.open(path)
        if path.parent==SRC and path.name!='Supplementary_Figures.pdf':
            p=d[0]
            for b in p.get_text('dict')['blocks']:
                if b['type']!=0:continue
                for l in b['lines']:
                    for s in l['spans']:
                        if s['text'].strip() in list('abcdefghijkl') and s['size']>=7:
                            p.add_redact_annot(fitz.Rect(s['bbox']),fill=False,cross_out=False)
            try:p.apply_redactions(images=0,graphics=0)
            except TypeError:p.apply_redactions(images=0)
        CACHE[path]=d
    return CACHE[path]

def newpage(w,h,title):
    d=fitz.open();p=d.new_page(width=w,height=h)
    p.insert_text((18,23),title,fontsize=14,fontname='hebo',color=(.08,.12,.19))
    return d,p

def place(p,letter,title,box,key=None,path=None,clip=None):
    x,y,w,h=box
    if letter:p.insert_text((x,y+13),letter,fontsize=14,fontname='hebo')
    if title:p.insert_text((x+18,y+12),title,fontsize=9.8,fontname='hebo',color=(.15,.18,.22))
    if key:path,clip=source(key)
    d=getdoc(path);r=fitz.Rect(clip) if clip is not None else d[0].rect
    scale=min(w/r.width,(h-23)/r.height)
    p.show_pdf_page(fitz.Rect(x,y+23,x+r.width*scale,y+23+r.height*scale),d,0,clip=r,keep_proportion=True)

def finish(d,n,supp=False):
    pdfdir=OUT/('supplementary_figures_pdf' if supp else 'main_figures_pdf');pngdir=OUT/('supplementary_figures_png' if supp else 'main_figures_png')
    pdfdir.mkdir(parents=True,exist_ok=True);pngdir.mkdir(parents=True,exist_ok=True)
    d.save(pdfdir/(n+'.pdf'),garbage=4,deflate=True)
    d[0].get_pixmap(dpi=300,alpha=False).save(pngdir/(n+'.png'))
    print('Produced',pdfdir/(n+'.pdf'),pngdir/(n+'.png'),flush=True)

def _fig3_clean_source():
    d=fitz.open(SRC/'Fig3.pdf')
    p=d[0]
    for clip in [fitz.Rect(CROPS['3a']),fitz.Rect(CROPS['3f']),fitz.Rect(CROPS['3h'])]:
        for block in p.get_text('dict',clip=clip)['blocks']:
            if block['type']!=0:
                continue
            for line in block['lines']:
                for s in line['spans']:
                    r=fitz.Rect(s['bbox'])
                    if r.intersects(clip):
                        r.x0-=.12;r.x1+=.12;r.y0-=.08;r.y1+=.08
                        p.add_redact_annot(r,fill=False,cross_out=False)
    try:
        p.apply_redactions(images=0,graphics=0)
    except TypeError:
        p.apply_redactions(images=0)
    return d

def _fig3_header(p,letter,title,x,y):
    p.insert_text((x,y+14),letter,fontsize=12,fontname='hebo')
    p.insert_text((x+17,y+12),title,fontsize=8.5,fontname='hebo',color=(.15,.18,.22))

def _fig3_show(p,doc,clip,box):
    x,y,w,h=box
    r=fitz.Rect(clip)
    scale=min(w/r.width,(h-23)/r.height)
    p.show_pdf_page(fitz.Rect(x,y+23,x+r.width*scale,y+23+r.height*scale),doc,0,clip=r,keep_proportion=True)
    return scale

def _fig3_putbox(p,rect,text,size=7,align=fitz.TEXT_ALIGN_LEFT,color=(.10,.10,.10),font='helv'):
    p.insert_textbox(fitz.Rect(rect),text,fontsize=size,fontname=font,align=align,color=color)

def build_fig3():
    d=fitz.open()
    p=d.new_page(width=547,height=950)
    clean=_fig3_clean_source()

    _fig3_header(p,'a','Representative program maps',8,40)
    aclip=fitz.Rect(CROPS['3a'])
    a_scale=_fig3_show(p,clean,aclip,(8,40,531,170))
    ax,ay=8,63
    p.insert_text((ax+9*a_scale,ay+5*a_scale),'Spatial program activity',fontsize=7.5,fontname='hebo')
    a_headers=[
        (46,'Excitatory\nP31\nActivity-dependent\nIEG'),
        (128,'Inhibitory (VIP)\nP21\nSynaptic-vesicle\nexocytosis (inh)'),
        (210,'Astrocyte\nP14 Astrocyte\nglutamate transport'),
        (292,'Oligodendrocyte\nP33\nOligodendrocyte/myelin'),
        (374,'Endothelial\nP51 Blood vessel\nmorphogenesis'),
        (456,'Microglia\nP36 Microglial\nimmune activation'),
    ]
    for cx,txt in a_headers:
        f='helv'
        _fig3_putbox(p,(ax+(cx-38)*a_scale,ay+10*a_scale,ax+(cx+38)*a_scale,ay+50*a_scale),txt,7,fitz.TEXT_ALIGN_CENTER,font=f)
        if cx==128:
            tw=fitz.get_text_length('Inhibitory (VIP)',fontname='helv',fontsize=7)
            xr=ax+cx*a_scale+tw/2+1
            yt=ay+10*a_scale+6
            p.draw_line((xr,yt-6),(xr,yt+1),color=(.10,.10,.10),width=.65)
            p.draw_line((xr-2.3,yt-2.2),(xr+2.3,yt-2.2),color=(.10,.10,.10),width=.65)
    for sx,txt in [(396,'low'),(427,'0'),(452,'high')]:
        p.insert_text((ax+sx*a_scale,ay+126*a_scale),txt,fontsize=7)
    p.insert_text((ax+331*a_scale,ay+137*a_scale),'standardized activity (per program)',fontsize=7)

    _fig3_header(p,'b','Same-section tissue domains',8,222)
    bpath,bclip=source('3b')
    _fig3_show(p,getdoc(bpath),bclip,(8,222,158,175))

    _fig3_header(p,'c','All-program domain profiles',174,222)
    cgraphic=fitz.Rect(4,284,223,520)
    cscale=1.20
    p.show_pdf_page(fitz.Rect(224,300,224+219*cscale,300+236*cscale),clean,0,clip=cgraphic,keep_proportion=True)
    def cy(src_y):
        return 300+(src_y-284)*cscale+6
    left_ids=[
        (293.6,'P32'),(310.8,'P12'),(319.5,'P41'),(328.1,'P47'),(336.8,'P14'),
        (345.4,'P34'),(354.0,'P39'),(362.7,'P38'),(371.3,'P7'),(380.0,'P46'),
        (388.6,'P20'),(397.2,'P17'),(405.9,'P18'),(414.5,'P2'),(423.2,'P3'),
        (449.1,'P30'),(466.4,'P43'),(475.0,'P42'),(483.6,'P27'),(492.3,'P35'),
        (500.9,'P44'),(509.6,'P6')
    ]
    right_ids=[
        (284.9,'P15'),(293.6,'P24'),(302.2,'P1'),(310.8,'P25'),(319.5,'P29'),
        (328.1,'P28'),(336.8,'P11'),(345.4,'P19'),(354.0,'P37'),(362.7,'P10'),
        (371.3,'P22'),(380.0,'P13'),(388.6,'P9'),(397.2,'P16'),(405.9,'P52'),
        (414.5,'P40'),(423.2,'P50'),(431.8,'P8'),(440.5,'P21'),(449.1,'P54'),
        (457.8,'P36'),(466.4,'P45'),(475.0,'P4'),(483.7,'P49'),(492.3,'P31'),
        (509.6,'P53')
    ]
    for sy,txt in left_ids:
        p.insert_text((291,cy(sy)),txt,fontsize=6.5)
    for sy,txt in right_ids:
        p.insert_text((489,cy(sy)),txt,fontsize=6.5)
    for sy,txt in [
        (285.5,'P23 Myelination'),(302.5,'P33 Myelin'),
        (432.0,'P5 L2/3 lipid'),(441.0,'P48 Astro'),(458.0,'P26 LAMP5')
    ]:
        _fig3_putbox(p,(176,cy(sy)-7,276,cy(sy)+10),txt,6.5,fitz.TEXT_ALIGN_RIGHT)
    _fig3_putbox(p,(489,cy(500.0)-7,538,cy(500.0)+10),'P51 Vessel',6.5,fitz.TEXT_ALIGN_LEFT)
    for sx,txt in [
        (73.8,'NOID'),(81.3,'L1'),(88.7,'L2'),(96.2,'L3'),(103.7,'L4'),(111.2,'L5'),(118.6,'L6'),(126.1,'WM'),
        (152.2,'NOID'),(159.7,'L1'),(167.2,'L2'),(174.6,'L3'),(182.1,'L4'),(189.6,'L5'),(197.1,'L6'),(204.5,'WM')
    ]:
        p.insert_text((224+(sx-4)*cscale,594),txt,fontsize=6.5,rotate=90)

    peakclip=fitz.Rect(248,328,288,385)
    p.show_pdf_page(fitz.Rect(224,600,264,657),clean,0,clip=peakclip,keep_proportion=True)
    p.insert_text((224,608),'peak cell',fontsize=7,fontname='hebo')
    for i,(l,r) in enumerate([('Exc','OPC'),('Inh','Imm'),('Ast','Vas'),('Olig','Other')]):
        yy=624+i*9
        p.insert_text((233,yy),l,fontsize=6.5)
        p.insert_text((264,yy),r,fontsize=6.5)
    zclip=fitz.Rect(248,394,284,474)
    p.show_pdf_page(fitz.Rect(380,614,416,674),clean,0,clip=zclip,keep_proportion=True)
    p.insert_text((380,608),'median z',fontsize=7,fontname='hebo')
    for yy,txt in [(624,'1.2'),(637,'0.6'),(650,'0'),(663,'-0.6'),(676,'-1.2')]:
        p.insert_text((397,yy),txt,fontsize=6.5)

    _fig3_header(p,'e','Across-section profile consistency',8,405)
    eclip=fitz.Rect(CROPS['3h'])
    escale=_fig3_show(p,clean,eclip,(8,405,158,210))
    ex,ey=8,428
    def eypos(src_y):
        return ey+(src_y-273)*escale
    p.insert_text((ex+(335.5-306)*escale,eypos(276.7)+7),'Cross-section profile consistency',fontsize=6.5)
    p.insert_text((8,506),'program',fontsize=7,rotate=90)
    for i,txt in enumerate(['P23','P32','P8','P33']):
        p.insert_text((14,452+i*8),txt,fontsize=6.5)
    for i,txt in enumerate(['P31','P51','P36','P48','P26','P5']):
        p.insert_text((14,488+i*8),txt,fontsize=6.5)
    for i,txt in enumerate(['P14','P34','P39','P47']):
        p.insert_text((14,540+i*8),txt,fontsize=6.5)
    for sx,txt in [(347.3,'-0.5'),(384.7,'0.0'),(421.3,'0.5'),(457.9,'1.0')]:
        p.insert_text((ex+(sx-306)*escale,585),txt,fontsize=6.5)
    p.insert_text((ex+(364.4-306)*escale,589),'correlation to mean layer profile',fontsize=6.5)
    p.insert_text((ex+(469.6-306)*escale,eypos(353.3)+7),'median r',fontsize=6.5)
    for sy,txt in [(359.7,'1.00'),(367.6,'0.95'),(375.5,'0.90'),(383.4,'0.85'),(391.3,'0.80')]:
        p.insert_text((ex+(484.7-306)*escale,eypos(sy)+6),txt,fontsize=6.5)

    _fig3_putbox(p,(8,620,166,678),'Inclusive mean reference; internal section-profile consistency.\n\nDistance is from high-OLIGO bins, not a reconstructed white-matter boundary.',7,fitz.TEXT_ALIGN_LEFT,color=(.25,.28,.32))

    _fig3_header(p,'d','Physical distance from high-OLIGO reference bins',8,690)
    dpanel=PAN/'fig4a_label_repair/Fig3_distance_clean.pdf'
    _fig3_show(p,getdoc(dpanel),getdoc(dpanel)[0].rect,(8,690,531,252))
    finish(d,'Fig3')

def build_fig1_original():
    d,p=newpage(600,845,'Figure 1 | Human program repertoire and internal stability')
    place(p,'a','Human-only study workflow',(10,38,580,125),path=PAN/'fig1_human_workflow/Fig1_human_workflow.pdf')
    place(p,'b','Program usage across subclasses',(10,177,320,635),path=PAN/'fig1b_label_repair/Fig1b.pdf')
    place(p,'c','Program-usage embedding',(342,177,245,215),key='1c')
    place(p,'d','Internal cNMF stability',(342,430,245,205),key='S2c')
    p.insert_textbox(fitz.Rect(346,671,581,737),'Names, leading genes and annotation confidence: Table S2.\nRank sensitivity remains in Fig. S1.',fontsize=9)
    p.insert_text((12,831),'Internal stability and rank sensitivity are not independent validation.',fontsize=9)
    finish(d,'Fig1')

def build_fig1_with_overall_box():
    d,p=newpage(600,970,'Figure 1 | Human program repertoire and internal stability')
    place(p,'a','Human-only study workflow',(10,38,580,125),path=PAN/'fig1_human_workflow/Fig1_human_workflow.pdf')
    place(p,'b','Program usage across subclasses',(10,177,320,635),path=PAN/'fig1b_label_repair/Fig1b.pdf')
    place(p,'c','Program-usage embedding',(342,177,245,215),key='1c')
    place(p,'d','Internal cNMF stability',(342,430,245,205),key='S2c')
    place(p,'e','Overall cell-pooled nearest-reference similarity',(342,671,245,145),path=OVERALL_BOX)
    p.insert_textbox(fitz.Rect(346,830,581,900),'Names, leading genes and annotation confidence: Table S2.\nRank sensitivity remains in Fig. S1.',fontsize=9)
    p.insert_text((12,955),'Internal stability and rank sensitivity are not independent validation.',fontsize=9)
    finish(d,'Fig1')

def build_fig1_local_fg():
    """Replace only the original f slot with native half-width f and g vectors."""
    source_doc = fitz.open(SRC/'Fig1.pdf')
    f_doc = fitz.open(F_HALF)
    box_doc = fitz.open(OVERALL_BOX)
    d = fitz.open()
    d.insert_pdf(source_doc, from_page=0, to_page=0)
    p = d[0]

    # The locked original page is retained as the full-page base.  Only the
    # bottom-right f slot is covered; a-e and the 547 x 736 pt page stay intact.
    f_region = fitz.Rect(189, 540, 547, 736)
    p.draw_rect(f_region, color=None, fill=(1, 1, 1), overlay=True)

    # Both source components are produced at the exact 175 x 196 pt half-slot
    # size; show_pdf_page therefore places each vector 1:1 without shrinkage.
    left_half = fitz.Rect(197, 540, 372, 736)
    right_half = fitz.Rect(372, 540, 547, 736)
    p.show_pdf_page(left_half, f_doc, 0, clip=f_doc[0].rect,
                    keep_proportion=True, overlay=True)
    p.show_pdf_page(right_half, box_doc, 0, clip=box_doc[0].rect,
                    keep_proportion=True, overlay=True)

    # Panel letters are added once at the new half-slot corners.
    p.insert_text((202, 556), 'f', fontsize=14, fontname='hebo', color=(0, 0, 0))
    p.insert_text((377, 556), 'g', fontsize=14, fontname='hebo', color=(0, 0, 0))

    finish(d, 'Fig1')
    box_doc.close()
    f_doc.close()
    source_doc.close()


def finish_fig2(d):
    """Serialize only this Fig2 branch without high-level garbage collection."""
    pdfdir = OUT / 'main_figures_pdf'
    pngdir = OUT / 'main_figures_png'
    pdfdir.mkdir(parents=True, exist_ok=True)
    pngdir.mkdir(parents=True, exist_ok=True)
    d.save(pdfdir / 'Fig2.pdf', garbage=1, deflate=True)
    d[0].get_pixmap(dpi=300, alpha=False).save(pngdir / 'Fig2.png')
    print('Produced', pdfdir / 'Fig2.pdf', pngdir / 'Fig2.png', flush=True)



def build_fig2_subclass_reorder():
    """Approved a/b/c/d/e/f order and current clear regional matrix in S11h."""
    main_pdf = OUT / "main_figures_pdf" / "Fig2.pdf"
    s11_pdf = OUT / "supplementary_figures_pdf" / "FigS11.pdf"
    bundle_pdf = OUT / "Supplementary_Figures.pdf"
    # Read all three retained sources before any destination can be overwritten.
    old_main = fitz.open(stream=main_pdf.read_bytes(), filetype="pdf")
    old_s11 = fitz.open(stream=s11_pdf.read_bytes(), filetype="pdf")
    old_bundle = fitz.open(stream=bundle_pdf.read_bytes(), filetype="pdf")
    support = fitz.open(PAN / "fig2_donor_evidence" / "Fig2_within_subclass_regional_support.pdf")

    def native(page, source, clip, destination):
        page.show_pdf_page(fitz.Rect(destination), source, 0,
                           clip=fitz.Rect(clip), keep_proportion=True)

    # Only the bottom g/h/i area is reflowed to retain h at readable native size.
    s11 = fitz.open()
    s = s11.new_page(width=1040, height=1810)
    native(s, old_s11, (0, 0, 1040, 1208), (0, 0, 1040, 1208))
    native(s, old_s11, (15, 1214, 335, 1240), (15, 1214, 335, 1240))
    native(s, old_s11, (15, 1375, 335, 1705), (15, 1250, 335, 1580))
    s.insert_text((350, 1233), "h", fontsize=14, fontname="hebo")
    s.insert_text((368, 1232), "Regional program activity (within-program z score)",
                  fontsize=9.8, fontname="hebo", color=(.15, .18, .22))
    # Uniform enlargement, preserving the current 14 x 38 matrix, order and marks.
    c_clip = fitz.Rect(8, 278, 539, 443)
    c_scale = 675 / c_clip.width
    c_bottom = 1250 + c_clip.height * c_scale
    native(s, old_main, c_clip, (350, 1250, 1025, c_bottom))
    s.insert_text((350 + (521-8)*c_scale, 1245), "z-score",
                  fontsize=6.5*c_scale, fontname="helv", color=(.13, .15, .17))
    s.insert_text((350, 1480), "Black outline: each region's top positive regional z score.",
                  fontsize=8.5, fontname="helv", color=(.2, .23, .26))
    s.insert_text((350, 1492), "Descriptive area profiles; not donor-adjusted significance.",
                  fontsize=8.5, fontname="helv", color=(.2, .23, .26))
    s.insert_text((350, 1528), "i", fontsize=14, fontname="hebo")
    s.insert_text((368, 1527), "Original lobe-level rankings",
                  fontsize=9.8, fontname="hebo", color=(.15, .18, .22))
    native(s, old_s11, (735, 1395, 1040, 1650), (350, 1540, 655, 1795))

    # Assemble the new main figure from unchanged native retained blocks.
    result = fitz.open()
    p = result.new_page(width=547, height=729)
    native(p, old_main, (0, 0, 547, 248), (0, 0, 547, 248))
    native(p, old_main, (8, 444, 539, 592), (8, 252, 539, 400))
    native(p, old_main, (0, 596, 547, 729), (0, 404, 547, 537))
    native(p, support, support[0].rect, (8, 539, 539, 729))
    spans = []
    for block in old_main[0].get_text("dict")["blocks"]:
        if block.get("type") == 0:
            for line in block["lines"]:
                spans.extend(line["spans"])
    for old_letter, new_letter, bounds in [
        ("d", "c", fitz.Rect(8, 444, 27, 465)),
        ("e", "d", fitz.Rect(8, 596, 27, 619)),
        ("f", "e", fitz.Rect(297, 596, 316, 619)),
    ]:
        mark = next(span for span in spans
                    if span["text"].strip() == old_letter and span["size"] >= 12
                    and fitz.Rect(span["bbox"]).intersects(bounds))
        rect = fitz.Rect(mark["bbox"]) + (-.7, -192.7, .7, -191.3)
        p.draw_rect(rect, color=None, fill=(1, 1, 1), overlay=True)
        ox, oy = mark["origin"]
        p.insert_text((ox, oy-192), new_letter, fontsize=13, fontname="hebo",
                      color=(0, 0, 0), overlay=True)

    combined = fitz.open()
    combined.insert_pdf(old_bundle, from_page=0, to_page=9)
    combined.insert_pdf(s11)
    combined.insert_pdf(old_bundle, from_page=11, to_page=len(old_bundle)-1)

    # Retain the clear original C in its requested formal destination first.
    s11.save(s11_pdf, garbage=1, deflate=True)
    s11_png = OUT / "supplementary_figures_png" / "FigS11.png"
    s.get_pixmap(dpi=300, alpha=False).save(s11_png)
    print("Produced", s11_pdf, s11_png, flush=True)
    combined.save(bundle_pdf, garbage=1, deflate=True)
    print("Produced", bundle_pdf, flush=True)
    result.save(main_pdf, garbage=1, deflate=True)
    main_png = OUT / "main_figures_png" / "Fig2.png"
    p.get_pixmap(dpi=300, alpha=False).save(main_png)
    print("Produced", main_pdf, main_png, "successful Fig2 production 6", flush=True)
    for document in (combined, result, s11, support, old_bundle, old_s11, old_main):
        document.close()
def _fig5_old_source():
    d = fitz.open(SRC / 'Fig5.pdf')
    p = d[0]
    allowed = [fitz.Rect(287, 0, 532, 216),
               fitz.Rect(0, 504, 125, 620)]
    redactions = []
    for block in p.get_text('dict')['blocks']:
        if block.get('type') != 0:
            continue
        for line in block.get('lines', []):
            for span in line.get('spans', []):
                rect = fitz.Rect(span['bbox'])
                txt = span['text'].strip()
                keep = any(region.contains(rect) for region in allowed)
                if not keep or (txt in list('abcdefghijkl') and span['size'] >= 7):
                    redactions.append(rect)
    for rect in redactions:
        p.add_redact_annot(rect, fill=False, cross_out=False)
    try:
        p.apply_redactions(images=0, graphics=0)
    except TypeError:
        p.apply_redactions(images=0)
    return d


def _fig5_place(p, letter, title, box, doc, clip=None, content_height=None,
                title_font='hebo', title_file=None):
    x, y, w, h = box
    p.insert_text((x, y + 13), letter, fontsize=14, fontname='hebo')
    title_kwargs = dict(fontsize=9.8, fontname=title_font,
                         color=(0.15, 0.18, 0.22))
    if title_file is not None:
        title_kwargs['fontfile'] = title_file
    p.insert_text((x + 18, y + 12), title, **title_kwargs)
    source_rect = fitz.Rect(clip) if clip is not None else doc[0].rect
    max_height = content_height if content_height is not None else h - 23
    scale = min(w / source_rect.width, max_height / source_rect.height)
    dest = fitz.Rect(x, y + 23, x + source_rect.width * scale,
                     y + 23 + source_rect.height * scale)
    p.show_pdf_page(dest, doc, 0, clip=source_rect, keep_proportion=True)
    return dest, scale


def _fig5_restore_negative_axis(p, dest, scale):
    tick_font = fitz.Font(fontfile=FIG5_SANS)
    ticks = [('-1.0', 52.3), ('-0.5', 85.8), ('0.0', 119.3)]
    baseline = dest.y1 + 10.0
    for label, source_x in ticks:
        width = tick_font.text_length(label, fontsize=6.8)
        p.insert_text((dest.x0 + source_x * scale - width / 2.0, baseline),
                      label, fontsize=6.8, fontname='fig5sans',
                      fontfile=FIG5_SANS, color=(0.2, 0.2, 0.2))
    axis = 'Spatial co-org log2 g (25 µm)'
    axis_width = tick_font.text_length(axis, fontsize=7.0)
    p.insert_text((dest.x0 + (dest.width - axis_width) / 2.0, dest.y1 + 21.0),
                  axis, fontsize=7.0, fontname='fig5sans',
                  fontfile=FIG5_SANS, color=(0.2, 0.2, 0.2))


def build_fig5():
    old = _fig5_old_source()
    same_bin = fitz.open(PAN / 'fig5_true_donor_evidence/Fig5_P16_same_bin.pdf')
    distance = fitz.open(PAN / 'fig5_true_donor_evidence/Fig5_P16_distance_support.pdf')
    d, p = newpage(547, 820,
                   'Figure 5 | Local program relationships and real-donor sensitivity')
    _fig5_place(p, 'a', 'Relationship to gene-loading similarity',
                (10, 42, 282, 230), old, clip=fitz.Rect(287, 0, 532, 216))
    b_dest, b_scale = _fig5_place(
        p, 'b', 'Negative relationships', (303, 42, 234, 230), old,
        clip=fitz.Rect(0, 504, 125, 620), content_height=185)
    _fig5_restore_negative_axis(p, b_dest, b_scale)
    _fig5_place(p, 'c', 'P16: same-bin effects and whole-donor omission',
                (10, 287, 527, 247), same_bin)
    _fig5_place(p, 'd', 'P16: separate 0–500 µm effects and whole-donor omission',
                (10, 549, 527, 247), distance, title_font='fig5bold',
                title_file=FIG5_BOLD)
    p.insert_text((12, 811),
                  'Donor coverage: 1 / 10 / 13 / 16 / 4 sections. LOO ranges describe influence, not confidence intervals.',
                  fontsize=7.5, fontname='helv', color=(0.2, 0.2, 0.2))
    finish(d, 'Fig5')
    old.close()
    same_bin.close()
    distance.close()



def build_fig6():
    import math

    source_doc = fitz.open(SRC / 'Fig6.pdf')
    source_page = source_doc[0]
    a_full = fitz.Rect(0, 0, 253, 247)
    b_full = fitz.Rect(0, 325, 247, 505)

    def collect_spans(clip):
        found = []
        for block in source_page.get_text('dict')['blocks']:
            if block.get('type') != 0:
                continue
            for line in block.get('lines', []):
                for span in line.get('spans', []):
                    if fitz.Rect(span['bbox']).intersects(clip):
                        item = dict(span)
                        item['dir'] = line['dir']
                        found.append(item)
        return found

    a_spans = collect_spans(a_full)
    b_spans = collect_spans(b_full)
    a_labels = sorted(
        [s for s in a_spans
         if s['bbox'][0] < 76 and 18 < s['origin'][1] < 238
         and s['text'].strip().startswith('P')],
        key=lambda s: s['origin'][1])
    b_labels = sorted(
        [s for s in b_spans
         if s['bbox'][0] < 72 and s['origin'][1] > 338
         and s['text'].strip().startswith('P')],
        key=lambda s: s['origin'][1])
    a_axis = sorted(
        [s for s in a_spans
         if s['bbox'][1] > 238 and 75 < s['origin'][0] < 225],
        key=lambda s: s['origin'][0])
    b_axis = sorted(
        [s for s in b_spans
         if s['bbox'][1] > 460 and 55 < s['origin'][0] < 240],
        key=lambda s: s['origin'][0])
    a_legend = [s for s in a_spans
                if s['bbox'][0] >= 225 and s['bbox'][1] < 200]
    b_legend = [s for s in b_spans
                if s['bbox'][0] >= 220 and s['bbox'][1] < 460]

    redactions = []
    for block in source_page.get_text('dict')['blocks']:
        if block.get('type') != 0:
            continue
        for line in block.get('lines', []):
            for span in line.get('spans', []):
                rect = fitz.Rect(span['bbox'])
                if rect.intersects(a_full) or rect.intersects(b_full):
                    redactions.append(rect)
                if span['text'].strip() == 'V1' and rect.intersects(
                        fitz.Rect(220, 315, 245, 335)):
                    redactions.append(rect)
    for rect in redactions:
        source_page.add_redact_annot(rect, fill=False, cross_out=False)
    try:
        source_page.apply_redactions(images=0, graphics=0)
    except TypeError:
        source_page.apply_redactions(images=0)

    d = fitz.open()
    p = d.new_page(width=547, height=630)
    text_color = (0.08, 0.08, 0.08)
    label_font = fitz.Font('helv')

    def mapped_point(span, clip, x0, y0, scale):
        return fitz.Point(
            x0 + (span['origin'][0] - clip.x0) * scale,
            y0 + (span['origin'][1] - clip.y0) * scale,
        )

    def put_rotated(span, clip, x0, y0, scale, fontsize):
        text = span['text'].strip()
        if text == '−log10 FDR':
            text = '-log10 FDR'
        if not text:
            return
        pt = mapped_point(span, clip, x0, y0, scale)
        angle = math.degrees(math.atan2(-span['dir'][1], span['dir'][0]))
        p.insert_text(pt, text, fontsize=fontsize, fontname='helv',
                      morph=(pt, fitz.Matrix(angle)), color=text_color)

    # a: source disease matrix and keys remain vector marks; all visible
    # labels are reflowed from the same source into a 7 pt label gutter.
    a_graph = fitz.Rect(75, 15, 245, 238)
    a_x, a_y, a_scale = 284.0, 35.0, 1.35
    a_dest = fitz.Rect(a_x, a_y,
                       a_x + a_graph.width * a_scale,
                       a_y + a_graph.height * a_scale)
    p.show_pdf_page(a_dest, source_doc, 0, clip=a_graph,
                    keep_proportion=True, overlay=True)
    p.insert_text((8, 24), 'a', fontsize=10, fontname='hebo', color=text_color)
    p.insert_text((25, 23), 'Brain-disease enrichment of cortical programs',
                  fontsize=8, fontname='hebo', color=(0.15, 0.18, 0.22))
    for span in a_labels:
        text = span['text'].strip()
        size = 6.5
        baseline = a_y + (span['origin'][1] - a_graph.y0) * a_scale
        width = label_font.text_length(text, fontsize=size)
        p.insert_text((a_x - 9.0 - width, baseline), text,
                      fontsize=size, fontname='helv', color=text_color)
    for span in a_axis + a_legend:
        put_rotated(span, a_graph, a_x, a_y, a_scale, 7.0)

    # b: independent aging matrix and independent FDR/OR keys; the
    # source abbreviations are replaced by source-bound 7 pt text.
    b_graph = fitz.Rect(55, 325, 247, 465)
    b_x, b_y, b_scale = 300.0, 405.0, 1.20
    b_dest = fitz.Rect(b_x, b_y,
                       b_x + b_graph.width * b_scale,
                       b_y + b_graph.height * b_scale)
    p.show_pdf_page(b_dest, source_doc, 0, clip=b_graph,
                    keep_proportion=True, overlay=True)
    p.insert_text((8, 402), 'b', fontsize=10, fontname='hebo', color=text_color)
    p.insert_text((25, 401), 'Aging gene-set enrichment',
                  fontsize=8, fontname='hebo', color=(0.15, 0.18, 0.22))
    for span in b_labels:
        text = span['text'].strip()
        size = 7.0
        baseline = b_y + (span['origin'][1] - b_graph.y0) * b_scale
        width = label_font.text_length(text, fontsize=size)
        p.insert_text((b_x - 9.0 - width, baseline), text,
                      fontsize=size, fontname='helv', color=text_color)
    for span in b_axis + b_legend:
        put_rotated(span, b_graph, b_x, b_y, b_scale, 7.0)

    finish(d, 'Fig6')
    source_doc.close()

def build_fig4_native():
    """Compose Fig4 directly at final width from the locked native a-d panels."""
    d=fitz.open()
    p=d.new_page(width=547,height=838)
    p.insert_text((12,18),
                  'Figure 4 | Subclass-program enrichment and separation in local tissue',
                  fontsize=9.5,fontname='hebo',color=(.08,.12,.19))

    native_panels=[
        ('a',PAN/'fig4a_label_repair/Fig4a.pdf',12,28,523,300),
        ('b',PAN/'fig4a_label_repair/Fig4b.pdf',12,340,255,215),
        ('c',PAN/'fig4a_label_repair/Fig4c.pdf',280,340,255,215),
        ('d',PAN/'fig4a_label_repair/Fig4d.pdf',12,568,255,240),
    ]
    for letter,path,x,y,w,h in native_panels:
        panel=fitz.open(path)
        p.show_pdf_page(fitz.Rect(x,y,x+w,y+h),panel,0,
                        keep_proportion=True,overlay=True)
        p.insert_text((x-9,y+13),letter,fontsize=9,fontname='hebo')
        panel.close()

    # e: use the safe old4g crop and reflow only its complete Upper-IT title.
    # No old4e/old4f content is used and the two maps/legend remain native.
    e_x,e_y=280.0,568.0
    p.insert_text((e_x-9,e_y+13),'e',fontsize=9,fontname='hebo')
    p.insert_text((e_x+9,e_y+13),
                  'V1: OLIGO–P41 enrichment /',
                  fontsize=7,fontname='helv',color=(.15,.18,.22))
    p.insert_text((e_x+9,e_y+22),
                  'upper-IT–P41 separation',
                  fontsize=7,fontname='helv',color=(.15,.18,.22))
    old=fitz.open(SRC/'Fig4.pdf')
    old_page=old[0]
    for block in old_page.get_text('dict')['blocks']:
        if block.get('type')!=0:
            continue
        for line in block.get('lines',[]):
            for span in line.get('spans',[]):
                if span['text'].strip()=='Upper-IT weight':
                    old_page.add_redact_annot(fitz.Rect(span['bbox']),
                                              fill=(1,1,1),cross_out=False)
    try:
        old_page.apply_redactions(images=0,graphics=0)
    except TypeError:
        old_page.apply_redactions(images=0)
    clip=fitz.Rect(312,430,536,531)
    scale=min(255/clip.width,205/clip.height)
    p.show_pdf_page(fitz.Rect(e_x,e_y+35,
                              e_x+clip.width*scale,
                              e_y+35+clip.height*scale),
                    old,0,clip=clip,keep_proportion=True,overlay=True)
    upper_x=e_x+5
    upper_y=e_y+35+(528.8-430)*scale
    p.insert_text((upper_x,upper_y),'Upper-IT weight',
                  fontsize=7,fontname='helv',color=(.10,.10,.10))
    old.close()

    p.insert_text((12,819),
                  'Same-bin marks, not cell contact. The 13-effect and 10-sensitivity sets are distinct;',
                  fontsize=7,fontname='helv',color=(.25,.28,.32))
    p.insert_text((12,831),'P7 weakening is retained.',
                  fontsize=7,fontname='helv',color=(.25,.28,.32))
    finish(d,'Fig4')

def build_s17():
    d,p=newpage(1010,1480,'Figure S17 | Supporting external disease and aging annotations')
    specs=[('a','6b',(15,42,475,305)),('b','6c',(515,42,478,245)),('c','6d',(515,307,478,220)),('d','6g',(15,375,475,320)),('e','6h',(15,735,475,300)),('f','6i',(515,735,478,300)),('g','6e',(15,1080,475,360)),('h','6j',(515,1080,478,360))]
    for l,k,box in specs:place(p,l,'Existing external annotation',box,path=PAN/'fig4a_label_repair'/('FigS17'+l+'_clean.pdf'))
    finish(d,'FigS17',True)

def affected_supplements():
    d,p=newpage(760,1010,'Figure S13 | Additional subclass-program tissue examples')
    place(p,'a','Original three-region maps',(15,42,728,640),path=PAN/'fig4a_label_repair/FigS13a_clean.pdf')
    place(p,'b','P47 / upper-IT examples',(15,714,350,275),path=PAN/'fig4a_label_repair/FigS13b_clean.pdf')
    place(p,'c','OPC / P8 example',(388,714,355,275),path=PAN/'fig4a_label_repair/FigS13c_clean.pdf')
    finish(d,'FigS13',True)
    d,p=newpage(1000,1080,'Figure S14 | Regional relationships and section sign consistency')
    place(p,'a','Regional summaries',(15,42,296,425),key='4j');place(p,'b','Regional clustering',(330,42,306,425),key='4k');place(p,'c','Regional effect ranges',(660,42,322,425),key='4l')
    place(p,'d','Existing section-level effect and sign consistency',(15,515,950,525),path=PAN/'fig4a_label_repair/Fig4e.pdf');finish(d,'FigS14',True)
    d,p=newpage(930,1620,'Figure S15 | Supporting program-program spatial views')
    place(p,'a','Myelin and neuropil fields',(15,42,440,585),key='5c');place(p,'b','Original 331-edge descriptive network',(478,42,435,355),key='5d')
    place(p,'c','P36–P26: one V1 section',(478,430,435,397),key='5e');place(p,'d','Same single-section bin relation',(15,665,440,435),key='5f')
    place(p,'e','Complete original same-bin matrix',(15,1130,895,460),key='5a');finish(d,'FigS15',True)
    d,p=newpage(940,1300,'Figure S16 | Regional relationships and original section consistency')
    place(p,'a','Representative regional effects',(15,42,909,305),key='5i');place(p,'b','Regional clustering',(15,375,443,510),key='5j');place(p,'c','Regional ranges',(480,375,445,510),key='5k')
    place(p,'d','Section sign consistency; no old donor-LOO panel',(15,920,905,345),path=PAN/'fig4a_label_repair/Fig5_section_hist_clean.pdf');finish(d,'FigS16',True)
    build_s17()
    combine_supplements()

def combine_supplements():
    combined=fitz.open()
    for i in range(1,20):
        z=fitz.open(OUT/'supplementary_figures_pdf'/('FigS'+str(i)+'.pdf'));combined.insert_pdf(z)
    combined.save(OUT/'Supplementary_Figures.pdf',garbage=1,deflate=True)
    print('Produced',OUT/'Supplementary_Figures.pdf',flush=True)


def repair_fig2_subclass_typography():
    """Repair only added text after the retained blocks have already been reordered."""
    main_path = OUT / "main_figures_pdf" / "Fig2.pdf"
    s11_path = OUT / "supplementary_figures_pdf" / "FigS11.pdf"
    bundle_path = OUT / "Supplementary_Figures.pdf"
    old_main = fitz.open(stream=main_path.read_bytes(), filetype="pdf")
    s11 = fitz.open(stream=s11_path.read_bytes(), filetype="pdf")
    old_bundle = fitz.open(stream=bundle_path.read_bytes(), filetype="pdf")
    support = fitz.open(PAN / "fig2_donor_evidence" / "Fig2_within_subclass_regional_support.pdf")
    result = fitz.open()
    p = result.new_page(width=547, height=729)
    p.show_pdf_page(fitz.Rect(0, 0, 547, 539), old_main, 0,
                    clip=fitz.Rect(0, 0, 547, 539), keep_proportion=True)
    p.show_pdf_page(fitz.Rect(8, 539, 539, 729), support, 0, keep_proportion=True)

    def text(page, origin, value, size, bold=False, color=(.13, .15, .17)):
        name = "WithinSansBold" if bold else "WithinSans"
        font = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else FIG5_SANS
        page.insert_font(fontname=name, fontfile=font)
        page.insert_text(origin, value, fontsize=size, fontname=name, color=color, overlay=True)

    spans = [span for block in old_main[0].get_text("dict")["blocks"]
             if block.get("type") == 0 for line in block["lines"] for span in line["spans"]]
    for letter, bounds in [
        ("c", fitz.Rect(8, 250, 27, 272)),
        ("d", fitz.Rect(8, 404, 27, 425)),
        ("e", fitz.Rect(297, 404, 317, 425)),
    ]:
        mark = [span for span in spans if span["text"].strip() == letter
                and span["size"] >= 12 and fitz.Rect(span["bbox"]).intersects(bounds)][-1]
        p.draw_rect(fitz.Rect(mark["bbox"])+(-.5, -.5, .5, .5), color=None, fill=(1, 1, 1), overlay=True)
        text(p, mark["origin"], letter, 13, True, (0, 0, 0))

    s = s11[0]
    for rect in [(348, 1214, 998, 1240), (998, 1232, 1040, 1249),
                 (348, 1467, 1040, 1497), (348, 1510, 1040, 1537)]:
        s.draw_rect(fitz.Rect(rect), color=None, fill=(1, 1, 1), overlay=True)
    text(s, (350, 1233), "h", 14, True, (0, 0, 0))
    text(s, (368, 1232), "Regional program activity (within-program z score)", 9.8, True)
    text(s, (350+(521-8)*675/531, 1245), "z-score", 6.5*675/531)
    text(s, (350, 1480), "Black outline: each region's top positive regional z score.", 8.5)
    text(s, (350, 1492), "Descriptive area profiles; not donor-adjusted significance.", 8.5)
    text(s, (350, 1528), "i", 14, True, (0, 0, 0))
    text(s, (368, 1527), "Original lobe-level rankings", 9.8, True)

    combined = fitz.open()
    combined.insert_pdf(old_bundle, from_page=0, to_page=9)
    combined.insert_pdf(s11)
    combined.insert_pdf(old_bundle, from_page=11, to_page=len(old_bundle)-1)
    s11.save(s11_path, garbage=1, deflate=True)
    s11_png = OUT / "supplementary_figures_png" / "FigS11.png"
    s.get_pixmap(dpi=300, alpha=False).save(s11_png)
    print("Produced", s11_path, s11_png, flush=True)
    combined.save(bundle_path, garbage=1, deflate=True)
    print("Produced", bundle_path, flush=True)
    result.save(main_path, garbage=1, deflate=True)
    main_png = OUT / "main_figures_png" / "Fig2.png"
    p.get_pixmap(dpi=300, alpha=False).save(main_png)
    print("Produced", main_path, main_png, "successful Fig2 production 7", flush=True)
    for document in (combined, result, support, old_bundle, s11, old_main):
        document.close()

def repair_fig2_label_mask():
    """Narrow this producer's d-label mask without redrawing the retained caption."""
    path = OUT / "main_figures_pdf" / "Fig2.pdf"
    document = fitz.open(stream=path.read_bytes(), filetype="pdf")
    page = document[0]
    changed = False
    for xref in page.get_contents():
        content = document.xref_stream(xref)
        if b"8 304 19 21 re" in content and b"1 1 1 rg" in content:
            document.update_stream(xref, content.replace(b"8 304 19 21 re", b"8 304 11 21 re"))
            changed = True
    if not changed:
        raise RuntimeError("The task-owned d-label mask was not found in the current Fig2.")
    document.save(path, garbage=1, deflate=True)
    png = OUT / "main_figures_png" / "Fig2.png"
    page.get_pixmap(dpi=300, alpha=False).save(png)
    print("Produced", path, png, "successful Fig2 production 8", flush=True)
    document.close()

def refresh_fig2_support_panel():
    """Replace only the existing f footprint; retain all a-e native content."""
    path = OUT / "main_figures_pdf" / "Fig2.pdf"
    original = fitz.open(stream=path.read_bytes(), filetype="pdf")
    panel = fitz.open(PAN / "fig2_donor_evidence" / "Fig2_within_subclass_regional_support.pdf")
    result = fitz.open()
    page = result.new_page(width=547, height=729)
    page.show_pdf_page(fitz.Rect(0, 0, 547, 539), original, 0,
                       clip=fitz.Rect(0, 0, 547, 539), keep_proportion=True)
    page.show_pdf_page(fitz.Rect(8, 539, 539, 729), panel, 0, keep_proportion=True)
    finish_fig2(result)
    print("Successful Fig2 production 9: f palette only", flush=True)
    result.close()
    panel.close()
    original.close()


def gigascience_labels():
    """Edit retained figure objects only; never run the legacy composition route."""
    import io, math, re, zlib, json, base64, sys
    from PIL import Image, ImageDraw
    baselines = json.load(sys.stdin) if "--baseline-stdin" in sys.argv else {}
    baseline_supp = fitz.open(stream=base64.b64decode(baselines["Supplementary_Figures"]), filetype="pdf") if "Supplementary_Figures" in baselines else None

    def spans(page):
        return [dict(s, direction=line['dir'])
                for block in page.get_text('dict', flags=fitz.TEXTFLAGS_DICT & ~fitz.TEXT_PRESERVE_IMAGES)['blocks'] if block.get('type') == 0
                for line in block['lines'] for s in line['spans']]

    def color(value):
        return tuple(((value >> shift) & 255) / 255 for shift in (16, 8, 0))

    def replace_text(page, changes):
        for span, text in changes:
            page.add_redact_annot(fitz.Rect(span['bbox']), fill=None, cross_out=False)
        if changes:
            page.apply_redactions(images=0, graphics=0, text=0)
        for span, text in changes:
            if not text:
                continue
            direction = span['direction']
            rotation = 90 if direction[1] < -.5 else 270 if direction[1] > .5 else 0
            font = 'hebo' if 'Bold' in span['font'] else 'helv'
            page.insert_text(span['origin'], text, fontsize=span['size'], fontname=font,
                             rotate=rotation, color=color(span['color']), overlay=True)

    def repair_fig4_samebin(page):
        current = spans(page)
        edits = []
        additions = []
        for span in current:
            text = span['text']
            if 'near-ring spatial markcorr' in text:
                edits.append(span)
                additions.append((span['origin'], span['size'], True,
                    'Cell type x program: same-bin g([0,25) µm); dot = retained after shared-gene removal'))
            elif text.strip() == '(25 um ring)':
                edits.append(span)
                additions.append((span['origin'], span['size'], True, '([0,25) µm)'))
        if any('r=25' in s['text'] and 390 < s['bbox'][1] < 405 for s in current):
            axis = [s for s in current if 400 < s['bbox'][0] < 475 and 390 < s['bbox'][1] < 405]
            edits.extend(axis)
            label = 'log2 median g ([0,25) µm)'
            size = 6.3047
            center = (405.093689 + 469.535706) / 2
            origin = (center-fitz.get_text_length(label, fontname='helv', fontsize=size)/2, 400.180206)
            additions.append((origin, size, False, label))
        if not edits:
            return
        for span in edits:
            page.add_redact_annot(fitz.Rect(span['bbox']), fill=None, cross_out=False)
        page.apply_redactions(images=0, graphics=0, text=0)
        page.insert_font(fontname='SameBinSans', fontbuffer=fitz.Font('helv').buffer)
        page.insert_font(fontname='SameBinSansBold', fontbuffer=fitz.Font('hebo').buffer)
        for origin,size,bold,label in additions:
            page.insert_text(origin,label,fontsize=size,
                fontname='SameBinSansBold' if bold else 'SameBinSans',color=(0,0,0))

    def repair_fig4_colorbar_block(page):
        # Reuse the original Type0 font and its CID map directly. A subset font
        # cannot safely be reinterpreted as a fresh Unicode font.
        original=next(f for f in page.get_fonts(full=True)
                      if 'LiberationSans-Bold' in f[3] and f[4]=='C2_1')
        doc=page.parent
        unicode_ref=int(doc.xref_get_key(original[0],'ToUnicode')[1].split()[0])
        cmap=doc.xref_stream(unicode_ref).decode('ascii')
        block=re.search(r'beginbfchar(.*?)endbfchar',cmap,re.S).group(1)
        reverse={bytes.fromhex(value).decode('utf-16-be'):code
                 for code,value in re.findall(r'<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>',block)}
        page.add_redact_annot(fitz.Rect(455,214,539,240),fill=None,cross_out=False)
        page.apply_redactions(images=0,graphics=0,text=0)
        doc.xref_set_key(page.xref,'Resources/Font/ColorbarSourceBold',str(original[0])+' 0 R')
        page.insert_text((477.2344055,223),'([0,25) µm)',fontsize=4.9814,
                         fontname='SameBinSansBold',color=(0,0,0))
        commands=[]
        for text,top in [('log2 g',223),('red=coloc / blue=avoid',235)]:
            encoded=''.join(reverse[c] for c in text)
            commands.append('q 0 g BT /ColorbarSourceBold 4.9814 Tf 1 0 0 1 457.2344055 '
                +str(page.rect.height-top)+' Tm <'+encoded+'> Tj ET Q')
        xref=doc.get_new_xref();doc.update_object(xref,'<<>>')
        doc.update_stream(xref,'\n'.join(commands).encode('ascii'))
        contents=page.get_contents()+[xref]
        doc.xref_set_key(page.xref,'Contents','['+' '.join(str(n)+' 0 R' for n in contents)+']')

    def repair_s11_panel_edges(doc):
        # The lower panels reuse the original a labels through one shared Form.
        # Complete a's label edge with only its original text operators, never
        # expose the adjacent ranking title or density curve under that edge.
        items=doc[0].get_xobjects()
        outer=next(x[0] for x in items if x[2]==0 and x[1]=='fzFrm0')
        owner=next(x[0] for x in items if x[2]==outer and x[1]=='fullpage')
        a_ref=next(x[0] for x in items if x[2]==owner and x[1]=='fzFrm0')
        source_ref=next(x[0] for x in items if x[2]==a_ref and x[1]=='fullpage')
        raw=doc.xref_stream(source_ref)
        marker=b'0 5 -5 0 49.5313 450.9561 Tm'
        begin=raw.index(marker)
        after=begin+len(marker)
        boundary=re.search(rb'(?:[-+\d.]+(?:\s+[-+\d.]+){5}\s+Tm)|\bET\b',raw[after:])
        label_ops=raw[begin:after+boundary.start()]
        kind,reference=doc.xref_get_key(a_ref,'Resources/XObject/RotatedProgramLabels')
        if kind=='xref':
            label_ref=int(reference.split()[0])
        else:
            label_ref=doc.get_new_xref()
            doc.update_object(label_ref,'<< /Type /XObject /Subtype /Form /BBox [0 0 547 661] >>')
            doc.xref_set_key(label_ref,'Resources',doc.xref_get_key(source_ref,'Resources')[1])
            doc.xref_set_key(a_ref,'Resources/XObject/RotatedProgramLabels',str(label_ref)+' 0 R')
        doc.update_stream(label_ref,b'q 0 g BT /TT0 1 Tf 0 Tc 0 Tw '+label_ops+b' ET Q')
        doc.update_stream(a_ref,
            b'q 0 450 547 211 re W n /fullpage Do Q '
            b'q 0 446.5 547 3.5 re W n /RotatedProgramLabels Do Q')
        for name,index,value in [('fzFrm0',1,446.5),('fzFrm1',3,447.5),
                                  ('fzFrm2',3,447.5),('fzFrm3',3,447.5),('fzFrm4',3,446.5)]:
            xref=next(x[0] for x in items if x[2]==owner and x[1]==name)
            box=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',doc.xref_get_key(xref,'BBox')[1])]
            box[index]=value
            doc.xref_set_key(xref,'BBox','['+' '.join(map(str,box))+']')
        # Exclude neighboring graphical material at the left edge of d/e.
        # The legitimate vertical axis labels cross these boundaries, so retain
        # their exact original text operators in a narrow label-only strip.
        for name,resource,marker,left,body_left,right,top in [
            ('fzFrm3','RegionAxisLabel',b'0 4.2285 -4.2285 0 233.5439 369.5781 Tm',228.6,235,385,447.5),
            ('fzFrm4','PCAAxisLabel',b'0 4.2988 -4.2988 0 390.3213 365.6064 Tm',385.3,390,547,446.5)]:
            wrapper=next(x[0] for x in items if x[2]==owner and x[1]==name)
            begin=raw.index(marker);after=begin+len(marker)
            boundary=re.search(rb'(?:[-+\d.]+(?:\s+[-+\d.]+){5}\s+Tm)|\bET\b',raw[after:])
            operations=raw[begin:after+boundary.start()]
            kind,reference=doc.xref_get_key(wrapper,'Resources/XObject/'+resource)
            if kind=='xref':axis_ref=int(reference.split()[0])
            else:
                axis_ref=doc.get_new_xref()
                doc.update_object(axis_ref,'<< /Type /XObject /Subtype /Form /BBox [0 0 547 661] >>')
                doc.xref_set_key(axis_ref,'Resources',doc.xref_get_key(source_ref,'Resources')[1])
                doc.xref_set_key(wrapper,'Resources/XObject/'+resource,str(axis_ref)+' 0 R')
            doc.update_stream(axis_ref,b'q 0 g BT /TT0 1 Tf 0 Tc 0 Tw '+operations+b' ET Q')
            commands=(f'q {body_left} 301 {right-body_left} {top-301} re W n /fullpage Do Q '
                      f'q {left} 301 {body_left-left} {top-301} re W n /{resource} Do Q')
            doc.update_stream(wrapper,commands.encode('ascii'))
            doc.xref_set_key(wrapper,'BBox',f'[{left} 301 {right} {top}]')

    def repair_s12_colorbar_labels(page):
        current=spans(page)
        targets=[s for s in current if s['bbox'][0]>440 and 818<s['bbox'][1]<845
                 and s['text'] in ('median z lo','0','hi')]
        if not any(s['text']=='median z lo' for s in targets):
            return
        ticks=[w for w in page.get_text('words') if w[0]>440 and 818<w[1]<845 and w[4] in ('lo','0','hi')]
        for span in targets:
            page.add_redact_annot(fitz.Rect(span['bbox']),fill=None,cross_out=False)
        page.apply_redactions(images=0,graphics=0,text=0)
        page.insert_font(fontname='MiniMapSans',fontbuffer=fitz.Font('helv').buffer)
        size=8.63386058807373
        for word in ticks:
            page.insert_text((word[0],847.5),word[4],fontsize=size,fontname='MiniMapSans',color=(0,0,0))
        label='median z'
        center=(507.674316+804.012939)/2
        page.insert_text((center-fitz.get_text_length(label,fontname='helv',fontsize=size)/2,863),
                         label,fontsize=size,fontname='MiniMapSans',color=(0,0,0))

    def repair_s10_network_edge(doc,page):
        outer=next(x[0] for x in page.get_xobjects() if x[2]==0 and x[1] in ('Fm3','fzFrm2'))
        matrix=fitz.Matrix(*[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',doc.xref_get_key(outer,'Matrix')[1])])
        page=doc.reload_page(page)
        # The uncropped original Form places the green circle's left stroke at
        # x=197.136 and the blue circle's at x=197.874. The x=196 viewport leaves
        # a margin for both. Remove only the old residual "Ne" text immediately
        # left of the already repaired full label; graphics and edges stay intact.
        old_ne=fitz.Rect(189.13299560546875,32.42755126953125,
                        199.11416625976562,43.3134765625)
        old_ne=old_ne*matrix*page.transformation_matrix
        page.add_redact_annot(old_ne,fill=None,cross_out=False)
        page.apply_redactions(images=0,graphics=0,text=0)
        return page

    def change_clip(doc, name, bounds):
        page = doc[0]
        item = next(x for x in page.get_xobjects() if x[2] == 0 and x[1] in (name, 'Fm'+str(int(name[-1])+1)))
        xref = item[0]
        old_box = [float(x) for x in re.findall(r'-?\d+(?:\.\d+)?', doc.xref_get_key(xref, 'BBox')[1])]
        for index, value in bounds.items():
            old_box[index] = value
        doc.xref_set_key(xref, 'BBox', '[' + ' '.join(str(x) for x in old_box) + ']')

    def caption(page, prefix, text):
        blocks = [b for b in page.get_text('blocks') if b[4].replace('\xa0', ' ').strip().startswith(prefix)]
        if not blocks:
            return
        y = min(b[1] for b in blocks)
        for block in blocks:
            page.add_redact_annot(fitz.Rect(block[:4]), fill=None, cross_out=False)
        page.apply_redactions(images=0, graphics=0, text=0)
        page.insert_textbox(fitz.Rect(40, y, page.rect.width - 40, page.rect.height - 35),
                            text, fontsize=8.2, fontname='helv', lineheight=1.15)

    def image_stream(doc, xref, im):
        doc.xref_set_key(xref, 'Width', str(im.width))
        doc.xref_set_key(xref, 'Height', str(im.height))
        doc.xref_set_key(xref, 'BitsPerComponent', '8')
        doc.xref_set_key(xref, 'DecodeParms', 'null')
        doc.update_stream(xref, zlib.compress(im.tobytes(), 9), compress=False)
        doc.xref_set_key(xref, "Filter", "/FlateDecode")

    def compress_images(doc):
        # Match dimensions conservatively, not by image/content identifiers.
        # No image-info hash option or image-identity scan is used.
        placements = {}
        for page in doc:
            for b in page.get_text('dict')['blocks']:
                if b.get('type') != 1:
                    continue
                rect = fitz.Rect(b['bbox'])
                key = (b['width'], b['height'])
                prior = placements.get(key, (0, 0))
                placements[key] = (max(prior[0], rect.width), max(prior[1], rect.height))
        used = set()
        for page in doc:
            for item in page.get_images(full=True):
                xref, mask, width, height, bpc = item[:5]
                if xref in used or bpc != 8:
                    continue
                used.add(xref)
                extent = placements.get((width, height))
                if extent is None:
                    continue
                factor = min(1., max(extent[0] * 300 / 72 / width,
                                     extent[1] * 300 / 72 / height))
                if factor >= .95:
                    continue
                data = doc.extract_image(xref)
                if data.get('colorspace') not in (1, 3):
                    continue
                im = Image.open(io.BytesIO(data['image']))
                if im.mode not in ('RGB', 'L'):
                    continue
                size = (max(1, math.ceil(width * factor)), max(1, math.ceil(height * factor)))
                image_stream(doc, xref, im.resize(size, Image.Resampling.LANCZOS))
                if mask:
                    alpha = Image.open(io.BytesIO(doc.extract_image(mask)['image'])).convert('L')
                    image_stream(doc, mask, alpha.resize(size, Image.Resampling.LANCZOS))

    def save(doc, name, supplement=False):
        pdfdir = OUT / ('supplementary_figures_pdf' if supplement else 'main_figures_pdf')
        pngdir = OUT / ('supplementary_figures_png' if supplement else 'main_figures_png')
        if name != 'FigS10':
            compress_images(doc)
        payload = doc.tobytes(garbage=0, deflate=True, use_objstms=1)
        (pdfdir / (name + '.pdf')).write_bytes(payload)
        rendered = fitz.open(stream=payload, filetype='pdf')
        rendered[0].get_pixmap(dpi=300, alpha=False).save(pngdir / (name + '.png'))
        rendered.close()
        print('Produced', pdfdir / (name + '.pdf'), pngdir / (name + '.png'), flush=True)

    names = ['Fig2', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'FigS2', 'FigS3', 'FigS9',
             'FigS10', 'FigS11', 'FigS12', 'FigS15', 'FigS16']
    selection = next((a.split('=',1)[1] for a in sys.argv if a.startswith('--figures=')), None)
    if selection is not None:
        names = selection.split(',')
    for name in names:
        supplement = name.startswith('FigS')
        directory = 'supplementary_figures_pdf' if supplement else 'main_figures_pdf'
        if name in baselines:
            doc = fitz.open(stream=base64.b64decode(baselines[name]), filetype='pdf')
        elif supplement and baseline_supp is not None:
            doc = fitz.open()
            number = int(name[4:]) - 1
            doc.insert_pdf(baseline_supp, from_page=number, to_page=number)
        else:
            doc = fitz.open(OUT / directory / (name + '.pdf'))
        page = doc[0]
        if name == 'FigS11':
            # Preserve existing text coordinates and vector paths. The original
            # panel uses split WinAnsi Tj strings; change only these label strings.
            for xref in range(1, doc.xref_length()):
                if not doc.xref_is_stream(xref) or doc.xref_get_key(xref, 'Subtype')[1] != '/Form':
                    continue
                raw = doc.xref_stream(xref)
                if b'P42' not in raw and b'5034322a' not in raw:
                    continue
                changed = raw.replace(b'(P42*)Tj', b'(P42)Tj').replace(b'<5034322a>', b'<503432>')
                changed = re.sub(rb'\(P42 Pos\. reg\. va\)Tj(\s*)\(s\)Tj(\s*7\.891 0 Td\s*)(?:<6f852a>|\(o\\205\*\))Tj',
                                 rb'(P42 NDNF synaptic)Tj\1()Tj\2()Tj', changed)
                changed = re.sub(rb'\(P42 Pos\. reg\. vasoconstrictio\)Tj(\s*)\(n\)Tj(\s*13\.503 0 Td\s*)\(\*\)Tj',
                                 rb'(P42 NDNF synaptic signaling)Tj\1()Tj\2()Tj', changed)
                if changed != raw:
                    doc.update_stream(xref, changed)
            repair_s11_panel_edges(doc)
            page = doc.reload_page(page)
        if name == 'FigS10':
            change_clip(doc, 'fzFrm0', {0: 360, 1: 188})
            change_clip(doc, 'fzFrm1', {2: 197})
            change_clip(doc, 'fzFrm2', {0: 196, 3: 186})
            page=repair_s10_network_edge(doc,page)
        elif name == 'FigS15':
            change_clip(doc, 'fzFrm0', {2: 262})
            change_clip(doc, 'fzFrm1', {0: 284})
            change_clip(doc, 'fzFrm2', {0: 267})
            change_clip(doc, 'fzFrm3', {0: 430})
        elif name == 'FigS16':
            change_clip(doc, 'fzFrm1', {2: 301, 3: 791.136 - 629})
            change_clip(doc, 'fzFrm2', {0: 304.5})
        elif name == 'FigS12':
            change_clip(doc, 'fzFrm2', {3: 531})

        changes = []
        for span in spans(page):
            text = span['text']
            new = text
            if name == 'FigS10' and text.strip() in ('loading', 'L2', 'Ne'):
                new = ''
            if 'P42' in text:
                if 'vaso' in text:
                    new = 'P42 NDNF synaptic signaling' if '…' not in text else 'P42 NDNF synaptic...'
                elif text.strip() == 'P42*':
                    new = 'P42'
            new = new.replace('r = 25 µm', '[0,25) µm').replace('(25 µm)', '([0,25) µm)')
            new = new.replace('dot = survives leave-one-gene-out', 'dot = retained after shared-gene removal')
            new = new.replace('Headline niches (LOGO-surviving', 'Shared-gene-retained niches')
            new = new.replace('LOGO gene-removal retention', 'Shared-gene removal retention')
            if name == 'Fig5':
                new = new.replace('donor-LOO stability', 'legacy-group omission (not donors)')
                new = new.replace('bar=donor-LOO min-max', 'bar=legacy-group range')
            if name == 'FigS16' and span['bbox'][1] < 500 and (
                    'donor-LOO' in text or 'Spatial co-org log2' in text):
                new = ''
            if name == 'FigS2':
                if 'Selected cNMF diagnostic K values from the full integer' in text:
                    new = 'Selected cNMF rank diagnostics for the K=60 analysis resolution'
                elif '100 replicates per K' in text:
                    new = 'Formal backbone: 20 starts; independent torch-MU diagnostics: 100 starts.'
                elif 'Criteria use the full integer' in text:
                    new = '- Selected grids: IC, K=20-200; silhouette, K=30-200.'
            if name == 'FigS15' and ('Exemplar co-orgs' in text or 'sufficient spatial coverage' in text):
                new = ''
            if name == 'FigS15' and span['bbox'][1] < 100 and text.startswith('Neurofilament cytoskeleton') and text != 'Neurofilament cytoskeleton':
                new = ''
            if name == 'FigS16' and span['bbox'][0] < 25 and 130 < span['bbox'][1] < 250:
                new = ''
            if name == 'FigS16' and 420 < span['bbox'][0] < 505 and 480 < span['bbox'][1] < 525 and ('log2' in text or text.strip() == '2 g'):
                new = ''
            if new != text:
                changes.append((span, new))
        replace_text(page, changes)
        if name == 'Fig4':
            repair_fig4_samebin(page)
            repair_fig4_colorbar_block(page)
        if name == 'FigS12':
            repair_s12_colorbar_labels(page)

        if name == 'FigS15' and not any(s['text'] == '(pan-neuronal, P8)' and s['bbox'][1] < 100 for s in spans(page)):
            page.insert_text((297.5, 72.5), 'Neurofilament cytoskeleton', fontsize=7.97, fontname='helv')
            page.insert_text((297.5, 81.4), '(pan-neuronal, P8)', fontsize=7.97, fontname='helv')
        if name == 'FigS16' and not any(s['text'] == 'Log2 g' and 480 < s['bbox'][1] < 525 for s in spans(page)):
            page.insert_text((437.3, 504.4), 'Log2 g', fontsize=10, fontname='helv')

        if name == 'Fig3' and not any(s['text'].strip() == 'a' and s['bbox'][0] < 20 and s['bbox'][1] < 20 for s in spans(page)):
            page.insert_text((3, 11), 'a', fontsize=12, fontname='hebo')
        if name == 'FigS10' and not any(s['text']=='Neurofilament cytoskeleton (pan-neuronal)' for s in spans(page)):
            broken = [s for s in spans(page) if 'filament cytoskeleton' in s['text'] and s['bbox'][1] > 550]
            replace_text(page, [(s, '') for s in broken])
            page.insert_text((90, 850), 'Neurofilament cytoskeleton (pan-neuronal)', fontsize=12, fontname='helv')
        if name == 'FigS12' and not any(s['text'] == 'Median program standardized activity' for s in spans(page)):
            broken = [s for s in spans(page) if 'median program standardized' in s['text']]
            replace_text(page, [(s, '') for s in broken])
            page.insert_text((111, 407), 'Median program standardized activity', fontsize=9, fontname='helv')
        if name == 'FigS2':
            caption(page, 'Fig. S2.',
                'Fig. S2. Internal rank and solver diagnostics. a, Information criteria for the selected ranks shown. '
                'b, Concordance of regional-variability diagnostics at selected ranks. c, Native cNMF stability '
                'from the 20-initialization backbone. d, Synthesis of the diagnostic results. The fixed K=60 '
                'reference comes from the 20-initialization backbone; independent torch-MU runs used 100 starts. '
                'The retained IC grid spans K=20-200, whereas the silhouette/distinctness grid spans K=30-200; '
                'both are sparse selected grids, not evaluations of every integer K. Six source components were '
                'excluded from K=60, leaving 54 retained programs.')
        if name == 'FigS9' and not any('P42 NDNF interneuron' in s['text'] for s in spans(page)):
            # This existing panel is already a raster. Erase only the old P42
            # label pixels on its white margin, then write the replacement as
            # real PDF text; no plotted cell, coefficient or scatter point changes.
            image_item = page.get_images(full=True)[0]
            xref = image_item[0]
            im = Image.open(io.BytesIO(doc.extract_image(xref)['image'])).convert('RGB')
            pixels = im.load()
            rows = []
            for y in range(721, 746):
                if sum(min(pixels[x, y]) < 190 for x in range(25, 285)) > 3:
                    rows.append(y)
            groups = []
            for y in rows:
                if not groups or y > groups[-1][-1] + 1:
                    groups.append([y])
                else:
                    groups[-1].append(y)
            band = min(groups, key=lambda g: abs((g[0]+g[-1])/2 - 733))
            ImageDraw.Draw(im).rectangle((25, band[0], 285, band[-1]), fill='white')
            doc.xref_set_key(xref, 'ColorSpace', '/DeviceRGB')
            image_stream(doc, xref, im)
            box = next(b['bbox'] for b in page.get_text('dict')['blocks'] if b['type'] == 1)
            pixel_scale = (box[2]-box[0]) / im.width
            label = 'P42 NDNF interneuron synaptic signaling'
            size = 3.5
            right = box[0] + 285 * pixel_scale
            baseline = box[1] + (band[-1] + .4) * pixel_scale
            page.insert_text((right-fitz.get_text_length(label, fontname='helv', fontsize=size), baseline),
                             label, fontsize=size, fontname='helv')
        if name == 'FigS9':
            caption(page, 'Fig. S9.',
                'Fig. S9. Regional program profiles and Neurosynth functional maps. a, Spearman correlations '
                'for 54 programs and 12 terms. Marked cells have q<0.1. b, Default-mode ranking. c, Examples '
                'across the 14 sampled cortical regions. d, Term-map values in the matched Harvard-Oxford ROIs. '
                'Two-sided P values use 10,000 ordinary region-label permutations (seed 0; plus-one correction); '
                'BH correction is applied across all 648 program-term tests. These permutations do not preserve '
                'spatial autocorrelation, so the associations are exploratory. P36-default mode has r=0.8896 '
                'and q=0.0648, not q<0.05.')
        save(doc, name, supplement)
        doc.close()

    if not any(name.startswith('FigS') for name in names):
        return
    combined = fitz.open()
    for i in range(1,20):
        part = fitz.open(OUT/'supplementary_figures_pdf'/('FigS'+str(i)+'.pdf'))
        combined.insert_pdf(part)
        part.close()
    (OUT/'Supplementary_Figures.pdf').write_bytes(combined.tobytes(garbage=0, deflate=True, use_objstms=1))
    print('Produced', OUT/'Supplementary_Figures.pdf', flush=True)


S14_CURRENT_CAPTION = "Fig. S14. Spatial distributions and associations of cortical programs. a, Myelination (P41), oligodendrocyte development (P12*) and pan-neuronal neurofilament cytoskeleton expression (P8) in DLPFC, V1 and M1. b, Network of 331 positive same-bin pairs with full-loading cosine <0.25 and log2 median g >0.32. c, Microglial immune activation (P36) and postsynaptic glutamate-receptor expression (P26*) in one V1 section. d, Their bin-level correlation, r = 0.79. e, The complete 54-program same-bin association matrix. f, P26*–P36 profiles across ten distance rings. Thin lines show the 14 cortical-area summaries; the dark line and shading show log2 of the median g and log2-transformed first and third quartiles across all 44 sections. Area section counts are given in the key. g, Same-bin P26*–P36 effects for five biological donors and five whole-donor omissions. The dashed reference marks the all-section effect; the lower point and range give the all-44-section estimate and the minimum and maximum of the omission estimates. h, Extracellular-matrix organization (P9) and blood-vessel morphogenesis expression (P51) in M1, with the SCT-derived spatial program score, tissue mask and individual display ranges. Each scale bar is 1 mm. i, Negative associations between neuronal and oligodendroglial/myelin programs. Each segment joins zero to its pair's effect, and point color indicates the section same-direction fraction. Asterisks denote lower functional-annotation confidence."


def restore_fig5_supplements():
    """Restore only approved S14 panels and the S15 caption from current pages."""
    import sys, io, json, base64
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import render_fig5_revision as native

    payload = json.load(sys.stdin)
    base14 = fitz.open(stream=base64.b64decode(payload['FigS14']), filetype='pdf')
    base15 = fitz.open(stream=base64.b64decode(payload['FigS15']), filetype='pdf')
    old5 = fitz.open(stream=base64.b64decode(payload['Fig5']), filetype='pdf')
    destination = ROOT / 'analysis/fig5_revision'
    caption_font_path = '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
    caption_font = fitz.Font(fontfile=caption_font_path)

    def remove_caption(doc, prefix):
        page = doc[0]
        blocks = [b for b in page.get_text('blocks')
                  if b[4].replace('\xa0', ' ').strip().startswith(prefix)]
        if not blocks:
            raise ValueError('Required current-page caption is missing: ' + prefix)
        box = fitz.Rect(blocks[0][:4])
        for block in blocks:
            box |= fitz.Rect(block[:4])
            page.add_redact_annot(fitz.Rect(block[:4]), fill=False, cross_out=False)
        page.apply_redactions(images=0, graphics=0, text=0)
        return box

    def caption_lines(text, width, size):
        lines, current = [], ''
        for word in text.split():
            candidate = (current + ' ' + word).strip()
            if current and caption_font.text_length(candidate, fontsize=size) > width:
                lines.append(current)
                current = word
            else:
                current = candidate
        if current:
            lines.append(current)
        return lines

    def write_caption(page, x, top, lines, size):
        page.insert_font(fontname='SupplementCaption', fontfile=caption_font_path)
        step = size * 1.18
        for index, line in enumerate(lines):
            page.insert_text((x, top + size + index * step), line,
                             fontsize=size, fontname='SupplementCaption', color=(.08, .08, .08))

    def save_pair(doc, name):
        pdf = destination / (name + '.pdf')
        png = destination / (name + '.png')
        doc.save(pdf, garbage=4, deflate=True)
        pix = doc[0].get_pixmap(dpi=300, alpha=False)
        pix.save(png)
        print('Produced', pdf, pdf.stat().st_size, 'bytes;', doc[0].rect, flush=True)
        print('Produced', png, png.stat().st_size, 'bytes;',
              pix.width, 'x', pix.height, 'px;', pix.n, 'channels', flush=True)

    cap14 = remove_caption(base14, 'Fig. S14.')
    cap15 = remove_caption(base15, 'Fig. S15.')
    width = float(base14[0].rect.width)
    appendix_height = 850.
    appendix_y = cap14.y0 - 4.

    # These changes are process-local helper geometry, not edits to main Fig5.
    native.W, native.H = width, appendix_height
    native.FS, native.SMALL, native.TITLE, native.TAG = 12., 10., 14., 14.
    native.FINAL_WIDTH_MM = width / 72 * 25.4
    plt.rcParams.update({'font.size':12, 'axes.labelsize':11,
                         'xtick.labelsize':10, 'ytick.labelsize':10})
    ann = {r['new_P']:r for r in native.read(ROOT/'tables/TableS3_program_annotation.tsv')}
    pair = ('P26', 'P36')
    cross = native.CROSS
    with np.load(cross/'markcorr_v2/final/progprog_median_iqr.npz', allow_pickle=True) as data:
        ia = list(data['A_names']).index('program_' + ann[pair[0]]['cnmf_component'])
        ib = list(data['B_names']).index('program_' + ann[pair[1]]['cnmf_component'])
        global_profile = tuple(data[field][ia,ib,:].copy()
                               for field in ['log2_median_g','log2_q1','log2_q3'])
        edges = data['ring_edges_um'].copy()
    with np.load(cross/'markcorr_v2/final/progprog_byarea_median_iqr.npz', allow_pickle=True) as data:
        ia = list(data['A_names']).index('program_' + ann[pair[0]]['cnmf_component'])
        ib = list(data['B_names']).index('program_' + ann[pair[1]]['cnmf_component'])
        area_names = [str(v) for v in data['area_names']]
        area_counts = data['n_chips_area'].tolist()
        area_edges = data['ring_edges_um'][1:].copy()
        regional_profiles = data['log2_median_g'][:,ia,ib,:].copy()
    profiles = {pair:(global_profile, regional_profiles, area_edges)}
    donor_rows = {r['donor']:r for r in native.read(
        ROOT/'analysis/02_true_donor_spatial_support/same_bin_per_donor_effects.tsv')
        if r['mode']=='progprog' and {r['A_label'],r['B_label']}==set(pair)}
    loo_rows = {r['omitted_donor']:r for r in native.read(
        ROOT/'analysis/02_true_donor_spatial_support/same_bin_true_donor_loo.tsv')
        if r['mode']=='progprog' and {r['A_label'],r['B_label']}==set(pair)}
    donors = native.DONORS
    summary = loo_rows[donors[0]]
    full = float(summary['all_sections_log2_median_g'])
    lower = float(summary['loo_min_log2_median_g'])
    upper = float(summary['loo_max_log2_median_g'])

    fig = plt.figure(figsize=(width/72, appendix_height/72), facecolor='none')
    fig.patch.set_alpha(0)
    native.tag(fig,15,20,'f')
    native.tx(fig,33,20,'P26*–P36: ten-ring distance profiles',14,'bold')
    ax = native.distance(fig,65,60,pair,profiles,edges,w=350,h=185)
    ax.xaxis.label.set_size(11)
    ax.yaxis.label.set_size(11)
    ax.tick_params(labelsize=10)
    native.tx(fig,30,282,'Region (n sections); all 44 sections: median and IQR',11)
    handles = [Line2D([],[],color=native.AREA_COLOURS[i],
               linestyle=native.AREA_STYLES[i%len(native.AREA_STYLES)],
               lw=.95,label=f'{name} (n={int(area_counts[i])})')
               for i,name in enumerate(area_names)]
    handles.append(Line2D([],[],color=native.GLOBAL_COLOUR,lw=1.6,label='All 44 (IQR)'))
    fig.legend(handles=handles,loc='upper left',
               bbox_to_anchor=(25/width,1-294/appendix_height),ncol=4,
               frameon=False,fontsize=10,handlelength=1.4,handletextpad=.35,
               columnspacing=.65,labelspacing=.35,borderaxespad=0,borderpad=0)

    native.tag(fig,478,20,'g')
    native.tx(fig,496,20,'P26*–P36: five donors and donor omission',14,'bold')
    native.tx(fig,590,53,'Five biological donors',11,'bold',ha='center')
    native.tx(fig,840,53,'One donor omitted',11,'bold',ha='center')
    axd = native.axat(fig,520,70,150,170)
    axo = native.axat(fig,765,70,150,170)
    for index, donor in enumerate(donors):
        axd.plot(float(donor_rows[donor]['donor_log2_median_g']),index,'o',
                 color=native.DC[index],ms=6,mew=0)
        axo.plot(float(loo_rows[donor]['loo_log2_median_g']),index,'D',
                 color='#333333',ms=4.6,mew=0)
    for panel in [axd,axo]:
        panel.axvline(full,color='#7b8790',lw=.8,ls='--',zorder=0)
        panel.set_xlim(.75,1.35)
        panel.set_ylim(4.6,-.6)
        panel.set_xticks([.8,1.,1.2])
        panel.set_xlabel('log2 median g',fontsize=11,labelpad=4)
        panel.spines[['top','right']].set_visible(False)
        panel.tick_params(labelsize=10,length=2.5,pad=3)
        panel.grid(axis='x',color='#e8ebee',lw=.45)
    axd.set_yticks(range(5),[f"{d.removeprefix('Donor')} (n={donor_rows[d]['n_sections_donor']})"
                           for d in donors])
    axo.set_yticks(range(5),[d.removeprefix('Donor')+' omitted' for d in donors])
    native.tx(fig,480,291,'n = sections; dashed line = all 44 sections',11)
    native.tx(fig,480,319,'All 44 / donor-omission range',11,'bold')
    axr = native.axat(fig,765,297,150,26)
    axr.hlines(0,lower,upper,color='#333333',lw=1)
    axr.vlines([lower,upper],-.13,.13,color='#333333',lw=1)
    axr.plot(full,0,'D',color='#173e5c',ms=5,mew=0)
    axr.set_xlim(.75,1.35);axr.set_ylim(-.4,.4)
    axr.set_xticks([.8,1.,1.2]);axr.set_yticks([])
    axr.spines[['top','right','left']].set_visible(False)
    axr.tick_params(labelsize=10,length=2.5,pad=3)
    native.tx(fig,480,353,
              f'All 44: {full:.3f}; omission range: {lower:.3f}–{upper:.3f}',11)

    native.tag(fig,15,404,'h')
    native.tx(fig,33,404,'M1: ECM and vessel programs',14,'bold')
    for cx,pid in [(122,'P9'),(342,'P51')]:
        colour=native.TYPE_TEXT[native.TYPE_BY_PROGRAM[pid]]
        native.tx(fig,cx,430,native.SHORT[pid],12,'bold',ha='center',color=colour)
        native.tx(fig,cx,447,pid+' · '+ann[pid]['dominant_subclass'],11,
                  ha='center',color=colour)
    # Keep the adopted three-section common physical span, while restoring only M1.
    native.PAIRS=[('P9','P51')]
    spatial,geometry,span=native.load_spatial(ann)
    native.tissue(fig,32,465,'A01186A4','P9',ann,spatial,geometry,span,side=180)
    native.tissue(fig,252,465,'A01186A4','P51',ann,spatial,geometry,span,side=180)
    bar=native.axat(fig,90,692,280,6)
    bar.imshow(np.linspace(0,1,256)[None,:],cmap='viridis',aspect='auto')
    bar.set_axis_off()
    native.tx(fig,90,687,'low',10)
    native.tx(fig,370,687,'high',10,ha='right')
    native.tx(fig,230,718,'Program score',12,ha='center')
    native.tag(fig,478,404,'i')
    native.tx(fig,496,404,'Original negative program-pair effects',14,'bold')

    stream=io.BytesIO()
    fig.savefig(stream,format='pdf',dpi=300,transparent=True,facecolor='none')
    plt.close(fig)
    appended=fitz.open(stream=stream.getvalue(),filetype='pdf')

    # Isolate only the retained original g in memory, preserving its native marks.
    old_page=old5[0]
    # The old adjacent histogram starts beyond the retained zero-axis label.
    clip=fitz.Rect(0,504,124.5,628)
    for block in old_page.get_text('dict')['blocks']:
        for line in block.get('lines',[]):
            for span_text in line.get('spans',[]):
                rect=fitz.Rect(span_text['bbox'])
                if span_text['text'].strip()=='g' and rect.x0<20 and rect.y0<528 and rect.intersects(clip):
                    old_page.add_redact_annot(rect,fill=False,cross_out=False)
    bounds=old_page.rect
    for rect in [fitz.Rect(0,0,bounds.width,clip.y0),
                 fitz.Rect(clip.x1,clip.y0,bounds.width,clip.y1),
                 fitz.Rect(0,clip.y1,bounds.width,bounds.height)]:
        old_page.add_redact_annot(rect,fill=False,cross_out=False)
    old_page.apply_redactions(images=1,graphics=1,text=0)
    old_page.set_cropbox(clip)
    isolated=fitz.open(stream=old5.tobytes(garbage=4,deflate=True),filetype='pdf')

    text14=S14_CURRENT_CAPTION
    font14=15.63025
    lines14=caption_lines(text14,805,font14)
    top14=appendix_y+appendix_height+22
    height14=top14+len(lines14)*font14*1.18+65
    out14=fitz.open()
    page14=out14.new_page(width=width,height=height14)
    page14.show_pdf_page(base14[0].rect,base14,0)
    page14.show_pdf_page(fitz.Rect(0,appendix_y,width,appendix_y+appendix_height),appended,0)
    page14.show_pdf_page(fitz.Rect(495,appendix_y+430,915,appendix_y+831),isolated,0)
    write_caption(page14,62.52,top14,lines14,font14)
    save_pair(out14,'FigS14')

    text15=(
        'Fig. S15. Regional patterns and section-direction distributions of program associations. '
        'a, Area-level P33-P34, P8-P32, P36-P26 and P16-P43 effects. '
        'b, Regional profiles of 40 program pairs with the retained clustering. '
        'c, Regional minima (blue) and maxima (red) for the 22 retained pairs. '
        'd, Across 44 sections, the fraction whose log2(g) sign matches the sign of log2(median g), '
        'for selected pairs among the 54 retained programs. Existing selection used '
        'Stouffer-based BH q <0.05, |median(log2 g)| >0.32, and a fraction >=0.85 of section Z_i '
        'signs matching the sign of median(log2 g). The panels expand Fig. 5f-i.'
    )
    font15=15.79832
    lines15=caption_lines(text15,813.55,font15)
    height15=max(base15[0].rect.height,cap15.y0+len(lines15)*font15*1.18+65)
    out15=fitz.open()
    page15=out15.new_page(width=base15[0].rect.width,height=height15)
    page15.show_pdf_page(base15[0].rect,base15,0)
    write_caption(page15,63.1933,cap15.y0,lines15,font15)
    save_pair(out15,'FigS15')
    print('Retained S14 a-e; appended f=P26-P36 curves, g=true-donor/LOO, '
          'h=M1 P9/P51, i=original Fig5g. S15 a-d unchanged; caption only.',flush=True)
    for doc in [out14,out15,appended,isolated,old5,base14,base15]:
        doc.close()



def repair_s14_caption_only():
    """Synchronize only the retained S14 native caption and program-score label."""
    import sys
    document=fitz.open(stream=sys.stdin.buffer.read(),filetype='pdf')
    page=document[0]
    blocks=page.get_text('dict')['blocks']
    caption_blocks=[b for b in blocks if b.get('type')==0
                    and ' '.join(s['text'] for line in b['lines'] for s in line['spans']).replace('\xa0',' ').startswith('Fig. S14.')]
    if not caption_blocks:
        raise ValueError('Current S14 caption was not found.')
    spans=[s for b in caption_blocks for line in b['lines'] for s in line['spans']]
    labels=[s for b in blocks if b.get('type')==0 for line in b['lines'] for s in line['spans']
            if s['text'].strip() in ('GEP score','Program score')]
    if len(labels)!=1:
        raise ValueError('The unique retained S14 program-score label was not found.')
    label=labels[0]
    text=S14_CURRENT_CAPTION
    first=spans[0]
    x,baseline=first['origin']
    size=first['size']
    font_path='/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
    label_font_path='/usr/share/fonts/truetype/liberation/LiberationSansNarrow-Regular.ttf'
    font=fitz.Font(fontfile=font_path)
    label_font=fitz.Font(fontfile=label_font_path)
    lines=[];line=''
    for word in text.split():
        candidate=(line+' '+word).strip()
        if line and font.text_length(candidate,fontsize=size)>805:
            lines.append(line);line=word
        else:
            line=candidate
    if line:
        lines.append(line)
    step=size*1.18
    # Keep the original font size and use only the existing blank caption band.
    last_baseline=page.rect.height-8-max(0.,-font.descender)*size
    if len(lines)>1:
        step=min(step,(last_baseline-baseline)/(len(lines)-1))
    if step<size*1.05:
        raise ValueError('The current caption band cannot retain the approved font size.')
    for block in caption_blocks:
        page.add_redact_annot(fitz.Rect(block['bbox']),fill=False,cross_out=False)
    page.add_redact_annot(fitz.Rect(label['bbox']),fill=False,cross_out=False)
    page.apply_redactions(images=0,graphics=0,text=0)
    page.insert_font(fontname='S14CaptionRevision',fontfile=font_path)
    for index,line in enumerate(lines):
        page.insert_text((x,baseline+index*step),line,fontsize=size,
                         fontname='S14CaptionRevision',color=(.08,.08,.08))
    page.insert_font(fontname='S14ProgramScore',fontfile=label_font_path)
    label_text='Program score'
    label_x=(label['bbox'][0]+label['bbox'][2])/2-label_font.text_length(label_text,fontsize=label['size'])/2
    label_color=tuple(((label['color']>>shift)&255)/255 for shift in (16,8,0))
    page.insert_text((label_x,label['origin'][1]),label_text,fontsize=label['size'],
                     fontname='S14ProgramScore',color=label_color)
    directory=ROOT/'analysis/fig5_revision'
    pdf=directory/'FigS14.pdf';png=directory/'FigS14.png'
    document.save(pdf,garbage=4,deflate=True)
    # Saving with object compaction can invalidate the pre-save page handle.
    rendered=fitz.open(pdf)
    pix=rendered[0].get_pixmap(dpi=300,alpha=False)
    pix.save(png)
    print('Produced',pdf,pdf.stat().st_size,'bytes;',rendered[0].rect,flush=True)
    print('Produced',png,png.stat().st_size,'bytes;',pix.width,'x',pix.height,'px',flush=True)
    print('S14 caption and program-score label synchronized; native a-i and page geometry retained;',
          'caption lines',len(lines),'font size',size,'line step',step,'S15 unchanged.',flush=True)
    rendered.close()
    document.close()


def clean_fig2_native_content():
    """Remove hidden/out-of-placement content while retaining Fig2 native vector marks."""
    import io
    import pymupdf as fitz
    from collections import defaultdict, deque
    from pypdf import PdfReader
    from pypdf.generic import (ArrayObject, ByteStringObject, ContentStream,
                              DictionaryObject, IndirectObject, NameObject,
                              TextStringObject)

    figure_path=OUT/'main_figures_pdf/Fig2.pdf'
    figure_source=figure_path.read_bytes()
    reader=PdfReader(io.BytesIO(figure_source))
    document=fitz.open(stream=figure_source,filetype='pdf')
    page=reader.pages[0]
    width,height=float(page.mediabox.width),float(page.mediabox.height)
    identity=fitz.Matrix(1,0,0,1,0,0)
    pdfdoc=fitz._as_pdf_document(document)
    mupdf=fitz.mupdf
    counters=defaultdict(int)
    parsed={}
    hebo=next((item[0] for item in document[0].get_fonts(full=True)
               if item[3]=='Helvetica-Bold' and item[4]=='hebo'),None)
    if hebo is None:
        raise ValueError('Fig2 is missing its retained Helvetica-Bold resource.')

    def dictionary(value):
        result=DictionaryObject()
        for key,item in value.items():result[key]=item
        return result

    def original_codes(value):
        # Preserve original encoded glyph bytes, including CID/CFF subsets.
        if isinstance(value,TextStringObject):return ByteStringObject(value.original_bytes)
        if isinstance(value,ArrayObject):return ArrayObject([original_codes(v) for v in value])
        return value

    def content_bytes(operations):
        stream=ContentStream(None,reader)
        stream.operations=[([original_codes(v) for v in args],op) for args,op in operations]
        return stream.get_data()

    def write_dictionary(value,target=document):
        buffer=io.BytesIO();value.write_to_stream(buffer)
        xref=target.get_new_xref()
        target.update_object(xref,buffer.getvalue().decode('latin1'))
        return xref

    def normalized_box(value):
        x0,y0,x1,y1=map(float,value)
        return fitz.Rect(min(x0,x1),min(y0,y1),max(x0,x1),max(y0,y1))

    class ClipFilter(mupdf.PdfSanitizeFilterOptions2):
        def __init__(self,clip):
            super().__init__();self.clip=clip;self.use_virtual_culler()
        def culler(self,ctx,bbox,kind):
            # Do not discard a clipping path or alter a crossing subpath.
            if kind in (0,4,5,6,7):return 0
            rect=fitz.Rect(bbox.x0,bbox.y0,bbox.x1,bbox.y1)
            remove=(rect & self.clip).is_empty
            if remove:counters['outside_glyphs' if kind==8 else 'outside_graphical_subpaths']+=1
            return int(remove)

    def rebuild(operations,resources,world,clip):
        revised=dictionary(resources)
        xobjects=resources.get('/XObject')
        xobjects=xobjects.get_object() if xobjects is not None else DictionaryObject()
        retained=DictionaryObject();result=[];stack=[];marks=[]
        matrix=fitz.Matrix(world);current_clip=fitz.Rect(clip)
        path_bounds=None;pending_clip=False
        def extend(points):
            nonlocal path_bounds
            for point in points:
                point=fitz.Point(*point)*matrix
                if path_bounds is None:path_bounds=[point.x,point.y,point.x,point.y]
                else:
                    path_bounds[0]=min(path_bounds[0],point.x)
                    path_bounds[1]=min(path_bounds[1],point.y)
                    path_bounds[2]=max(path_bounds[2],point.x)
                    path_bounds[3]=max(path_bounds[3],point.y)
        for args,op in operations:
            if op in (b'BDC',b'BMC'):
                optional=bool(args and args[0]=='/OC');marks.append(optional)
                if optional:counters['removed_optional_wrappers']+=1;continue
            elif op==b'EMC' and marks:
                if marks.pop():continue
            if op==b'q':stack.append((fitz.Matrix(matrix),fitz.Rect(current_clip)))
            elif op==b'Q' and stack:matrix,current_clip=stack.pop()
            elif op==b'cm':matrix=fitz.Matrix(*map(float,args))*matrix
            elif op in (b'm',b'l'):extend([(float(args[0]),float(args[1]))])
            elif op==b're':
                x,y,w,h=map(float,args);extend([(x,y),(x+w,y),(x+w,y+h),(x,y+h)])
            elif op in (b'c',b'v',b'y'):
                extend([(float(args[i]),float(args[i+1])) for i in range(0,len(args),2)])
            elif op in (b'W',b'W*'):pending_clip=True
            elif op in (b'n',b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*'):
                if pending_clip and path_bounds is not None:
                    current_clip &= fitz.Rect(path_bounds)
                path_bounds=None;pending_clip=False
            elif op==b'Do' and args[0] in xobjects:
                reference=xobjects.raw_get(args[0]);obj=reference.get_object()
                if obj.get('/Subtype')=='/Form':
                    xref=instance(reference,resources,matrix,current_clip)
                    if xref is None:
                        counters['outside_form_calls']+=1
                        continue
                    name=NameObject('/Native'+str(len(retained)))
                    retained[name]=IndirectObject(xref,0,reader)
                    result.append(([name],op));continue
                retained[args[0]]=reference
            result.append((args,op))
        if xobjects:revised[NameObject('/XObject')]=retained
        needed={args[0] for args,op in result if op==b'Tf'}
        fonts=resources.get('/Font')
        fonts=fonts.get_object() if fonts is not None else DictionaryObject()
        if '/hebo' in needed and '/hebo' not in fonts:
            fonts=dictionary(fonts);fonts[NameObject('/hebo')]=IndirectObject(hebo,0,reader)
            revised[NameObject('/Font')]=fonts;counters['font_bindings']+=1
        properties=resources.get('/Properties')
        if properties is not None:
            kept=DictionaryObject()
            for name,reference in properties.get_object().items():
                if reference.get_object().get('/Type') not in ('/OCG','/OCMD'):kept[name]=reference
            if kept:revised[NameObject('/Properties')]=kept
            else:revised.pop('/Properties',None)
        return result,revised

    def instance(reference,inherited,parent,clip):
        obj=reference.get_object()
        matrix=fitz.Matrix(*map(float,obj.get('/Matrix',(1,0,0,1,0,0))))
        world=matrix*parent
        box=normalized_box(obj['/BBox'])
        effective=clip & (box*world)
        if effective.is_empty:return None
        key=(reference.idnum,reference.generation)
        if key not in parsed:parsed[key]=ContentStream(obj,reader).operations
        resources=obj.get('/Resources',inherited).get_object()
        operations,revised=rebuild(parsed[key],resources,world,effective)
        xref=document.get_new_xref()
        document.update_object(xref,document.xref_object(reference.idnum,compressed=False))
        resource_xref=write_dictionary(revised)
        document.xref_set_key(xref,'Resources',f'{resource_xref} 0 R')
        document.xref_set_key(xref,'BBox','['+' '.join(format(v,'.17g') for v in box)+']')
        if obj.get('/OC') is not None:document.xref_set_key(xref,'OC','null')
        document.update_stream(xref,content_bytes(operations),compress=True)
        document.xref_set_key(xref,'DecodeParms','null')
        # Per-placement filtering; fonts and visible paths stay native.
        options=ClipFilter(effective)
        filters=fitz._make_PdfFilterOptions(recurse=0,instance_forms=0,sanitize=1,sopts=options)
        filtered=mupdf.pdf_filter_xobject_instance(
            mupdf.pdf_new_indirect(pdfdoc,xref,0),
            mupdf.pdf_new_indirect(pdfdoc,resource_xref,0),
            mupdf.FzMatrix(*parent),filters,mupdf.PdfCycleList())
        counters['retained_form_instances']+=1
        return mupdf.pdf_to_num(filtered)

    root_operations,root_resources=rebuild(
        ContentStream(page['/Contents'],reader).operations,page['/Resources'].get_object(),
        identity,fitz.Rect(0,0,width,height))
    contents=document.get_new_xref();document.update_object(contents,'<<>>')
    document.update_stream(contents,content_bytes(root_operations),compress=True)
    resource_xref=write_dictionary(root_resources)
    document.xref_set_key(document[0].xref,'Contents',f'{contents} 0 R')
    document.xref_set_key(document[0].xref,'Resources',f'{resource_xref} 0 R')

    # Remove earlier glyphs completely covered by later opaque white rectangles.
    # These are observed old labels, not all text or all repeated marks.
    native_page=document[0]
    masks=[];clips={};groups={}
    for drawing in native_page.get_drawings(extended=True):
        level=drawing['level']
        clips={k:v for k,v in clips.items() if k<level}
        groups={k:v for k,v in groups.items() if k<level}
        if drawing['type']=='clip':clips[level]=fitz.Rect(drawing['scissor']);continue
        if drawing['type']=='group':groups[level]=drawing.get('opacity',1);continue
        fill=drawing.get('fill')
        if not fill or min(fill)<.999999 or drawing.get('fill_opacity',0)<.999999:continue
        if any(opacity<.999999 for opacity in groups.values()):continue
        for item in drawing['items']:
            if item[0]!='re':continue
            rect=fitz.Rect(item[1]) & native_page.rect
            for bound in clips.values():rect &= bound
            if not rect.is_empty:masks.append((drawing['seqno'],rect))
    glyphs=defaultdict(deque)
    for span in native_page.get_texttrace():
        for code,gid,origin,bbox in span['chars']:
            rect=fitz.Rect(bbox)
            covered=any(seq>span['seqno'] and mask.contains(rect) for seq,mask in masks)
            glyphs[(round(origin[0],2),round(origin[1],2))].append(covered)
    class CoveredTextFilter(mupdf.PdfSanitizeFilterOptions2):
        def __init__(self):super().__init__();self.use_virtual_text_filter()
        def text_filter(self,ctx,ucs,length,trm,ctm,bbox,render,ca,CA):
            a=fitz.Matrix(trm.a,trm.b,trm.c,trm.d,trm.e,trm.f)
            b=fitz.Matrix(ctm.a,ctm.b,ctm.c,ctm.d,ctm.e,ctm.f)
            origin=a*b;key=(round(origin.e,2),round(height-origin.f,2))
            queue=glyphs.get(key)
            if not queue:return 0
            remove=queue.popleft()
            if remove:counters['occluded_old_glyphs']+=1
            return int(remove)
    if any(any(queue) for queue in glyphs.values()):
        covered=CoveredTextFilter()
        options=fitz._make_PdfFilterOptions(recurse=0,instance_forms=1,sanitize=1,sopts=covered)
        mupdf.pdf_filter_page_contents(pdfdoc,fitz._as_pdf_page(native_page),options)

    figure_payload=document.tobytes(garbage=2,clean=False,deflate=True,
        deflate_images=False,deflate_fonts=False,use_objstms=True,
        compression_effort=100,no_new_id=True,preserve_metadata=True)
    document.close()
    figure_path.write_bytes(figure_payload)
    print("Produced",figure_path,"bytes",figure_path.stat().st_size,
          "native vector/text retained",flush=True)


def _retained_form(document, marker):
    """Find the existing native form that owns a scoped layout operation."""
    seen = set()
    for xref, name, parent, box in document[0].get_xobjects():
        if xref in seen:
            continue
        seen.add(xref)
        stream = document.xref_stream(xref)
        if stream and marker in stream:
            return xref, stream
    raise RuntimeError("The expected retained native PDF form was not found.")


def _translate_native_x(fragment, delta):
    """Translate absolute native coordinates without changing fonts or scales."""
    import re
    number = rb"-?(?:\d+(?:\.\d*)?|\.\d+)"
    line = re.compile(rb"(?m)^(?:" + number + rb"[ \t]+)+(?:m|l|c|re|cm|Tm)[ \t]*$")

    def shifted(match):
        parts = match.group(0).split()
        operator = parts[-1]
        indices = {b"m": (0,), b"l": (0,), b"c": (0, 2, 4),
                   b"re": (0,), b"cm": (4,), b"Tm": (4,)}[operator]
        for index in indices:
            parts[index] = f"{float(parts[index]) + delta:.7f}".rstrip("0").rstrip(".").encode()
        return b" ".join(parts)

    return line.sub(shifted, fragment)


def _move_s9_row_labels(document):
    import re
    xref, stream = _retained_form(document, b"(ARACHNOID) Tj")
    moved = 0
    block = re.compile(rb"\bBT\b.*?\bET\b", re.S)
    matrix = re.compile(rb"(?m)^5 0 0 -5 ([\d.]+) ([\d.]+) Tm$")

    def reposition(match):
        nonlocal moved
        text = match.group(0)
        location = matrix.search(text)
        if location is None:
            return text
        x, y = map(float, location.groups())
        if 25 <= y <= 170 and (60 <= x <= 120 or 315 <= x <= 362):
            moved += 1
            new = f"5 0 0 -5 {x - 50:.7f} {y:.7f} Tm".encode()
            return text[:location.start()] + new + text[location.end():]
        return text

    changed = block.sub(reposition, stream)
    if moved != 30:
        raise RuntimeError(f"Expected 22 subclass and 8 domain labels; found {moved}.")
    document.update_stream(xref, changed)
    return "30 native row labels moved left; dendrograms and font operators unchanged"


def _move_s16_pair_panel(document):
    import re
    xref, stream = _retained_form(document, b"301.91166 40.248 m")
    start = stream.index(b"301.91166 40.248 m")
    boundary = b"\nQ\nQ\nq\n1 1 1 rg\n1 1 1 RG\n301.91166 188.5397"
    end = stream.index(boundary, start)
    delta = 48.0
    changed = stream[:start] + _translate_native_x(stream[start:end], delta) + stream[end:]
    # Move the existing e title and its existing replacement masks with the plot.
    # These are scoped coordinates, not a new white overlay or retyped label.
    for x, y in [("301.91166", "188.5397"), ("301.91166", "190.4367"),
                 ("360.89634", "20.75238"), ("360.19633", "19.91235"),
                 ("350.6814", "22.10736")]:
        old = f"{x} {y}".encode()
        if changed.count(old) != 1:
            raise RuntimeError("The retained S16e title/mask coordinates changed unexpectedly.")
        new = f"{float(x) + delta:.7f}".rstrip("0").rstrip(".").encode() + b" " + y.encode()
        changed = changed.replace(old, new, 1)
    document.update_stream(xref, changed)
    # The parent frame originally ended at the old figure's right edge. Expand
    # its clip into existing page whitespace; do not scale the figure or fonts.
    document.xref_set_key(xref, "BBox", "[0 0 566.4 619.2]")
    for parent, name, owner, box in document[0].get_xobjects():
        if parent == xref:
            continue
        resource_kind, resource = document.xref_get_key(parent, "Resources")
        if resource_kind == "xref":
            resource = document.xref_object(int(resource.split()[0]))
        if re.search(rf"/fullpage\s+{xref}\s+0\s+R", resource):
            document.xref_set_key(parent, "BBox", "[0 0 566.4 619.2]")
    return "S16e translated right without resizing; native d heatmap retained"


def repair_retained_supplement_layouts(names):
    """Produce only the owner's requested S9/S16 retained-layout repairs."""
    repairs = {"FigS9": _move_s9_row_labels, "FigS16": _move_s16_pair_panel}
    if not names or any(name not in repairs for name in names):
        raise RuntimeError("Select only FigS9 and/or FigS16 for this retained repair.")
    for name in names:
        pdf = ROOT / "source_figure_pdfs/supplementary_figures_pdf" / f"{name}.pdf"
        png = ROOT / "figures_png/supplementary_figures" / f"{name}.png"
        document = fitz.open(pdf)
        description = repairs[name](document)
        payload = document.tobytes(garbage=0, clean=False, deflate=True, no_new_id=True)
        document.close()
        pdf.write_bytes(payload)
        print("Produced", pdf, description, flush=True)
        with fitz.open(pdf) as rendered:
            rendered[0].get_pixmap(dpi=300, alpha=False).save(png)
        print("Produced", png, flush=True)
        number = int(name[4:]) + 1
        attachment = (ROOT / "supplementary_data/gigascience_supplementary_material"
                      / f"Additional_file_{number}_supplementary material_{name}.pdf")
        attachment.write_bytes(payload)
        print("Produced", attachment, flush=True)


def main():
    import sys
    if '--retained-supplement-layouts' in sys.argv:
        position = sys.argv.index('--retained-supplement-layouts')
        repair_retained_supplement_layouts(sys.argv[position + 1:])
        return
    if '--fig2-clean-native-only' in sys.argv:
        clean_fig2_native_content()
        return
    if '--fig5-s14-caption-only' in sys.argv:
        repair_s14_caption_only()
        return
    if '--fig5-supplements-only' in sys.argv:
        restore_fig5_supplements()
        return
    if "--gigascience-labels" in sys.argv:
        gigascience_labels()
        return
    if '--fig2-support-panel-only' in sys.argv:
        refresh_fig2_support_panel()
        return
    if '--fig2-label-mask-only' in sys.argv:
        repair_fig2_label_mask()
        return
    if '--fig2-subclass-typography-only' in sys.argv:
        repair_fig2_subclass_typography()
        return
    if '--fig2-subclass-reorder' in sys.argv:
        build_fig2_subclass_reorder()
        return
    if '--fig4-only' in sys.argv:
        build_fig4_native()
        return
    if '--fig3-only' in sys.argv:
        build_fig3()
        return
    if '--fig5-only' in sys.argv:
        build_fig5()
        return
    if '--fig6-only' in sys.argv:
        build_fig6()
        return
    if '--fig1-local-fg' in sys.argv:
        build_fig1_local_fg()
        return
    if '--fig1-original-only' in sys.argv:
        build_fig1_original()
        return
    if '--fig1-only' in sys.argv:
        build_fig1_with_overall_box()
        return
    if '--final-text-only' in sys.argv:
        build_fig3();build_s17();combine_supplements();return
    if "--combine-only" in sys.argv:
        combine_supplements();return
    return
if __name__=='__main__':main()
