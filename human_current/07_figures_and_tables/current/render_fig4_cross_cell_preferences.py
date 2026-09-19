#!/usr/bin/env python3
"""Directed 22-subclass program counts from retained Figure 4 association calls.

Program origin is the retained dominant_subclass annotation. The target is the
original spatial subclass. No new P values, FDR, pooled weights, or specificity
scores are estimated. Diagonal counts remain visible without colored fill.
"""
from pathlib import Path
import argparse
import csv

ROOT = Path((__import__("os").environ["NMF_WORK_ROOT"] + ""))
PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'analysis/fig4_cross_cell_preferences'
ANNOTATIONS = ROOT / 'tables/TableS3_program_annotation.tsv'
RELATIONS = ROOT / 'tables/TableS6_between_chip_colocalization.tsv'
COUNT_PATH = OUT / 'cell_type_pair_program_counts.tsv'

# The row organization of the retained original Figure 4a; also used for columns.
SUBCLASSES = [
    'OLIGO', 'SST', 'CHANDELIER', 'PVALB', 'VIP', 'NDNF', 'LAMP5', 'PAX6',
    'L6 CT', 'L6B', 'ET', 'NP', 'L6 IT', 'L6 CAR3', 'L2-L3 IT LINC00507',
    'L4-L5 IT RORB', 'L3-L4 IT RORB', 'ENDO', 'VLMC', 'MICRO', 'AST', 'OPC',
]
DISPLAY = {'OLIGO': 'Oligo', 'ENDO': 'Endo', 'MICRO': 'Micro', 'AST': 'Ast'}
CLASS_ORDER = ['Exc', 'Inh', 'Ast', 'Oligo', 'OPC', 'Micro', 'Endo', 'VLMC']
CLASS_MEMBERS = {
    'Exc': ['ET', 'NP', 'L6 CT', 'L6B', 'L6 IT', 'L6 CAR3',
            'L2-L3 IT LINC00507', 'L3-L4 IT RORB', 'L4-L5 IT RORB'],
    'Inh': ['CHANDELIER', 'PVALB', 'SST', 'VIP', 'NDNF', 'LAMP5', 'PAX6'],
    'Ast': ['AST'], 'Oligo': ['OLIGO'], 'OPC': ['OPC'],
    'Micro': ['MICRO'], 'Endo': ['ENDO'], 'VLMC': ['VLMC'],
}
CELL_CLASS = {cell: broad for broad, cells in CLASS_MEMBERS.items() for cell in cells}
CLASS_COLORS = {
    'Exc': '#4C78A8', 'Inh': '#8F63B8', 'Ast': '#2E9D58',
    'Oligo': '#C65353', 'OPC': '#B7791F', 'Micro': '#E88945',
    'Endo': '#795548', 'VLMC': '#8C8C8C',
}


def program_order(values):
    return sorted(values, key=lambda value: int(value[1:]))


def build_counts():
    with ANNOTATIONS.open() as handle:
        annotation = {row['new_P']: row for row in csv.DictReader(handle, delimiter='\t')}
    positive = {(source, target): set() for source in SUBCLASSES for target in SUBCLASSES}
    negative = {(source, target): set() for source in SUBCLASSES for target in SUBCLASSES}
    with RELATIONS.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['mode'] != 'cellprog' or row['is_headline'].lower() != 'true':
                continue
            pid = row['B_new_P']
            source = annotation[pid]['dominant_subclass']
            target = row['A_label']
            effect = float(row['median_log2g'])
            if effect > 0:
                positive[(source, target)].add(pid)
            elif effect < 0:
                negative[(source, target)].add(pid)
    records = []
    for source in SUBCLASSES:
        source_programs = {pid for pid, row in annotation.items() if row['dominant_subclass'] == source}
        for target in SUBCLASSES:
            pos, neg = positive[(source, target)], negative[(source, target)]
            together, mixed = pos | neg, pos & neg
            records.append({
                'reference_preference_subclass': source,
                'spatial_target_subclass': target,
                'n_source_programs': len(source_programs),
                'n_programs_any': len(together),
                'n_programs_positive': len(pos),
                'n_programs_negative': len(neg),
                'n_programs_both_signs': len(mixed),
                'program_ids_any': ';'.join(program_order(together)),
                'program_ids_positive': ';'.join(program_order(pos)),
                'program_ids_negative': ';'.join(program_order(neg)),
                'program_ids_both_signs': ';'.join(program_order(mixed)),
                'positive_functional_programs': '; '.join(
                    f"{annotation[pid]['functional_name']} ({pid})" for pid in program_order(pos)),
                'negative_functional_programs': '; '.join(
                    f"{annotation[pid]['functional_name']} ({pid})" for pid in program_order(neg)),
            })
    OUT.mkdir(parents=True, exist_ok=True)
    with COUNT_PATH.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    print('RETAINED', COUNT_PATH, flush=True)
    print('COUNT_GRAIN 22 reference-preference subclasses x 22 original spatial targets; 484 directed cells', flush=True)
    return records


def draw_overview(records):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    from matplotlib.collections import PolyCollection, LineCollection
    from matplotlib.patches import Rectangle, Patch, Polygon
    from matplotlib.transforms import blended_transform_factory

    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 8,
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.facecolor': 'white', 'axes.linewidth': 0.5,
    })
    lookup = {(row['reference_preference_subclass'], row['spatial_target_subclass']): row
              for row in records}
    n = len(SUBCLASSES)
    source_n = {source: int(lookup[(source, source)]['n_source_programs']) for source in SUBCLASSES}
    positive = np.asarray([[int(lookup[(source, target)]['n_programs_positive'])
                            for target in SUBCLASSES] for source in SUBCLASSES])
    negative = np.asarray([[int(lookup[(source, target)]['n_programs_negative'])
                            for target in SUBCLASSES] for source in SUBCLASSES])
    warm = LinearSegmentedColormap.from_list('positive_count',
        ['#fff7ec', '#fee8c8', '#fdbb84', '#fc8d59', '#e34a33', '#b30000'])
    cool = LinearSegmentedColormap.from_list('negative_count',
        ['#effaf8', '#ccece6', '#9adbd3', '#64c4c2', '#2a9eaa', '#086b8c'])
    norm = Normalize(0, 5)
    fig = plt.figure(figsize=(240 / 25.4, 205 / 25.4), facecolor='white')
    ax = fig.add_axes([0.208, 0.213, 0.638, 0.747])
    ax.set_xlim(0, n)
    ax.set_ylim(n, 0)
    ax.set_aspect('equal')
    polygons, facecolors, diagonals = [], [], []
    for row in range(n):
        for column in range(n):
            top_left = [(column, row), (column + 1, row), (column, row + 1)]
            bottom_right = [(column + 1, row), (column + 1, row + 1), (column, row + 1)]
            polygons.extend([top_left, bottom_right])
            diagonal = row == column
            pos_color = (1, 1, 1, 1) if diagonal else warm(norm(positive[row, column]))
            neg_color = (1, 1, 1, 1) if diagonal else cool(norm(negative[row, column]))
            facecolors.extend([pos_color, neg_color])
            diagonals.append([(column, row + 1), (column + 1, row)])
            for value, x_offset, y_offset, fill in [
                (positive[row, column], 0.29, 0.29, pos_color),
                (negative[row, column], 0.71, 0.71, neg_color),
            ]:
                luminance = 0.299 * fill[0] + 0.587 * fill[1] + 0.114 * fill[2]
                color = '#222222' if diagonal else '#9b9b9b' if value == 0 else 'white' if luminance < 0.53 else '#222222'
                ax.text(column + x_offset, row + y_offset, str(value),
                        ha='center', va='center', fontsize=7.7,
                        color=color, fontweight='normal', zorder=3)
    ax.add_collection(PolyCollection(polygons, facecolors=facecolors, edgecolors='none', zorder=0))
    boundaries = [[(line, 0), (line, n)] for line in range(n + 1)]
    boundaries += [[(0, line), (n, line)] for line in range(n + 1)]
    ax.add_collection(LineCollection(boundaries, colors='#e2e2e2', linewidths=0.20, zorder=1))
    ax.add_collection(LineCollection(diagonals, colors='#bdbdbd', linewidths=0.28, zorder=2))
    ax.set_xticks(np.arange(n) + 0.5, [DISPLAY.get(value, value) for value in SUBCLASSES],
                  rotation=65, ha='right', rotation_mode='anchor', fontsize=8.0)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels([])
    ax.set_xlabel('Spatial target subclass', fontsize=10, labelpad=9)
    ax.tick_params(length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # A small aligned class swatch precedes each complete shared row label.
    row_transform = blended_transform_factory(fig.transFigure, ax.transData)
    for row, value in enumerate(SUBCLASSES):
        ax.add_patch(Rectangle((0.046, row + 0.38), 0.00695, 0.24,
                               transform=row_transform, clip_on=False,
                               facecolor=CLASS_COLORS[CELL_CLASS[value]], edgecolor='none'))
        ax.text(0.059, row + 0.5, f'{DISPLAY.get(value, value)} ({source_n[value]})',
                transform=row_transform, clip_on=False, ha='left', va='center',
                fontsize=8.0, color='#222222')
    position = ax.get_position()
    fig.text(0.018, position.y0 + position.height / 2,
             'Program preference\n(available programs)', rotation=90,
             fontsize=9.5, ha='center', va='center')
    fig.legend(handles=[Patch(facecolor=CLASS_COLORS[broad], edgecolor='none', label=broad)
                        for broad in CLASS_ORDER],
               loc='upper center', bbox_to_anchor=(0.54, 0.999), ncol=8, frameon=False,
               fontsize=8.1, handlelength=0.7, handleheight=0.7, handletextpad=0.4,
               columnspacing=0.85, borderaxespad=0)

    # Direction key follows exactly the same lower-left to upper-right split.
    key = fig.add_axes([0.873, 0.792, 0.095, 0.111])
    key.set_xlim(0, 1)
    key.set_ylim(1, 0)
    key.set_aspect('equal')
    key.add_patch(Polygon([(0, 0), (1, 0), (0, 1)], facecolor=warm(norm(3)), edgecolor='none'))
    key.add_patch(Polygon([(1, 0), (1, 1), (0, 1)], facecolor=cool(norm(3)), edgecolor='none'))
    key.plot([0, 1], [1, 0], color='#a4a4a4', linewidth=0.6)
    key.text(0.29, 0.27, 'Positive', fontsize=7.1, ha='center', va='center')
    key.text(0.71, 0.73, 'Negative', fontsize=7.1, ha='center', va='center')
    key.set_title('Direction', fontsize=9.0, pad=6)
    key.set_axis_off()
    for bottom, cmap, label in [(0.660, warm, 'Positive programs'), (0.540, cool, 'Negative programs')]:
        cax = fig.add_axes([0.880, bottom, 0.080, 0.015])
        colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                                cax=cax, orientation='horizontal', ticks=range(6))
        colorbar.ax.tick_params(labelsize=7.2, length=2, width=0.45, pad=2)
        colorbar.outline.set_visible(False)
        cax.set_xlabel(label, fontsize=8.1, labelpad=5)
        cax.xaxis.set_label_position('top')

    pdf_path = OUT / 'Fig4_cross_cell_preferences.pdf'
    png_path = OUT / 'Fig4_cross_cell_preferences.png'
    fig.savefig(pdf_path, format='pdf')
    print('RETAINED', pdf_path, flush=True)
    fig.savefig(png_path, format='png', dpi=300)
    plt.close(fig)
    print('RETAINED', png_path, flush=True)
    print('DISPLAY one directed 22x22 split-cell matrix; upper-left positive and lower-right negative; separate warm/cool count scales 0-5; identity diagonal white with both counts; row class swatches and complete eight-class key', flush=True)



def add_svg_spatial(fig, axes_mm, text_mm, path, x, y, width, titles):
    """Reuse the source SVG's embedded image layers, alpha and colorbar pixels."""
    import base64, io
    import xml.etree.ElementTree as ET
    import numpy as np
    from PIL import Image
    from matplotlib.lines import Line2D
    root=ET.parse(path).getroot()
    is_m=Path(path).name=='figA_m.svg'
    left,right=(48.82,332.01) if is_m else (51.33,329.50)
    scale=width/(right-left)
    W,H=fig.get_size_inches()*25.4
    for node in root.iter():
        kind=node.tag.rsplit('}',1)[-1]
        if kind=='image':
            ix=float(node.attrib['x']);iy=float(node.attrib['y'])
            if ix<left or iy<15:continue
            iw=float(node.attrib['width']);ih=float(node.attrib['height'])
            href=next(v for k,v in node.attrib.items() if k.endswith('href'))
            arr=np.asarray(Image.open(io.BytesIO(base64.b64decode(href.split(',',1)[1]))))
            ax=axes_mm(x+(ix-left)*scale,y+(iy-15.85)*scale,iw*scale,ih*scale)
            ax.set_zorder(7)
            ax.imshow(arr,aspect='auto',interpolation='none');ax.set_axis_off()
        elif kind=='polyline':
            points=[tuple(map(float,p.split(','))) for p in node.attrib['points'].split()]
            if all(py>=100 for px,py in points):
                fig.lines.append(Line2D([(x+(px-left)*scale)/W for px,py in points],
                    [1-(y+(py-15.85)*scale)/H for px,py in points],
                    transform=fig.transFigure,color='white',linewidth=.38,zorder=8))
    for tx,title in zip([92.58,190.42,288.25],titles):
        text_mm(x+(tx-left)*scale,y-9,title,ha='center',fontsize=8.0,linespacing=1.05,zorder=9)
    for sx,label in [(98.33,'Low'),(168.96,'High'),(188.14,'Low'),(258.77,'High')]:
        text_mm(x+(sx-left)*scale,y+(109-15.85)*scale,label,ha='center',fontsize=7.0,zorder=9)
    for sx,label in [(133.64,'Relative cell weight'),(223.46,'Relative SCT score')]:
        text_mm(x+(sx-left)*scale,y+(118-15.85)*scale,label,ha='center',fontsize=7.1,zorder=9)


def draw_reorganized(records, base_pdf, compact_candidate=False):
    """Local replacements on the original Figure 4 evidence, without refitting."""
    import io, copy, math, subprocess
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    from matplotlib.collections import PolyCollection, LineCollection
    from matplotlib.patches import Rectangle, Patch, Polygon
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory
    from pypdf import PdfReader, PdfWriter, Transformation
    from pypdf.generic import ContentStream, DecodedStreamObject, NameObject, RectangleObject, DictionaryObject

    with ANNOTATIONS.open() as handle:
        annotation={r['new_P']:r for r in csv.DictReader(handle,delimiter='\t')}
    with RELATIONS.open() as handle:
        relations={(r['A_label'],r['B_new_P']):r for r in csv.DictReader(handle,delimiter='\t') if r['mode']=='cellprog'}
    W,H,MM=(170.0,234.0,72/25.4) if compact_candidate else (169.4215,250.37,72/25.4)
    original_reader=PdfReader(io.BytesIO(base_pdf)); original=original_reader.pages[0]
    reader_holders=[original_reader]
    def mul(a,b):
        return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
                a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3],
                a[4]*b[0]+a[5]*b[2]+b[4],a[4]*b[1]+a[5]*b[3]+b[5])
    text_ops={b'BT',b'ET',b'Tf',b'Tm',b'Td',b'TD',b'T*',b'Tj',b'TJ',b'Tc',b'Tw',b'Tz',b'TL',b'Ts',b'Tr',b"'",b'"'}
    def variant(no_text=False,no_rotated=False,no_orange=False,no_example_labels=False):
        reader=PdfReader(io.BytesIO(base_pdf)); page=reader.pages[0]; reader_holders.append(reader)
        def orange(color):
            return len(color)==3 and all(abs(x-y)<.004 for x,y in zip(color,(.85098,.37255,.00784)))
        def rewrite(stream,initial=(1,0,0,1,0,0)):
            cs=ContentStream(stream,reader); result=[]; stack=[]
            ctm=initial; stroke=(0,0,0); fill=(0,0,0); buffer=None; rotated=False
            path_points=[]; path_closed=False; text_position=None
            for operands,operator in cs.operations:
                if operator==b'q': stack.append((ctm,stroke,fill))
                elif operator==b'Q' and stack: ctm,stroke,fill=stack.pop()
                elif operator==b'cm': ctm=mul(tuple(map(float,operands)),ctm)
                elif operator==b'RG': stroke=tuple(map(float,operands))
                elif operator==b'rg': fill=tuple(map(float,operands))
                elif operator==b'G': stroke=(float(operands[0]),)*3
                elif operator==b'g': fill=(float(operands[0]),)*3
                if operator in (b'm',b'l',b'c',b'v',b'y'):
                    if operator==b'm': path_points=[]; path_closed=False
                    vals=list(map(float,operands))
                    path_points.extend((vals[i]*ctm[0]+vals[i+1]*ctm[2]+ctm[4],
                                        721-(vals[i]*ctm[1]+vals[i+1]*ctm[3]+ctm[5])) for i in range(0,len(vals),2))
                elif operator==b're':
                    path_closed=True
                elif operator==b'h': path_closed=True
                if operator==b'BT': buffer=[(operands,operator)]; rotated=False; text_position=None; continue
                if buffer is not None:
                    buffer.append((operands,operator))
                    if operator==b'Tm':
                        tm=tuple(map(float,operands)); dx=tm[0]*ctm[0]+tm[1]*ctm[2]; dy=tm[0]*ctm[1]+tm[1]*ctm[3]
                        rotated=abs(dy)>3*abs(dx)
                        text_position=(tm[4]*ctm[0]+tm[5]*ctm[2]+ctm[4],721-(tm[4]*ctm[1]+tm[5]*ctm[3]+ctm[5]))
                    example_label=no_example_labels and text_position is not None and 340<text_position[0]<539 and 299<text_position[1]<383
                    if (example_label or (no_rotated and rotated)) and operator in (b'Tj',b'TJ',b"'",b'"'):
                        buffer.pop()
                    if False and no_rotated and rotated and operator in (b'Tj',b'TJ',b"'",b'"'):
                        buffer.pop()
                    if operator==b'ET':
                        if no_text:
                            result.extend(item for item in buffer if item[1] not in text_ops)
                        else: result.extend(buffer)
                        buffer=None
                    continue
                remove = no_orange and (
                    (operator in (b'S',b's') and orange(stroke)) or
                    (operator in (b'f',b'F',b'f*') and orange(fill)) or
                    (operator in (b'B',b'B*',b'b',b'b*') and (orange(fill) or orange(stroke))))
                if no_example_labels and operator in (b'S',b's') and not path_closed and len(path_points)>=2 and max(stroke)<.30:
                    x0=min(v[0] for v in path_points); x1=max(v[0] for v in path_points)
                    y0=min(v[1] for v in path_points); y1=max(v[1] for v in path_points)
                    remove=remove or (340<x0<x1<539 and 299<y0<y1<390 and x1-x0>1 and y1-y0>1)
                result.append(([],b'n') if remove else (operands,operator))
                if operator in (b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*',b'n'):
                    path_points=[]; path_closed=False
            cs.operations=result
            return cs
        seen=set()
        def forms(resources):
            objects=resources.get('/XObject',{})
            if hasattr(objects,'get_object'):objects=objects.get_object()
            for name,ref in list(objects.items()):
                obj=ref.get_object()
                if obj.get('/Subtype')!='/Form' or id(obj) in seen:continue
                seen.add(id(obj)); child=obj.get('/Resources',resources)
                if hasattr(child,'get_object'):child=child.get_object()
                forms(child); replacement=DecodedStreamObject()
                for key,value in obj.items():
                    if key not in ('/Length','/Filter','/DecodeParms'):replacement[key]=value
                replacement.set_data(rewrite(obj,tuple(map(float,obj.get('/Matrix',[1,0,0,1,0,0])))).get_data())
                objects[name]=replacement
        forms(page['/Resources']); page[NameObject('/Contents')]=rewrite(page.get_contents())
        return page
    matrix_page=variant(no_text=True)

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.55})
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='none')
    def axes_mm(x,y,w,h):
        return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],facecolor='none')
    def text_mm(x,y,text,**kwargs):
        opts=dict(ha='left',va='top',fontsize=7.5);opts.update(kwargs)
        return fig.text(x/W,1-y/H,text,**opts)
    def white_mm(x,y,w,h):
        fig.patches.append(Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,
                                     facecolor='white',edgecolor='none',zorder=6))
    panel_positions=([('a',.5,.1),('b',.8,113.7),('c',105.8,113.7),
                      ('d',.8,150.2),('e',85.8,150.2),('f',.8,186.2),
                      ('g',85.8,186.2),('h',.8,221.1)]
                     if compact_candidate else
                     [('a',.5,.3),('b',.9,132.1),('c',109,132.1),('d',.9,168.2),('e',85.6,168.2),('f',.9,202.2),('g',85.6,202.2),('h',.9,236.0)])
    for letter,x,y in panel_positions:
        text_mm(x,y,letter,fontsize=8.5,fontweight='bold',zorder=30)
    # Owner-requested light grouping boundaries: one per evidence group,
    # never a second frame around individual scientific plots.
    group_boxes=([(.5,150.0,84,34.0),(85.5,150.0,84,34.0),
                  (.5,184.5,84,34.0),(85.5,184.5,84,34.0),
                  (.5,113.5,104,36.0),(105.5,113.5,64,36.0)]
                 if compact_candidate else
                 [(.5,168,83.7,33.3),(85.2,168,83.7,33.3),
                  (.5,202,83.7,33.3),(85.2,202,83.7,33.3),
                  (.5,132,107.0,35.0),(108.5,132,60.4,35.0)])
    for xx,yy,ww,hh in group_boxes:
        fig.patches.append(Rectangle((xx/W,1-(yy+hh)/H),ww/W,hh/H,
            transform=fig.transFigure,facecolor='none' if compact_candidate else '#F7F8FA',edgecolor='none' if compact_candidate else '#D4D8E0',
            linewidth=.45,zorder=-10))
    crops=[]
    def native(page,box,x,y,width,clockwise=False,display_height=None):
        if clockwise:height=width*(box[2]-box[0])/(box[3]-box[1])
        else:height=width*(box[3]-box[1])/(box[2]-box[0])
        crops.append((page,box,x,y,width,clockwise,display_height))
        return height


    # a is the original matrix and both native trees rotated clockwise 90°.
    # It is not a transpose/re-cluster computation: original vector marks and
    # their geometry are preserved. The two old thin annotation strips are
    # outside these clips, rather than drawn again beside the named bands.
    original_ox,original_sx=(47.0,.730) if compact_candidate else (60.0,.555)
    OX,SY1=(46.4,215.0) if compact_candidate else (54.5,215.0)
    # Reclaim only the reference-cell band/gap; keep the matrix's right edge
    # and the program tree in place while enlarging the actual matrix.
    SX=original_sx+(original_ox-OX)/(SY1-76.75)
    col_ids=[47,34,12,23,33,32,41,21,8,27,38,30,25,35,20,18,17,29,26,43,37,19,1,11,10,15,13,3,6,
             28,52,22,24,46,9,7,16,2,5,42,44,50,51,53,36,45,49,54,31,4,39,40,14,48]
    col_starts=[51.751,59.943,65.791,71.626,77.471,83.307,89.143,97.333,103.182,111.373,117.208,
                125.404,131.240,137.083,142.920,148.756,154.600,160.436,166.279,172.115,177.958,
                183.795,189.631,195.479,201.314,207.158,212.994,218.830,224.674,232.869,238.705,
                244.549,252.733,260.928,266.771,274.963,280.803,288.994,297.189,305.381,311.217,
                317.053,322.896,328.732,336.928,342.771,348.607,354.442,362.639,370.830,379.024,
                387.221,393.057,398.899]

    # Keep every original row and its gap/order; compress only decorative
    # inter-group spacing. Names are set in final physical points, not scaled.
    row_scale=1.75/5.836 if compact_candidate else .345
    gap_scale=.025 if compact_candidate else .065
    row_tops=[4.2 if compact_candidate else 6.5]
    for left,right in zip(col_starts[:-1],col_starts[1:]):
        extra=max(0.0,right-left-5.836)
        row_tops.append(row_tops[-1]+5.836*row_scale+extra*gap_scale)
    row_anchor=np.asarray(col_starts+[col_starts[-1]+5.836])
    row_display=np.asarray(row_tops+[row_tops[-1]+5.836*row_scale])
    def program_y(x):
        if x<row_anchor[0]:return row_display[0]+(x-row_anchor[0])*row_scale
        if x>row_anchor[-1]:return row_display[-1]+(x-row_anchor[-1])*row_scale
        if compact_candidate:
            i=min(len(col_starts)-1,int(np.searchsorted(col_starts,x,side='right'))-1)
            offset=x-col_starts[i]
            return row_tops[i]+min(offset,5.836)*row_scale+max(0,offset-5.836)*gap_scale
        return float(np.interp(x,row_anchor,row_display))

    # Extract the original heatmap/tree paths in their exact page coordinates.
    # Cell colors, borders, all mark centers and both tree topologies are reused;
    # only their plotting coordinates change. Original selection centres
    # identify supported cells; only cross-subclass cells receive open frames.
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    original_paths=[]
    def source_paths(stream,resources,initial=(1,0,0,1,0,0),initial_style=None):
        state={'ctm':initial,'fill':(0,0,0),'stroke':(0,0,0),'width':1.0,'alpha':1.0}
        if initial_style:state.update(initial_style);state['ctm']=initial
        stack=[];vertices=[];codes=[];has_curve=False;substart=None;last=None
        def point(x,y):
            m=state['ctm'];return (x*m[0]+y*m[2]+m[4],721-(x*m[1]+y*m[3]+m[5]))
        for operands,op in ContentStream(stream,original_reader).operations:
            if op==b'q':stack.append(state.copy())
            elif op==b'Q' and stack:state=stack.pop()
            elif op==b'cm':state['ctm']=mul(tuple(map(float,operands)),state['ctm'])
            elif op==b'rg':state['fill']=tuple(map(float,operands))
            elif op==b'RG':state['stroke']=tuple(map(float,operands))
            elif op==b'g':state['fill']=(float(operands[0]),)*3
            elif op==b'G':state['stroke']=(float(operands[0]),)*3
            elif op==b'w':state['width']=float(operands[0])
            elif op==b'gs':
                states=resources.get('/ExtGState',{})
                if hasattr(states,'get_object'):states=states.get_object()
                gs=states.get(operands[0])
                if gs:state['alpha']=float(gs.get_object().get('/ca',state['alpha']))
            elif op==b'Do':
                objects=resources.get('/XObject',{})
                if hasattr(objects,'get_object'):objects=objects.get_object()
                obj=objects[operands[0]].get_object()
                if obj.get('/Subtype')=='/Form':
                    child=obj.get('/Resources',resources)
                    if hasattr(child,'get_object'):child=child.get_object()
                    source_paths(obj,child,mul(tuple(map(float,obj.get('/Matrix',[1,0,0,1,0,0]))),state['ctm']),state)
            elif op==b'm':
                last=point(*map(float,operands));substart=last
                vertices.append(last);codes.append(MplPath.MOVETO)
            elif op==b'l':
                last=point(*map(float,operands));vertices.append(last);codes.append(MplPath.LINETO)
            elif op==b'c':
                vals=list(map(float,operands));pts=[point(vals[i],vals[i+1]) for i in (0,2,4)]
                vertices.extend(pts);codes.extend([MplPath.CURVE4]*3);last=pts[-1];has_curve=True
            elif op==b'v':
                vals=list(map(float,operands));pts=[last,point(*vals[:2]),point(*vals[2:])]
                vertices.extend(pts);codes.extend([MplPath.CURVE4]*3);last=pts[-1];has_curve=True
            elif op==b'y':
                vals=list(map(float,operands));pts=[point(*vals[:2]),point(*vals[2:]),point(*vals[2:])]
                vertices.extend(pts);codes.extend([MplPath.CURVE4]*3);last=pts[-1];has_curve=True
            elif op==b're':
                x,y,w,h=map(float,operands)
                pts=[point(x,y),point(x+w,y),point(x+w,y+h),point(x,y+h),point(x,y)]
                vertices.extend(pts);codes.extend([MplPath.MOVETO]+[MplPath.LINETO]*3+[MplPath.CLOSEPOLY]);last=pts[0];substart=last
            elif op==b'h' and substart is not None:
                vertices.append(substart);codes.append(MplPath.CLOSEPOLY);last=substart
            elif op in (b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*',b'n'):
                if vertices and op!=b'n':
                    vv=np.asarray(vertices);x0,y0=vv.min(axis=0);x1,y1=vv.max(axis=0)
                    kind=None
                    if 49.95<=x0<=x1<=405.55 and 76.75<=y0<=y1<=215.05:kind='matrix'
                    elif 49.95<=x0<=x1<=405.55 and 42.45<=y0<=y1<=69.95:kind='program_tree'
                    elif 17.95<=x0<=x1<=43.65 and 76.75<=y0<=y1<=215.05:kind='target_tree'
                    if kind:
                        fill=state['fill'] if op in (b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*') else None
                        stroke=state['stroke'] if op in (b'S',b's',b'B',b'B*',b'b',b'b*') else None
                        original_paths.append((kind,vv,list(codes),fill,stroke,state['width'],state['alpha'],has_curve))
                vertices=[];codes=[];has_curve=False;substart=None;last=None
    source_paths(original.get_contents(),original['/Resources'])
    row_ranges=[(77.823,83.280),(85.628,91.081),(91.081,96.538),(96.546,101.999),
                (101.999,107.452),(107.452,112.909),(112.909,118.362),(118.362,123.815),
                (126.175,131.628),(131.628,137.081),(137.081,142.534),(142.534,147.991),
                (150.347,155.804),(155.804,161.257),(161.257,166.710),(166.710,172.163),
                (172.163,177.620),(179.976,185.429),(185.429,190.886),(193.233,198.686),
                (201.046,206.499),(208.858,214.311)]
    dot_count=0;cross_frame_count=0;same_subclass_count=0
    for kind,vv,codes,fill,stroke,lw,alpha,curved in original_paths:
        def display_point(v):
            if compact_candidate:
                y=program_y(v[0]) if kind!='target_tree' else .5+(v[0]-18)/(43.6-18)*2.7
                x=149.0+(69.9-v[1])/(69.9-42.5)*7.5 if kind=='program_tree' else OX+(SY1-v[1])*SX
                return (x,y)
            y=program_y(v[0]) if kind!='target_tree' else .7+(v[0]-18)/(43.6-18)*4.3
            x=original_ox+(SY1-v[1])*original_sx if kind=='program_tree' else OX+(SY1-v[1])*SX
            return (x,y)
        bounds=vv.max(axis=0)-vv.min(axis=0)
        is_dot=kind=='matrix' and curved and fill is not None and max(fill)<.12 and max(bounds)<3.0
        if is_dot:
            dot_count+=1
            source_x,source_y=(vv.max(axis=0)+vv.min(axis=0))/2
            program_matches=[(pid,sx) for pid,sx in zip(col_ids,col_starts)
                             if sx-.04<=source_x<=sx+5.836+.04]
            target_matches=[(target,y0,y1) for target,(y0,y1) in zip(SUBCLASSES,row_ranges)
                            if y0-.04<=source_y<=y1+.04]
            if len(program_matches)!=1 or len(target_matches)!=1:
                raise ValueError('Original selected cell cannot be mapped uniquely: '+
                                 str((float(source_x),float(source_y))))
            pid,sx=program_matches[0];target,y0,y1=target_matches[0]
            preference=annotation['P'+str(pid)]['dominant_subclass']
            if preference not in SUBCLASSES:
                raise ValueError('Unknown formal dominant_subclass for P'+str(pid)+': '+str(preference))
            if preference==target:
                same_subclass_count+=1
                continue
            # Compare the 22 formal subclass identities, not broad classes or
            # shortened display labels. The source selection mask is untouched.
            xx=OX+(SY1-y1)*SX;ww=(y1-y0)*SX
            yy=program_y(sx);hh=program_y(sx+5.836)-yy
            inset=.10
            fig.patches.append(Rectangle(((xx+inset)/W,1-(yy+hh-inset)/H),
                (ww-2*inset)/W,(hh-2*inset)/H,transform=fig.transFigure,
                facecolor='none',edgecolor='black',linewidth=.65,zorder=5))
            cross_frame_count+=1
        else:
            transformed=[(display_point(v)[0]/W,1-display_point(v)[1]/H) for v in vv]
            fig.patches.append(PathPatch(MplPath(transformed,codes),transform=fig.transFigure,
                facecolor=fill if fill is not None else 'none',
                edgecolor=stroke if stroke is not None else 'none',
                linewidth=max(.12,min(.60,lw*.78)),alpha=alpha,zorder=3,capstyle='butt',joinstyle='miter'))
    print('ORIGINAL_A_VECTOR_PATHS',len(original_paths),
          'source_selected_cells',dot_count,'cross_subclass_frames',cross_frame_count,
          'same_subclass_unmarked',same_subclass_count,flush=True)
    original_low_confidence={12,32,27,38,25,35,20,18,29,26,43,37,11,10,15,6,24,2,50,4}
    name_x,pref_x,pref_w,pref_center,group_x=(33.4,34,12.0,40.0,156.8) if compact_candidate else (40,40.5,13.5,47.25,157)
    text_mm(name_x,.2 if compact_candidate else 2.2,'Program',ha='right',fontsize=7 if compact_candidate else 6.2)
    text_mm(pref_center,.1 if compact_candidate else .7,'Ref. cell\npreference' if compact_candidate else 'Reference cell\npreference',ha='center',fontsize=5.05 if compact_candidate else 5.8,linespacing=1.0)
    text_mm(group_x,.1 if compact_candidate else .7,'Functional\ngroup',fontsize=5.05 if compact_candidate else 5.0,linespacing=1.0)
    # Current owner style: readable word-level abbreviations, never a fixed
    # character cutoff or whole-phrase ellipsis. Table S3 keeps the full meaning.
    import re
    word_abbreviations={
        'oligodendrocyte':'oligo.','development':'dev.','myelination':'myelinat.',
        'reactive':'react.','astrocyte':'astro.','vascular':'vasc.',
        'neurofilament':'neurofil.','cytoskeleton':'cytoskel.','neuronal':'neuron.',
        'synaptic':'syn.','postsynaptic':'postsyn.','vesicle':'ves.',
        'exocytosis':'exocyt.','secretion':'secret.','acetylcholine':'ACh',
        'receptor':'recept.','receptors':'recept.','voltage':'volt.',
        'channel':'chan.','channels':'chan.','adhesion':'adhes.','assembly':'asm.',
        'signaling':'sig.','glutamate':'Glu','ionotropic':'ionotrop.',
        'calcium':'Ca','homeostasis':'homeost.','guidance':'guid.','neural':'neur.',
        'splicing':'splic.','release':'rel.','outgrowth':'outgr.',
        'neuropeptide':'neuropept.','specificity':'spec.','organization':'org.',
        'interneuron':'IN','sprouting':'sprout.','angiogenesis':'angiogen.',
        'chemotaxis':'chemotax.','vessel':'vess.','morphogenesis':'morphog.',
        'microglial':'microgl.','immune':'immun.','activation':'activ.',
        'complement':'compl.','chemokine':'chemok.','activity':'act.',
        'dependent':'depend.','proteoglycan':'proteoglyc.','transport':'transp.'}
    word_pattern=re.compile(r'(?<![A-Za-z])('+
        '|'.join(sorted(word_abbreviations,key=len,reverse=True))+r')(?![A-Za-z.])',re.I)
    def abbreviated_program_name(pid,name):
        def shorten_word(match):
            original_word=match.group(0)
            short=word_abbreviations[original_word.lower()]
            return short[0].upper()+short[1:] if original_word[0].isupper() else short
        name=word_pattern.sub(shorten_word,name)
        if pid==8:
            name='Pan-neuron. neurofil. cytoskel.'
        return name
    def a_cell_display(cell):
        return {'L2-L3 IT LINC00507':'L2-L3 IT',
                'L3-L4 IT RORB':'L3-L4 IT',
                'L4-L5 IT RORB':'L4-L5 IT'}.get(cell,DISPLAY.get(cell,cell))
    shortened_rows=[];unchanged_rows=[]
    for pid,sx in zip(col_ids,col_starts):
        functional=annotation['P'+str(pid)]['functional_name']
        display_name=abbreviated_program_name(pid,functional)
        if display_name!=functional:shortened_rows.append((pid,functional,display_name))
        else:unchanged_rows.append(pid)
        label=display_name+' (P'+str(pid)+')'+('*' if pid in original_low_confidence else '')
        text_mm(name_x,program_y(sx+2.918),label,ha='right',va='center',fontsize=5.05 if compact_candidate else 6.0)
    print('A_WORD_ABBREVIATIONS','changed',len(shortened_rows),'unchanged',len(unchanged_rows),
          'rows',shortened_rows,flush=True)
    print('A_IT_DISPLAY','L2-L3 IT LINC00507 -> L2-L3 IT; L3-L4 IT RORB -> L3-L4 IT; L4-L5 IT RORB -> L4-L5 IT; internal keys unchanged',flush=True)
    runs=[]
    for pid,sx in zip(col_ids,col_starts):
        preference=annotation['P'+str(pid)]['dominant_subclass']
        if runs and runs[-1][0]==preference and abs(sx-runs[-1][2])<.03:
            runs[-1]=(preference,runs[-1][1],sx+5.836)
        else:runs.append((preference,sx,sx+5.836))
    for preference,x0,x1 in runs:
        yy=program_y(x0);hh=program_y(x1)-yy
        fig.patches.append(Rectangle((pref_x/W,1-(yy+hh)/H),pref_w/W,hh/H,transform=fig.transFigure,
                          facecolor=CLASS_COLORS[CELL_CLASS[preference]],edgecolor='white',linewidth=.18,zorder=7))
        foreground='#222222' if CELL_CLASS[preference] in ['Micro','VLMC'] else 'white'
        label=a_cell_display(preference)
        text_mm(pref_center,yy+hh/2,label,ha='center',va='center',fontsize=5.05 if compact_candidate else 5.1,
                color=foreground,zorder=8,linespacing=.95)
    functional_groups=[('Astro-react',1),('Oligo',6),('IT-neuropil',2),('Secretory',2),
        ('Neuron',18),('Neuropeptide',3),('Mixed',1),('ECM',2),('Inh-neuron',2),
        ('Cilia',1),('Lipid',1),('Vascular',5),('Micro',4),('Neuron-IEG',1),
        ('Mito',1),('OPC/ECM',1),('Astro',3)]
    cursor=0
    for label,n in functional_groups:
        center=(program_y(col_starts[cursor])+program_y(col_starts[cursor+n-1]+5.836))/2
        text_mm(group_x,center,label,va='center',fontsize=5.05 if compact_candidate else (4.6 if label=='Neuron-IEG' else 5.0));cursor+=n

    # Clockwise rotation makes target columns run from original bottom to top.
    # The complete subclass names are inside this single named color band.
    target_y,target_h=(100.8,11.7) if compact_candidate else (118.7,10.0)
    for target,(y0,y1) in zip(SUBCLASSES,row_ranges):
        xx=OX+(SY1-y1)*SX;ww=(y1-y0)*SX
        fig.patches.append(Rectangle((xx/W,1-(target_y+target_h)/H),ww/W,target_h/H,
                          transform=fig.transFigure,facecolor=CLASS_COLORS[CELL_CLASS[target]],
                          edgecolor='white',linewidth=.2,zorder=7))
        foreground='#222222' if CELL_CLASS[target] in ['Micro','VLMC'] else 'white'
        text_mm(xx+ww/2,target_y+target_h/2,a_cell_display(target),rotation=90,
                ha='center',va='center',fontsize=5.05 if compact_candidate else 4.6,color=foreground,zorder=8)
    if compact_candidate:
        text_mm(45.4,106.65,'Target subclass',ha='center',va='center',rotation=90,fontsize=5.05)
        key_x,key_y,key_w,key_h=149.0,102.8,2.4,8.8
    else:
        text_mm(OX+(SY1-(77.823+214.311)/2)*SX,129.3,'Spatial target subclass',ha='center',fontsize=5.8)
        key_x,key_y,key_w,key_h=141.0,119.8,2.3,10.8
    native(matrix_page,(457.234,243.031,466.648,294.800),key_x,key_y,key_w,display_height=key_h)
    for label,source_y in [('1',243.934),('0.5',256.045),('0',268.160),('-0.5',280.279),('-1',292.390)]:
        text_mm(key_x+key_w+.7,key_y+(source_y-243.031)/(294.800-243.031)*key_h,label,
                va='center',fontsize=5.05 if compact_candidate else 5.5)
    if compact_candidate:
        text_mm(156,100.8,'log₂ median g\n([0,25] µm)',fontsize=5.05,linespacing=1.0)
        text_mm(156,107.0,'red: co-occ.\nblue: avoid.',fontsize=5.05,linespacing=1.0)
    else:
        text_mm(149,119.5,r'$\log_2$ median $g$'+'\n'+r'([0,25] $\mu$m)',fontsize=5.8,linespacing=1.05)
        text_mm(149,126,'red: co-occ.\nblue: avoid.',fontsize=5.2,linespacing=1.15)

    # c: the original effect/consistency coordinates are redrawn below in a
    # taller plotting frame; type and point glyphs retain ordinary proportions.

    source_root=Path((__import__("os").environ["CORTEX_PROGRAM_ROOT"] + ""))
    spatial_root=source_root/'results/crossregion_v1'
    svg_root=source_root/'scripts/figmarkcorr_A/svg_panels'
    # Existing global, full-precision curve fields; no aggregation is performed.
    curves=np.load(spatial_root/'markcorr_v2/final/cellprog_median_iqr.npz')
    a_names=curves['A_names'].tolist();b_names=curves['B_names'].tolist()
    edges=curves['ring_edges_um'][1:]

    # Approved retained per-area product: no regional aggregation is performed.
    area_curves=np.load(spatial_root/'markcorr_v2/final/cellprog_byarea_median_iqr.npz')
    area_names=area_curves['area_names'].tolist()
    area_counts=area_curves['n_chips_area'].tolist()
    area_a=area_curves['A_names'].tolist();area_b=area_curves['B_names'].tolist()
    area_edges=area_curves['ring_edges_um'][1:]
    area_colors=['#7b9cc7','#d5a06e','#87ad8a','#c98289','#a99bc4','#a78471',
                 '#c595b5','#939ba3','#b4b77c','#7eafb5','#8290b7','#c9b58e',
                 '#83a9a0','#b49d9d']
    area_styles=['-',(0,(4,2)),(0,(1,1.4)),(0,(5,1.5,1,1.5))]
    global_color='#173e5c'

    # Original scatter source: figA_svg_panels.R uses the FINAL TSV's
    # log2(median g), fraction, and set.seed(7) runif display jitter. The TSV is
    # written from this approved NPZ in A-major/B-minor order at .6g/.4f.
    # TableS6 median(log2 g) is a different estimator and is not substituted.
    npoints=len(a_names)*len(b_names)
    rcode=('options(future.globals.maxSize = 50 * 1024^3); '
           'set.seed(7); writeLines(sprintf("%.17g", '
           f'runif({npoints}, -0.012, 0.012)))')
    jitter_text=subprocess.run(['/usr/local/bin/Rscript','--vanilla','-e',rcode],
                               check=True,capture_output=True,text=True).stdout
    jitter=np.fromstring(jitter_text,sep=' ')

    # Retain the original full-axis coordinates and seed-7 jitter first.
    # Panel c alone displays the 54 canonical programs in the current Table S3;
    # filtering after jitter preserves each retained point's original offset.
    retained_b_names={'program_'+row['cnmf_component'] for row in annotation.values()}
    missing_b_names=retained_b_names-set(b_names)
    if len(retained_b_names)!=54 or missing_b_names:
        raise ValueError('Panel c needs the 54 retained Table S3 raw components: '+
                         str(sorted(missing_b_names)))
    retained_b_mask=np.asarray([name in retained_b_names for name in b_names],dtype=bool)
    if int(retained_b_mask.sum())!=54:
        raise ValueError('Panel c raw axis does not contain exactly 54 retained components')
    c_display_mask=np.tile(retained_b_mask,len(a_names))
    # No estimators, tiers, original point order, or coordinates are refitted.
    effects=np.asarray([float(format(float(v),'.6g')) for v in curves['log2_median_g'][:,:,0].ravel()])
    fractions=np.asarray([float(format(float(v),'.4f')) for v in curves['frac_same_sign'][:,:,0].ravel()])
    shown_fractions=fractions+jitter
    tier_colors=np.where(np.abs(effects)>=.32,np.where(effects>0,'#B2182B','#2166AC'),'#CCCCCC')
    axc=axes_mm(112,119.0,55,24) if compact_candidate else axes_mm(116,136,50,22)
    axc.set_axisbelow(True)
    axc.grid(True,color='#e5e5e5',linewidth=.45)
    for boundary in [-.32,.32]:
        axc.axvline(boundary,color='#999999',linewidth=.6,linestyle=(0,(3,3)),zorder=1)
    axc.scatter(effects[c_display_mask],shown_fractions[c_display_mask],
                c=tier_colors[c_display_mask],s=.8,alpha=.7,linewidths=0,zorder=2)
    print('C_DISPLAY_SCOPE','retained_programs',int(retained_b_mask.sum()),
          'subclasses',len(a_names),'display_points',int(c_display_mask.sum()),
          'original_jitter_axis_points',npoints,
          'excluded_components',[name for name in b_names if name not in retained_b_names],flush=True)
    # Bounds are the original PDF's native plotting rectangle, converted using
    # its observed -1/0/+1 and .5/1 grid coordinates, not refitted data limits.
    axc.set_xlim((341.6719-436.8828)/50.6523,(532.0899-436.8828)/50.6523)
    axc.set_ylim(1+(310.1479-386.0352)/141.6026,1+(310.1479-304.7582)/141.6026)
    axc.set_xticks([-1,0,1])
    axc.set_yticks([.5,.6,.7,.8,.9,1.0])
    axc.set_xlabel('log₂ median g ([0,25) µm)' if compact_candidate else r'$\log_2$ median $g$ ([0,25) µm)',fontsize=6 if compact_candidate else 5.6,labelpad=.7 if compact_candidate else 2)
    axc.set_ylabel('Same-direction fraction',fontsize=6 if compact_candidate else 5.8,labelpad=.7 if compact_candidate else 2)
    axc.tick_params(length=0,labelsize=5.05 if compact_candidate else 5.6,pad=.7 if compact_candidate else 2)
    for spine in axc.spines.values():spine.set_visible(False)
    axc.legend(handles=[
        Line2D([],[],marker='o',markersize=1.8,linestyle='',color='#2166AC',label='Negative'),
        Line2D([],[],marker='o',markersize=1.8,linestyle='',color='#B2182B',label='Positive'),
        Line2D([],[],marker='o',markersize=1.8,linestyle='',color='#CCCCCC',label='Below threshold')],
        loc='upper center',bbox_to_anchor=(.5,1.16 if compact_candidate else -.3045),ncol=3,frameon=False,
        fontsize=5.05 if compact_candidate else 5.2,handlelength=.8,columnspacing=.65 if compact_candidate else 1.0,borderaxespad=0)
    highlights=[
        ('d','P49','L6 IT','o',.97,.36,'right'),
        ('e','P52','L6 IT','D',.04,.29,'left'),
        ('f','P47','CHANDELIER','s',.97,.73,'right'),
        ('g','P47','L3-L4 IT RORB','^',.04,.78,'left')]
    for letter,pid,target,marker,tx,ty,align in highlights:
        ai=a_names.index(target);bi=b_names.index('program_'+annotation[pid]['cnmf_component'])
        index=ai*len(b_names)+bi
        if not c_display_mask[index]:
            raise ValueError('Existing highlight is outside the retained panel c scope: '+pid)
        effect=effects[index];fraction=fractions[index];shown_fraction=shown_fractions[index]
        identity=annotation[pid]['dominant_subclass']
        color=CLASS_COLORS[CELL_CLASS[identity]]
        axc.scatter([effect],[shown_fraction],s=2.9**2,marker=marker,facecolors='none',
                    edgecolors=color,linewidths=.65,zorder=5)
        axc.annotate(letter+' '+pid+'\n'+target.replace('L3-L4','L3–L4'),
            xy=(effect,shown_fraction),xycoords='data',xytext=(tx,ty),textcoords='axes fraction',
            ha=align,va='center',fontsize=5.5,fontweight='bold',color=color,zorder=6,
            bbox=dict(facecolor='none' if compact_candidate else 'white',edgecolor='none',pad=.7,alpha=.9),
            arrowprops=dict(arrowstyle='-',color=color,lw=.45,alpha=.9,shrinkA=3,shrinkB=3))
        print('ORIGINAL_SCATTER_CALLOUT',letter,pid,target,
              'log2_median_g',effect,'fraction',fraction,'display_jitter',float(jitter[index]),flush=True)
    def curve(x,y,w,h, selections):
        ax=axes_mm(x,y,w,h)
        ax.axhline(0,color='#999999',linewidth=.65,zorder=0)
        for pid,target,color,label in selections:
            raw=int(annotation[pid]['cnmf_component'])
            ai=a_names.index(target);bi=b_names.index('program_'+str(raw))
            middle=curves['log2_median_g'][ai,bi]
            lower=curves['log2_q1'][ai,bi];upper=curves['log2_q3'][ai,bi]
            ax.fill_between(edges,lower,upper,color=global_color,alpha=.055,linewidth=0,zorder=1)
            area_ai=area_a.index(target);area_bi=area_b.index('program_'+str(raw))
            # Plot every stored area, including original n=1/2 areas and all
            # extrema. Matplotlib's data bounds include all fourteen lines.
            for ri,region in enumerate(area_names):
                regional=area_curves['log2_median_g'][ri,area_ai,area_bi]
                ax.plot(area_edges,regional,color=area_colors[ri],
                        linestyle=area_styles[ri%len(area_styles)],
                        linewidth=.45,alpha=.65,zorder=2)
            # This is the original all-44-chip median, not a mean of area curves.
            ax.plot(edges,middle,'o-',color=global_color,linewidth=1.1,
                    markersize=1.8,label=label,zorder=4)
        ax.set_xlim(15,510)
        ax.set_xticks([25,250,500])
        ax.set_xlabel('Distance (µm)',fontsize=6 if compact_candidate else 5.8,labelpad=.7 if compact_candidate else 1.4)
        ax.set_ylabel('log₂ median g' if compact_candidate else r'$\log_2$ median $g$',fontsize=6 if compact_candidate else 5.6,labelpad=.7 if compact_candidate else 1.0)
        from matplotlib.ticker import MaxNLocator, FuncFormatter
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4,steps=[1,2,5,10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos: f'{value:g}'.replace('-','−')))
        ax.tick_params(labelsize=5.05 if compact_candidate else 5.6,length=1 if compact_candidate else 1.4,width=.45,pad=.7 if compact_candidate else 1.5)
        for side in ['top','right']:ax.spines[side].set_visible(False)
        # The shared region/all-44 legend names these same four curves.
    # Owner-locked V1 section for all spatial displays; join by bin identity, not row position.
    import pyarrow.parquet as pq
    chips=['D00865B3']
    meta=pq.read_table(spatial_root/'spatial_bin50_meta.parquet',
         filters=[('chip','in',chips)],columns=['bin','chip','x','y','region']).to_pandas()
    weights=pq.read_table(spatial_root/'spatial_bin50_rctd_weights.parquet',
         filters=[('chip','in',chips)],
         columns=['bin','L6 IT','CHANDELIER','L3-L4 IT RORB','OLIGO','OPC','rctd_pass_mask']).to_pandas()
    score=pq.read_table(spatial_root/'spatial_bin50_program_score_SCT.parquet',
         columns=['bin','bin_total_umi','program_58','program_51','program_54','program_45','program_43']).to_pandas()
    score=score[score['bin'].isin(meta['bin'])]
    maps=meta.merge(weights,on='bin',how='left').merge(score,on='bin',how='left')
    maps['valid']=maps['rctd_pass_mask'].fillna(False).astype(bool)&(maps['bin_total_umi'].fillna(0)>=200)

    # A single continuous blue scale for spatial fields; area-line colors are
    # categorical and shared across d-g, while the all-44 line stays dark blue.
    # Normalization remains per original chip/field; each quantity is named.
    continuous=LinearSegmentedColormap.from_list('spatial_blue',['#f3f7fb','#bdd5e9','#629ec7','#245c88','#10395c'])
    field_cache={}
    def tissue(chip,field,x,y,w,h,title,identity,program=False,title_size=6.2,tick_size=5.5,label_size=5.5,compact_key=True,title_offset=6.0):
        key=(chip,field)
        if key not in field_cache:
            sub=maps[maps['chip']==chip]
            values=sub[field].to_numpy(dtype=float)
            valid=sub['valid'].to_numpy()&np.isfinite(values)
            lo,hi=np.percentile(values[valid],[1,99]);lo=lo if lo<0 else 0.0
            if hi<=lo:hi=lo+1e-6
            display=np.clip((values-lo)/(hi-lo),0,1)
            xx=sub['x'].to_numpy(dtype=float)/2;yy=sub['y'].to_numpy(dtype=float)/2
            field_cache[key]=(xx,yy,valid,display)
        # In particular, both P47 comparisons reuse this identical cached field.
        xx,yy,valid,display=field_cache[key]
        ax=axes_mm(x,y,w,h)
        original_height=61.0 if chip=='C00841F3' else 49.0
        area=(.18*72/25.4*(h/original_height))**2
        ax.scatter(xx,yy,s=area,c='#dcdcdc',marker='s',linewidths=0,rasterized=True)
        ax.scatter(xx[valid],yy[valid],s=area,c=display[valid],vmin=0,vmax=1,
                   cmap=continuous,marker='s',linewidths=0,rasterized=True)
        ax.set_xlim(xx.min()-12.5,xx.max()+12.5);ax.set_ylim(yy.max()+12.5,yy.min()-12.5)
        ax.set_aspect('equal',adjustable='box')
        if compact_candidate:
            ax.set_box_aspect((yy.max()-yy.min()+25)/(xx.max()-xx.min()+25))
        ax.set_axis_off()
        # x and y have already been converted from DNB px to micrometres
        # above. Each field gets its own real 1000-um (=1 mm) scale bar,
        # placed in a reserved lane above the tissue, apart from its color key.
        bar_transform=blended_transform_factory(ax.transData,ax.transAxes)
        xmin,xmax=ax.get_xlim()
        bar_x=xmin+.05*(xmax-xmin)
        bar_y=1+.7/h if compact_candidate else 1.065
        ax.add_line(Line2D([bar_x,bar_x+1000],[bar_y,bar_y],
            transform=bar_transform,color='#202020',linewidth=1.2,
            solid_capstyle='butt',clip_on=False,zorder=12))
        ax.text(bar_x+1000+.025*(xmax-xmin),bar_y,'1 mm',
            transform=bar_transform,ha='left',va='center',fontsize=5.05 if compact_candidate else 5.0,
            color='#202020',clip_on=False,zorder=12)
        # Categorical title accents use the same hue as the identity bands,
        # but a light tint and dark text keep the scientific fields primary.
        from matplotlib.colors import to_rgb
        identity_rgb=to_rgb(CLASS_COLORS[CELL_CLASS[identity]])
        title_background=tuple(.84+.16*v for v in identity_rgb)
        title_foreground=tuple(.62*v for v in identity_rgb)
        if compact_candidate:
            if program:
                pid=next(k for k,v in annotation.items() if 'program_'+v['cnmf_component']==field)
                identity_label=DISPLAY.get(identity,identity)
                name=annotation[pid]['functional_name']
                if pid=='P49':name='Microglial complement/MHC'
                elif pid=='P52':name='Neuropeptide signaling (inh)'
                elif pid=='P47':name='Reactive astrocyte/\nvascular'
                text_mm(x+w/2,y-9.1,identity_label+' · '+pid,ha='center',fontsize=7,
                        color=title_foreground,fontweight='bold',zorder=10)
                text_mm(x+w/2,y-6.1,name,ha='center',fontsize=5.05 if '\n' in name else 5.5,
                        color=title_foreground,fontweight='bold',linespacing=1.0,zorder=10)
            else:
                text_mm(x+w/2,y-9.1,'Target cell',ha='center',fontsize=5.5,
                        color=title_foreground,zorder=10)
                text_mm(x+w/2,y-5.7,DISPLAY.get(identity,identity).replace('L3-L4','L3–L4'),
                        ha='center',fontsize=6 if identity=='L3-L4 IT RORB' else 7,color=title_foreground,fontweight='bold',zorder=10)
        else:
            text_mm(x+w/2,y-title_offset,title,ha='center',va='center',fontsize=title_size,fontweight='bold',
                    linespacing=1.08,color=title_foreground,zorder=10,
                    bbox=dict(facecolor=title_background,edgecolor='none',boxstyle='square,pad=.18'))
        cax=axes_mm(x+w*(.08 if compact_key else .12),y+h+(.25 if compact_candidate else .5),
                    w*(.84 if compact_key else .76),.35 if compact_candidate else .6)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=continuous),cax=cax,
                       orientation='horizontal',ticks=[0,1])
        cb.ax.set_xticklabels(['Low','High']);cb.ax.tick_params(length=0,labelsize=5.05 if compact_candidate else tick_size,pad=.65)
        cb.outline.set_visible(False)
        text_mm(x+w/2,y+h+(.95 if compact_candidate else (1.45 if compact_key else 3.6)),
                ('Program score' if program else 'Cell weight') if compact_key
                else ('Relative program score' if program else 'Relative cell weight'),
                ha='center',fontsize=5.05 if compact_candidate else label_size)
        print('SPATIAL_FIELD',chip,field,flush=True)
        return ax
    for pid,raw in [('P49',54),('P52',58),('P47',51),('P41',45),('P39',43)]:
        if int(annotation[pid]['cnmf_component'])!=raw:
            raise ValueError('The approved component identity changed for '+pid)
    # Each complete relation has exactly two continuous fields plus its own curve.
    # Header colors identify reference preference / target subclass, not magnitude.
    relation_groups=[
        ('d','D00865B3','V1','P49','program_54','L6 IT',175.5,
         'Microglial\ncomplement/MHC (P49)'),
        ('e','D00865B3','V1','P52','program_58','L6 IT',175.5,
         'Neuropeptide\nsignaling (inh, P52)'),
        ('f','D00865B3','V1','P47','program_51','CHANDELIER',209.0,
         'Reactive astrocyte/\nvascular (P47)'),
        ('g','D00865B3','V1','P47','program_51','L3-L4 IT RORB',209.0,
         'Reactive astrocyte/\nvascular (P47)')]
    for letter,chip,region,pid,field,target,yy,program_title in relation_groups:
        # Two categorical segments, not an arrow or pooled-class analysis:
        # upper = program reference-preferred class; lower = target class.
        class_pair=[CELL_CLASS[annotation[pid]['dominant_subclass']],CELL_CLASS[target]]
        group_x=0 if letter in ['d','f'] else (85.0 if compact_candidate else 84.7)
        if compact_candidate:
            yy=159.0 if letter in ['d','e'] else 195.0
            side_x,side_y,side_w,side_h=group_x+1.0,yy,1.8,23.2
        else:
            side_x,side_y,side_w,side_h=group_x+3.0,yy-5.8,2.2,29.2
        for segment,broad in enumerate(class_pair):
            top=side_y+segment*side_h/2
            fig.patches.append(Rectangle((side_x/W,1-(top+side_h/2)/H),
                side_w/W,(side_h/2)/H,transform=fig.transFigure,
                facecolor=CLASS_COLORS[broad],edgecolor='white',linewidth=.45,zorder=7))
            fg='#222222' if broad in ['Micro','VLMC'] else 'white'
            text_mm(side_x+side_w/2,top+side_h/4,broad,rotation=90,ha='center',va='center',
                    fontsize=6.2,fontweight='bold',color=fg,zorder=8)
        print('DISPLAY_CLASS_PAIR',letter,class_pair,flush=True)
        main_text=dict(title_size=6.2,tick_size=5.5,label_size=5.5,compact_key=True,title_offset=4.7)
        map_w=23.2 if compact_candidate else 24
        map_h=23.2 if compact_candidate else 21.4
        tissue(chip,field,group_x+(4.6 if compact_candidate else 6.7),yy,map_w,map_h,program_title,
               annotation[pid]['dominant_subclass'],program=True,**main_text)
        target_title='Target cell\n'+target.replace('L3-L4','L3–L4')
        tissue(chip,target,group_x+(30.6 if compact_candidate else 31.5),yy,map_w,map_h,target_title,target,**main_text)
        # Equal x/y physical units; square V1 fields are never stretched to
        # fill a non-square grid slot. Curves reflow without rescaling data.
        curve_w,curve_h=(20.6,23.5) if compact_candidate else (19.2,21.4)
        curve(group_x+(62 if compact_candidate else 64),yy-(3 if compact_candidate else 2.2),
              curve_w,curve_h,[(pid,target,'#245c88',pid+' – '+target)])
        print('HORIZONTAL_RELATION',letter,'curve_mm',[curve_w,curve_h],
              'field_box_mm',[map_w,map_h],'chip',chip,flush=True)

    # b: two anatomical-reference pairs in one horizontal row, both shown on
    # the owner-locked V1 section, not new mechanistic findings.
    # Only continuous fields are shown, with the original >=200 display mask
    # and per-chip/field q1/q99 recipe. No overlay thresholds or new curves.
    text_mm(99 if compact_candidate else 4.0,113.8 if compact_candidate else 132.4,'V1',fontsize=5.2)
    reference_text=dict(title_size=5.8,tick_size=5.0,label_size=5.2,compact_key=True,title_offset=5.2)
    by,bsize=(123.0,23.2) if compact_candidate else (141.0,21.2)
    bx=[2.2,27.0,53.2,79.0] if compact_candidate else [4.0,29.0,57.0,82.0]
    tissue('D00865B3','program_45',bx[0],by,bsize,bsize,'Myelination\n(P41)',
           annotation['P41']['dominant_subclass'],program=True,**reference_text)
    tissue('D00865B3','OLIGO',bx[1],by,bsize,bsize,'Target cell\nOligo','OLIGO',**reference_text)
    tissue('D00865B3','program_43',bx[2],by,bsize,bsize,'ECM/OPC proteoglycan\n(P39)',
           annotation['P39']['dominant_subclass'],program=True,**reference_text)
    tissue('D00865B3','OPC',bx[3],by,bsize,bsize,'Target cell\nOPC','OPC',**reference_text)

    # One shared area key for all four distance panels; counts are source
    # metadata. No source 'unstable' area is silently omitted.
    region_handles=[Line2D([],[],color=area_colors[i],
        linestyle=area_styles[i%len(area_styles)],linewidth=.65,
        label=f'{region} {int(area_counts[i])}' if compact_candidate else f'{region} (n={int(area_counts[i])})') for i,region in enumerate(area_names)]
    region_handles.append(Line2D([],[],color=global_color,linewidth=1.1,label='All 44 (IQR)'))
    text_mm(2 if compact_candidate else 3,100.4 if compact_candidate else 119.5,
            'd–g: regions (n)' if compact_candidate else 'd–g: regions (n sections)',fontsize=5.5 if compact_candidate else 5.3)
    fig.legend(handles=region_handles,loc='upper center',
        bbox_to_anchor=((22.7 if compact_candidate else 28)/W,1-(103 if compact_candidate else 121.6)/H),ncol=4,frameon=False,
        fontsize=5.05 if compact_candidate else 4.8,handlelength=1.1,handletextpad=.3,columnspacing=.4 if compact_candidate else .45,
        labelspacing=.05 if compact_candidate else .2,borderaxespad=0,borderpad=0)
    # This legend describes a display intersection with the unchanged original
    # support mask, not an additional statistical threshold or confidence star.
    frame_key_x,frame_key_y=(2.0,111.5) if compact_candidate else (3.0,130.2)
    fig.patches.append(Rectangle((frame_key_x/W,1-(frame_key_y+1.2)/H),1.7/W,1.2/H,
        transform=fig.transFigure,facecolor='none',edgecolor='black',linewidth=.65,zorder=7))
    text_mm(frame_key_x+2.4,frame_key_y-.2,'Original support; cross-subclass',fontsize=5.05)
    print('RETAINED_AREA_LINES',[(name,int(n)) for name,n in zip(area_names,area_counts)],flush=True)

    # h: read the original five relations and six displayed regions.
    # Owner-approved display only: centre each pair within those six regions
    # and divide by its sample SD (ddof=1), equivalent to R scale.
    # The original effect product, a, spatial scores and distance curves are
    # never overwritten or standardised by this display transformation.
    h_rows=[
        ('P39','OPC','OPC'),
        ('P28','L6B','L6B'),
        ('P19','L6 CT','L6 CT'),
        ('P30','PAX6','VIP'),
        ('P49','MICRO','L6 IT')]
    h_regions=['DLPFC','S1','AG','M1','V1','SMG']
    h_original=np.asarray([
        [area_curves['log2_median_g'][
            area_names.index(region),area_a.index(target),
            area_b.index('program_'+annotation[pid]['cnmf_component']),0]
         for region in h_regions]
        for pid,reference,target in h_rows],dtype=float)
    h_means=h_original.mean(axis=1,keepdims=True)
    h_sd=h_original.std(axis=1,ddof=1,keepdims=True)
    invalid_rows=(~np.isfinite(h_original).all(axis=1)) | (~np.isfinite(h_sd[:,0])) | (h_sd[:,0]<=0)
    if invalid_rows.any():
        invalid_pairs=[h_rows[i] for i in np.flatnonzero(invalid_rows)]
        raise ValueError('Within-pair display Z-score has missing or zero-SD source values: '+str(invalid_pairs))
    h_values=(h_original-h_means)/h_sd
    h_min=float(h_values.min());h_max=float(h_values.max())
    from matplotlib.colors import TwoSlopeNorm
    h_norm=TwoSlopeNorm(vmin=h_min,vcenter=0,vmax=h_max)
    h_cmap=LinearSegmentedColormap.from_list('within_pair_z',
                ['#2166AC','#F7F7F7','#B2182B'])
    hy,hrow_y,hrow_h=(221.5,224.1,1.75) if compact_candidate else (236.5,239.5,2.05)
    hname,href,htarget,hband_w,hmat,hcell=(3.8,44.0,60.0,14.0,77.0,10.4) if compact_candidate else (4,60,82.5,17.5,104,7.4)
    text_mm(hname,hy,'Program',fontsize=5.5 if compact_candidate else 5.8,fontweight='bold')
    text_mm(href+hband_w/2,hy,'Ref. pref.' if compact_candidate else 'Ref. preference',ha='center',fontsize=5.5 if compact_candidate else 5.6,fontweight='bold')
    text_mm(htarget+hband_w/2,hy,'Target',ha='center',fontsize=5.5 if compact_candidate else 5.6,fontweight='bold')
    for j,region in enumerate(h_regions):
        text_mm(hmat+(j+.5)*hcell,hy,region,ha='center',fontsize=5.05 if compact_candidate else 5.5)
    for i,(pid,reference,target) in enumerate(h_rows):
        yy=hrow_y+(i+.5)*hrow_h
        text_mm(hname,yy,annotation[pid]['functional_name']+' ('+pid+')',
                va='center',fontsize=5.05 if compact_candidate else 5.8)
        for xx,identity in [(href,reference),(htarget,target)]:
            broad=CELL_CLASS[identity]
            band_h=hrow_h-.1 if compact_candidate else 1.92
            fig.patches.append(Rectangle((xx/W,1-(yy+band_h/2)/H),hband_w/W,band_h/H,
                transform=fig.transFigure,facecolor=CLASS_COLORS[broad],
                edgecolor='none',zorder=7))
            fg='#222222' if broad in ['Micro','VLMC'] else 'white'
            text_mm(xx+hband_w/2,yy,DISPLAY.get(identity,identity),ha='center',va='center',
                    fontsize=5.05 if compact_candidate else 5.8,color=fg,zorder=8)
        for j,value in enumerate(h_values[i]):
            fig.patches.append(Rectangle(((hmat+j*hcell)/W,1-(hrow_y+(i+1)*hrow_h)/H),
                hcell/W,hrow_h/H,transform=fig.transFigure,
                facecolor=h_cmap(h_norm(value)),edgecolor='white',linewidth=.18,zorder=7))
    h_cax=axes_mm(154.5,225.9,3,6.0) if compact_candidate else axes_mm(156,241.3,3,7.4)
    h_ticks=[h_min,0.0,h_max]
    h_cb=fig.colorbar(plt.cm.ScalarMappable(norm=h_norm,cmap=h_cmap),
                     cax=h_cax,orientation='vertical',ticks=h_ticks)
    h_cb.ax.set_yticklabels([f'{value:.3f}' if value else '0' for value in h_ticks])
    h_cb.ax.tick_params(labelsize=5.05 if compact_candidate else 5.2,length=1,pad=1)
    h_cb.outline.set_visible(False)
    if compact_candidate:
        text_mm(144,221.5,'Within-pair Z-score',fontsize=5.5)
    else:
        text_mm(152,236.0,'Within-pair\nZ-score',fontsize=5.3,linespacing=1.05)
    print('REGIONAL_H_DISPLAY','within-pair Z-score across the six displayed regions',
          'sample_sd_ddof',1,'display_min',h_min,'display_max',h_max,
          'original_field','log2_median_g','first_ring_um',25,flush=True)

    area_curves.close()
    curves.close()

    overlay=io.BytesIO();fig.savefig(overlay,format='pdf',dpi=300,transparent=True);plt.close(fig);overlay.seek(0)
    writer=PdfWriter();page=writer.add_blank_page(width=W*MM,height=H*MM)
    # One native Form XObject per source variant; crops reuse its original
    # painting rather than reparsing/copying the entire page for every crop.
    xobjects=DictionaryObject(); form_names={}; commands=[]
    for source,box,x,y,width,clockwise,display_height in crops:
        identity=id(source)
        if identity not in form_names:
            name=NameObject('/Native'+str(len(form_names)))
            form=DecodedStreamObject()
            form[NameObject('/Type')]=NameObject('/XObject')
            form[NameObject('/Subtype')]=NameObject('/Form')
            form[NameObject('/BBox')]=RectangleObject(source.mediabox)
            form[NameObject('/Resources')]=source['/Resources'].clone(writer)
            form.set_data(source.get_contents().get_data())
            xobjects[name]=writer._add_object(form.flate_encode())
            form_names[identity]=name
        name=form_names[identity]
        x0,y0,x1,y1=box
        bottom=float(source.mediabox.height)-y1
        if clockwise:
            height=width*(x1-x0)/(y1-y0)
            scale=width*MM/(y1-y0)
            dx=x*MM;dy=(H-y-height)*MM
            transform=f'0 {-scale:.10f} {scale:.10f} 0 {dx-bottom*scale:.8f} {dy+height*MM+x0*scale:.8f}'
        else:
            height=width*(y1-y0)/(x1-x0) if display_height is None else display_height
            scale=width*MM/(x1-x0)
            scale_y=height*MM/(y1-y0)
            dx=x*MM;dy=(H-y-height)*MM
            transform=f'{scale:.10f} 0 0 {scale_y:.10f} {dx-x0*scale:.8f} {dy-bottom*scale_y:.8f}'
        commands.append(f'q {dx:.8f} {dy:.8f} {width*MM:.8f} {height*MM:.8f} re W n '
                        f'{transform} cm {name} Do Q\n')
    page[NameObject('/Resources')]=DictionaryObject({NameObject('/XObject'):xobjects})
    native_stream=DecodedStreamObject();native_stream.set_data(''.join(commands).encode('ascii'))
    page[NameObject('/Contents')]=writer._add_object(native_stream.flate_encode())
    page.merge_page(PdfReader(overlay).pages[0])
    page.compress_content_streams()
    output_stem='Fig4_compact_candidate' if compact_candidate else 'Fig4_reorganized'
    pdf_path=OUT/(output_stem+'.pdf');png_path=OUT/(output_stem+'.png')
    with pdf_path.open('wb') as handle:writer.write(handle)
    print('RETAINED',pdf_path,flush=True)
    png_command=(['pdftocairo','-png','-transp','-singlefile','-r','300']
                 if compact_candidate else ['pdftoppm','-singlefile','-png','-r','300'])
    subprocess.run(png_command+[str(pdf_path),str(png_path.with_suffix(''))],check=True)
    print('RETAINED',png_path,flush=True)
    if compact_candidate:
        print('CANDIDATE_LAYOUT','page_mm',[W,H],'dense_min_pt',5.05,
              'spatial_aspect','equal x/y data units; square V1 field',
              'background','native transparency; scientific whites unchanged',flush=True)
    print('CONTENT a clockwise full matrix with retained abbreviations/IT display; b four V1 reference fields in one row; c original-coordinate scatter with updated e=P52/L6 IT; d-g aligned programme/target/near-square-curve triptychs with lighter identity title tints; all twelve spatial fields use V1 D00865B3 with separate real 1-mm scale bars and program score labels; fourteen regional curves and original all-44 median/IQR unchanged; h within-pair six-region display Z-score with sample SD and actual display min/max around zero; original effects untouched; S13 untouched',flush=True)


def draw_removed_supplement(base_pdf):
    """One retained S13: authentic original support once, plus the count overview."""
    import io,subprocess
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from pypdf import PdfReader,PdfWriter
    from pypdf.generic import ContentStream,DecodedStreamObject,NameObject,RectangleObject,DictionaryObject
    W,H,MM=390.0,535.0,72/25.4
    original_reader=PdfReader(io.BytesIO(base_pdf));original=original_reader.pages[0]
    reader_holders=[original_reader]
    def mul(a,b):
        return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
                a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3],
                a[4]*b[0]+a[5]*b[2]+b[4],a[4]*b[1]+a[5]*b[3]+b[5])
    text_ops={b'BT',b'ET',b'Tf',b'Tm',b'Td',b'TD',b'T*',b'Tj',b'TJ',b'Tc',b'Tw',b'Tz',b'TL',b'Ts',b'Tr',b"'",b'"'}
    def variant(no_text=False,no_rotated=False,no_orange=False,no_example_labels=False):
        reader=PdfReader(io.BytesIO(base_pdf)); page=reader.pages[0]; reader_holders.append(reader)
        def orange(color):
            return len(color)==3 and all(abs(x-y)<.004 for x,y in zip(color,(.85098,.37255,.00784)))
        def rewrite(stream,initial=(1,0,0,1,0,0)):
            cs=ContentStream(stream,reader); result=[]; stack=[]
            ctm=initial; stroke=(0,0,0); fill=(0,0,0); buffer=None; rotated=False
            path_points=[]; path_closed=False; text_position=None
            for operands,operator in cs.operations:
                if operator==b'q': stack.append((ctm,stroke,fill))
                elif operator==b'Q' and stack: ctm,stroke,fill=stack.pop()
                elif operator==b'cm': ctm=mul(tuple(map(float,operands)),ctm)
                elif operator==b'RG': stroke=tuple(map(float,operands))
                elif operator==b'rg': fill=tuple(map(float,operands))
                elif operator==b'G': stroke=(float(operands[0]),)*3
                elif operator==b'g': fill=(float(operands[0]),)*3
                if operator in (b'm',b'l',b'c',b'v',b'y'):
                    if operator==b'm': path_points=[]; path_closed=False
                    vals=list(map(float,operands))
                    path_points.extend((vals[i]*ctm[0]+vals[i+1]*ctm[2]+ctm[4],
                                        721-(vals[i]*ctm[1]+vals[i+1]*ctm[3]+ctm[5])) for i in range(0,len(vals),2))
                elif operator==b're':
                    path_closed=True
                elif operator==b'h': path_closed=True
                if operator==b'BT': buffer=[(operands,operator)]; rotated=False; text_position=None; continue
                if buffer is not None:
                    buffer.append((operands,operator))
                    if operator==b'Tm':
                        tm=tuple(map(float,operands)); dx=tm[0]*ctm[0]+tm[1]*ctm[2]; dy=tm[0]*ctm[1]+tm[1]*ctm[3]
                        rotated=abs(dy)>3*abs(dx)
                        text_position=(tm[4]*ctm[0]+tm[5]*ctm[2]+ctm[4],721-(tm[4]*ctm[1]+tm[5]*ctm[3]+ctm[5]))
                    example_label=no_example_labels and text_position is not None and ((250<text_position[0]<303 and 519<text_position[1]<530) or (297<text_position[0]<338 and 524<text_position[1]<538))
                    if (example_label or (no_rotated and rotated)) and operator in (b'Tj',b'TJ',b"'",b'"'):
                        buffer.pop()
                    if False and no_rotated and rotated and operator in (b'Tj',b'TJ',b"'",b'"'):
                        buffer.pop()
                    if operator==b'ET':
                        if no_text:
                            result.extend(item for item in buffer if item[1] not in text_ops)
                        else: result.extend(buffer)
                        buffer=None
                    continue
                remove = no_orange and (
                    (operator in (b'S',b's') and orange(stroke)) or
                    (operator in (b'f',b'F',b'f*') and orange(fill)) or
                    (operator in (b'B',b'B*',b'b',b'b*') and (orange(fill) or orange(stroke))))
                if no_example_labels and operator in (b'S',b's') and not path_closed and len(path_points)>=2 and max(stroke)<.30:
                    x0=min(v[0] for v in path_points); x1=max(v[0] for v in path_points)
                    y0=min(v[1] for v in path_points); y1=max(v[1] for v in path_points)
                    remove=remove or (340<x0<x1<539 and 299<y0<y1<390 and x1-x0>1 and y1-y0>1)
                result.append(([],b'n') if remove else (operands,operator))
                if operator in (b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*',b'n'):
                    path_points=[]; path_closed=False
            cs.operations=result
            return cs
        seen=set()
        def forms(resources):
            objects=resources.get('/XObject',{})
            if hasattr(objects,'get_object'):objects=objects.get_object()
            for name,ref in list(objects.items()):
                obj=ref.get_object()
                if obj.get('/Subtype')!='/Form' or id(obj) in seen:continue
                seen.add(id(obj)); child=obj.get('/Resources',resources)
                if hasattr(child,'get_object'):child=child.get_object()
                forms(child); replacement=DecodedStreamObject()
                for key,value in obj.items():
                    if key not in ('/Length','/Filter','/DecodeParms'):replacement[key]=value
                replacement.set_data(rewrite(obj,tuple(map(float,obj.get('/Matrix',[1,0,0,1,0,0])))).get_data())
                objects[name]=replacement
        forms(page['/Resources']); page[NameObject('/Contents')]=rewrite(page.get_contents())
        return page
    original=variant(no_example_labels=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.linewidth':.55})
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='none')
    def axes_mm(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
    def text_mm(x,y,text,**kwargs):
        opts=dict(ha='left',va='top',fontsize=8);opts.update(kwargs)
        return fig.text(x/W,1-y/H,text,**opts)
    def white_mm(x,y,w,h):
        fig.patches.append(Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,
                                     facecolor='white',edgecolor='none',zorder=6))
    crops=[]
    def native(page,box,x,y,width):
        height=width*(box[3]-box[1])/(box[2]-box[0]);crops.append((page,box,x,y,width));return height
    scale=380/539;ox,oy=5,10
    native(original,(0,299,539,721),ox,oy,380)
    def erase(box):
        x0,y0,x1,y1=box
        white_mm(ox+x0*scale,oy+(y0-299)*scale,(x1-x0)*scale,(y1-y0)*scale)
    # The sole consistency display stays in the main figure, not duplicated here.
    erase((328,299,539,416))
    erase((320,318.9,328.1,385.2))
    # Recreate only the three J text fragments touched by the neighboring axis.
    for sx,sy,txt in [(314.77,325.3,'J=0.32'),(314.59,334.6,'J=0.39'),(312.13,343.8,'J=0.22')]:
        erase((sx-.3,sy-3,327.2,sy+3))
        text_mm(ox+sx*scale,oy+(sy-299)*scale,txt,fontsize=8.0,va='center',color='#666666',zorder=9)
    # Original panel letters and a tiny preceding label fragment only.
    for box in [(206.8,299,214,300.5),(0,407,5,420),(212,406,220,421),
                (213,526,221,541),(354,531,359,546),(0,588,5,602),
                (133,594,140,608),(280.5,601,286,616)]:
        erase(box)
    # Original V1 overlay headings abut the old consistency legend: reset labels,
    # not image marks, after excluding that old legend.
    erase((315,415.8,539,431))
    text_mm(ox+356*scale,oy+(416-299)*scale,'Oligo weight +\nmyelination-high region (P41)',
            ha='center',fontsize=9.1,zorder=9,linespacing=1.04)
    text_mm(ox+489*scale,oy+(416-299)*scale,'L2–L3 IT weight +\nmyelination-high region (P41)',
            ha='center',fontsize=9.1,zorder=9,linespacing=1.04)
    # The old merged-upper-IT figure contains no high-zone overlay in its actual
    # native display, despite an old caption. Do not repeat that incorrect label.
    erase((251,519,302,528))
    text_mm(ox+318*scale,oy+(524-299)*scale,'Upper-IT weight\n(L2–L3 + L3–L4)',
            ha='center',fontsize=8.0,zorder=9,linespacing=1.02)
    # Keep the sensitivity legend clear of the adjoining overlay headings.
    erase((248,410,318,421))
    from matplotlib.lines import Line2D
    for sx,mark,color,label in [(249,'o','#555555','Full'),(273,'^','#b31636','Retained'),(298,'^','#de9700','Weakened')]:
        xx=ox+sx*scale;yy=oy+(419-299)*scale
        fig.lines.append(Line2D([xx/W],[1-yy/H],transform=fig.transFigure,
                         marker=mark,color=color,markersize=4,linestyle='',zorder=9))
        text_mm(xx+2.5,yy,label,va='center',fontsize=7.2,zorder=9)
    # Current functional names/IDs for the four original regional references.
    for box in [(14,602,69,613),(75,602,130,613),(14,649.5,69,658),(75,649.5,130,658)]:
        erase(box)
    for sx,sy,label in [(42,603,'Ast · Glutamate\ntransport (P14)'),(104,603,'Oligo · Oligodendrocyte/\nmyelin (P33)'),
                       (42,650,'PVALB · Cation\nchannels (P7)'),(104,650,'SST · Glutamate-receptor\nsignaling (P17)')]:
        text_mm(ox+sx*scale,oy+(sy-299)*scale,label,ha='center',fontsize=7.7,zorder=9,linespacing=1.02)
    for letter,x,y in [('a',5,3),('b',84,3),('c',153,3),('d',246,3),
                       ('e',5,86),('f',227,86),('g',155,169),('h',256,176),
                       ('i',5,215),('j',100,219),('k',204,223),('l',5,323)]:
        text_mm(x,y,letter,fontsize=13,fontweight='bold',zorder=10)
    # P30/VIP source is supporting evidence, not another large main example.
    svg=Path((__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/scripts/figmarkcorr_A/svg_panels/figA_n.svg"))
    add_svg_spatial(fig,axes_mm,text_mm,svg,246,26,133,
        ['AG · VIP\nweight','Acetylcholine-receptor\nsignaling (P30)','VIP weight +\nP30-high region (blue)'])
    count_reader=PdfReader(str(OUT/'Fig4_cross_cell_preferences.pdf'))
    reader_holders.append(count_reader);count_page=count_reader.pages[0]
    native(count_page,(0,0,float(count_page.mediabox.width),float(count_page.mediabox.height)),77.5,327,235)
    overlay=io.BytesIO();fig.savefig(overlay,format='pdf',dpi=300,transparent=True);plt.close(fig);overlay.seek(0)
    writer=PdfWriter();page=writer.add_blank_page(width=W*MM,height=H*MM)
    # One native Form XObject per source variant; crops reuse its original
    # painting rather than reparsing/copying the entire page for every crop.
    xobjects=DictionaryObject(); form_names={}; commands=[]
    for source,box,x,y,width in crops:
        identity=id(source)
        if identity not in form_names:
            name=NameObject('/Native'+str(len(form_names)))
            form=DecodedStreamObject()
            form[NameObject('/Type')]=NameObject('/XObject')
            form[NameObject('/Subtype')]=NameObject('/Form')
            form[NameObject('/BBox')]=RectangleObject(source.mediabox)
            form[NameObject('/Resources')]=source['/Resources'].clone(writer)
            form.set_data(source.get_contents().get_data())
            xobjects[name]=writer._add_object(form.flate_encode())
            form_names[identity]=name
        name=form_names[identity]
        x0,y0,x1,y1=box;height=width*(y1-y0)/(x1-x0)
        bottom=float(source.mediabox.height)-y1
        scale=width*MM/(x1-x0)
        dx=x*MM;dy=(H-y-height)*MM
        commands.append(f'q {dx:.8f} {dy:.8f} {width*MM:.8f} {height*MM:.8f} re W n '
                        f'{scale:.10f} 0 0 {scale:.10f} {dx-x0*scale:.8f} {dy-bottom*scale:.8f} cm {name} Do Q\n')
    page[NameObject('/Resources')]=DictionaryObject({NameObject('/XObject'):xobjects})
    native_stream=DecodedStreamObject();native_stream.set_data(''.join(commands).encode('ascii'))
    page[NameObject('/Contents')]=writer._add_object(native_stream.flate_encode())
    page.merge_page(PdfReader(overlay).pages[0])
    page.compress_content_streams()
    pdf_path=OUT/'FigS13.pdf';png_path=OUT/'FigS13.png'
    with pdf_path.open('wb') as handle:writer.write(handle)
    print('RETAINED',pdf_path,flush=True)
    subprocess.run(['pdftoppm','-singlefile','-png','-r','300',str(pdf_path),str(png_path.with_suffix(''))],check=True)
    print('RETAINED',png_path,flush=True)
    print('CONTENT a/b/c original b/c/d; d P30-VIP; e/f/g/h original f/g/h/i; i/j/k original j/k/l; l unchanged directed preference-to-target program counts; no duplicated consistency or intermediate profiles',flush=True)

def display_labels_only():
    """Relabel only the retained native Fig4 legend; no data are loaded."""
    import fitz
    path = PROJECT_ROOT / 'source_figure_pdfs/main_figures_pdf/Fig4_reorganized.pdf'
    replacements = {
        'headline −': 'Negative',
        'headline +': 'Positive',
        'Negative association': 'Negative',
        'Positive association': 'Positive',
        'sub-threshold': 'Below threshold',
    }
    doc = fitz.open(path)
    page = doc[0]
    pending = []
    legend_lines = []
    for block in page.get_text('dict').get('blocks', []):
        if block.get('type') != 0:
            continue
        for line in block.get('lines', []):
            line_text = ''.join(span.get('text', '') for span in line.get('spans', []))
            if ('Positive associationBelow threshold' in line_text or
                    'Positive association Below threshold' in line_text or
                    any(key in line_text for key in ('headline −', 'headline +', 'sub-threshold'))):
                legend_lines.append(line)
                continue
            for span in line.get('spans', []):
                text = span.get('text', '')
                if text in replacements:
                    pending.append((span, replacements[text]))
    if legend_lines:
        spans = [span for line in legend_lines for span in line.get('spans', [])]
        legend_rect = fitz.Rect(legend_lines[0]['bbox'])
        for line in legend_lines[1:]:
            legend_rect |= fitz.Rect(line['bbox'])
        if spans:
            baseline = float(spans[0]['origin'][1])
            size = float(spans[0]['size'])
            c = int(spans[0].get('color', 0x1A1A1A))
            text_color = ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)
        else:
            baseline, size, text_color = 472.74, 5.2, (0.10, 0.10, 0.10)
        page.add_redact_annot(legend_rect, fill=(1, 1, 1), cross_out=False)
    else:
        baseline = size = None
        text_color = (0.10, 0.10, 0.10)
    for span, _ in pending:
        page.add_redact_annot(fitz.Rect(span['bbox']), fill=(1, 1, 1), cross_out=False)
    if legend_lines or pending:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=0, text=0)
        if legend_lines:
            # Recover the three existing marker centers from the native legend;
            # text is re-emitted at the original baseline with short labels.
            centers = []
            for drawing in page.get_drawings():
                r = drawing['rect']; fill = drawing.get('fill')
                if (fill is not None and 465 < r.y0 < 475 and r.width < 4 and r.height < 4):
                    centers.append((r.x0 + r.x1) / 2)
            centers = sorted(centers)[:3]
            if len(centers) != 3:
                centers = [337.43, 379.61, 421.78]
            for center, label in zip(centers, ('Negative', 'Positive', 'Below threshold')):
                page.insert_text((center + 6.24, baseline), label, fontname='helv',
                                 fontsize=size, color=text_color, overlay=True)
        for span, new_text in pending:
            page.insert_text(
                fitz.Point(span['origin']), new_text, fontname='helv',
                fontsize=float(span['size']), color=(0.10, 0.10, 0.10), overlay=True,
            )
    payload = doc.tobytes(garbage=4, deflate=True)
    doc.close()
    path.write_bytes(payload)
    rendered = fitz.open(stream=payload, filetype='pdf')
    rendered[0].get_pixmap(dpi=300, alpha=False).save(
        PROJECT_ROOT / 'figures_png/main_figures/Fig4_reorganized.png'
    )
    rendered.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--overview-only', action='store_true')
    parser.add_argument('--reorganized', action='store_true')
    parser.add_argument('--compact-candidate', action='store_true')
    parser.add_argument('--removed-supplement', action='store_true')
    parser.add_argument('--display-labels-only', action='store_true')
    args = parser.parse_args()
    if args.display_labels_only:
        display_labels_only()
        return
    if args.compact_candidate:
        import sys
        draw_reorganized(None, sys.stdin.buffer.read(),compact_candidate=True)
        return
    if args.removed_supplement:
        import sys
        draw_removed_supplement(sys.stdin.buffer.read())
        return
    if args.reorganized:
        import sys
        draw_reorganized(None, sys.stdin.buffer.read())
        return
    if args.overview_only:
        with COUNT_PATH.open() as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
    else:
        records = build_counts()
    draw_overview(records)


if __name__ == '__main__':
    main()
