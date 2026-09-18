#!/usr/bin/env python3
"""Figure 3: retained fixed-score distributions and donor-balanced descriptions.

No hypothesis tests, model fitting, source-score mutation, or new cell annotation.
Owner-authorized column z scores are heatmap-only display transformations.
Local inferential frames are read from the retained donor-contrast table.
The selected nine-panel PDF is a fixed input, never a regenerated output.
"""
from pathlib import Path
import os
import argparse
import io
import re
import textwrap

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

R = Path(os.environ["NMF_SOURCE_ROOT"])
P = Path(os.environ["NMF_ARCHIVE_SOURCE_ROOT"])
OUT = Path(os.environ["NMF_WORK_ROOT"]) / "analysis/fig3_layer_region"
SOURCE = Path(__file__).with_name("Fig3_selected_nine_panel.pdf")
ANNOTATION = R / "tables/TableS3_program_annotation.tsv"
IDENTITY = R / "inputs/cortex_nmf_program/archived/revisions/2026-08-02_v49_professor_review/human_validation/HUMAN_VALIDATION_spatial_section_by_layer_aggregates_all54.tsv"
CROSS = P / "results/crossregion_v1"
F2 = P / "scripts/fig2"
LAYERS = ["L1", "L2", "L3", "L4", "L5", "L6"]
DOMAINS = ["ARACHNOID"] + LAYERS + ["WM"]
REGIONS = ["ACC", "AG", "DLPFC", "FPPFC", "ITG", "M1", "PoCG", "S1",
           "S1E", "SMG", "SPL", "STG", "V1", "VLPFC"]
DONORS = ["DonorB", "DonorF", "DonorG", "DonorH", "DonorI"]
DONOR_COLORS = dict(zip(DONORS, ["#4779B5", "#D48B35", "#479B88", "#9B70AB", "#66727C"]))
DONOR_MARKERS = dict(zip(DONORS, ["o", "s", "^", "D", "v"]))
DOMAIN_COLORS = dict(zip(DOMAINS, ["#9E9E9E", "#3B4CC0", "#5A78D6", "#7DA0E0",
                                 "#36A66B", "#E8A93B", "#D65A3B", "#7B3294"]))
SCALE = "fixed_global_program_z_of_SCT_corrected_count_projection"
DIST_PATH = OUT / "layer_bin_distributions.tsv"
PROFILE_PATH = OUT / "layer_region_profiles.tsv"
CONTRAST_PATH = OUT / "fig3_local_donor_contrasts.tsv"


def inputs():
    ann = pd.read_csv(ANNOTATION, sep="\t", dtype=str)
    ann["number"] = ann.new_P.str.removeprefix("P").astype(int)
    ann = ann.sort_values("number").set_index("new_P")
    if len(ann) != 54:
        raise ValueError("The adopted repertoire must contain exactly 54 programs")
    ids = pd.read_csv(IDENTITY, sep="\t", usecols=["chip", "region", "donor"]).drop_duplicates()
    if ids.chip.duplicated().any():
        raise ValueError("The adopted chip-to-donor mapping is not unique")
    return ann, ids


def summarize_bins(ann, ids):
    raw_names = ["program_" + str(x) for x in ann.cnmf_component]
    programs = ann.index.tolist()
    meta = pq.read_table(CROSS / "spatial_bin50_meta.parquet",
                         columns=["bin", "chip", "region", "majorDomain"]).to_pandas()
    meta = meta.rename(columns={"majorDomain": "layer"}).set_index("bin")
    scores = pq.read_table(CROSS / "spatial_bin50_program_score_SCT.parquet",
                           columns=["bin"] + raw_names).to_pandas().set_index("bin")
    # This parquet is already the adopted global per-program z score, not rawdot.
    frame = scores.join(meta, how="left", validate="one_to_one")
    del scores, meta
    if frame.chip.isna().any():
        raise ValueError("A fixed-score bin is missing its adopted spatial metadata")
    frame = frame.loc[frame.layer.isin(LAYERS)]
    donor_for_chip = ids.set_index("chip").donor.to_dict()
    records = []

    def add(block, record_type, chip, region, layer):
        values = block[raw_names].to_numpy(dtype=np.float64)
        total = len(values)
        for j, pid in enumerate(programs):
            v = values[:, j]
            v = v[np.isfinite(v)]
            if len(v):
                lo, q25, median, q75, hi = np.quantile(v, [0, .25, .5, .75, 1],
                                                      method="linear")
            else:
                lo = q25 = median = q75 = hi = np.nan
            a = ann.loc[pid]
            records.append(dict(
                record_type=record_type, program=pid,
                raw_component=a.cnmf_component, functional_name=a.functional_name,
                reference_subclass=a.dominant_subclass, confidence=a.confidence,
                chip=chip, donor=donor_for_chip[chip] if chip else "",
                region=region, layer=layer, n_bins=total, n_finite=len(v),
                n_nonfinite=total-len(v), minimum=lo, q25=q25, median=median,
                q75=q75, maximum=hi, score_scale=SCALE))
        print("BIN_DISTRIBUTION", record_type, chip, region, layer, total, flush=True)

    for (chip, region, layer), block in frame.groupby(["chip", "region", "layer"], sort=True):
        if chip not in donor_for_chip:
            raise ValueError("No true donor identity for section " + chip)
        add(block, "chip_layer_bins", chip, region, layer)
    # Pooled-bin IQRs are display summaries, not biological-replicate intervals.
    for layer in LAYERS:
        add(frame.loc[frame.layer.eq(layer)], "pooled_layer_bins", "", "ALL", layer)
    del frame
    return pd.DataFrame(records)


def regional_profiles(ann, ids):
    raw_to_p = {"program_" + str(row.cnmf_component): pid for pid, row in ann.iterrows()}
    old = pd.read_csv(F2 / "prog_x_layer_per_chip.tsv", sep="\t")
    old = old.loc[old.program.isin(raw_to_p) & old.majorDomain.isin(LAYERS)].copy()
    old["program"] = old.program.map(raw_to_p)
    # Historical mean_z is a chip-layer bin MEDIAN. Do not substitute the
    # separate GEE mean response, or recompute scores/old statistical results.
    old = old.rename(columns={"majorDomain": "layer", "mean_z": "section_median_score"})
    old = old.merge(ids, on="chip", how="left", validate="many_to_one")
    if old[["region", "donor"]].isna().any().any():
        raise ValueError("A retained section-layer median has no true donor/region")
    donor = old.groupby(["program", "region", "layer", "donor"], sort=True).agg(
        score=("section_median_score", "mean"),
        n_sections=("chip", "nunique"), n_bins=("n", "sum")).reset_index()
    lookup = donor.set_index(["program", "region", "layer", "donor"])
    rows = []
    for pid in ann.index:
        a = ann.loc[pid]
        for region in REGIONS:
            for layer in LAYERS:
                obs = donor.loc[donor.program.eq(pid) & donor.region.eq(region) & donor.layer.eq(layer)]
                nd = obs.donor.nunique()
                common = dict(program=pid, raw_component=a.cnmf_component,
                              functional_name=a.functional_name,
                              reference_subclass=a.dominant_subclass,
                              confidence=a.confidence, region=region, layer=layer,
                              n_donors_observed=nd, single_donor_region=nd == 1,
                              score_scale=SCALE,
                              section_statistic="median_of_fixed_bin_z_scores",
                              donor_statistic="equal_section_mean_of_section_medians",
                              region_statistic="equal_observed_donor_mean")
                for person in DONORS:
                    key = (pid, region, layer, person)
                    if key in lookup.index:
                        z = lookup.loc[key]
                        value, ns, nb, present = float(z.score), int(z.n_sections), int(z.n_bins), True
                    else:
                        value, ns, nb, present = np.nan, 0, 0, False
                    rows.append(dict(**common, record_type="donor", donor=person,
                                     observed=present, score=value, n_sections=ns, n_bins=nb))
                rows.append(dict(**common, record_type="region_summary", donor="",
                                 observed=nd > 0, score=float(obs.score.mean()) if nd else np.nan,
                                 n_sections=int(obs.n_sections.sum()), n_bins=int(obs.n_bins.sum())))
    return pd.DataFrame(rows)


def examples(ann, distributions):
    pool = distributions.loc[distributions.record_type.eq("pooled_layer_bins")]
    matrix = pool.pivot(index="program", columns="layer", values="median").reindex(columns=LAYERS)
    # The retained first production showed nearly identical broad P8/P21
    # profiles. Use four non-redundant, higher-confidence functional examples
    # actually observed in the L1-L6 summaries, without forcing one per layer:
    # P5 has a more concentrated L2 pattern, P13 an L4 preference, P28 an L6
    # preference, and P8 is an explicitly broad neuronal comparator.
    selected = ["P5", "P13", "P28", "P8"]
    peaks = {pid: matrix.loc[pid].idxmax() for pid in selected}
    return selected, peaks, matrix


def render(ann, ids, distributions, profiles, local_significance=False):
    """Repaint retained values only; the two L1-L6 tables are never rewritten."""
    import fitz
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize, TwoSlopeNorm, LinearSegmentedColormap, to_rgba
    from matplotlib.ticker import MaxNLocator
    from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

    display_domains = LAYERS + ["WM"]
    case_programs = ["P16"] if local_significance else ["P16", "P40", "P48"]
    c_frames = set()
    if local_significance:
        contrasts = pd.read_csv(CONTRAST_PATH, sep="\t")
        sig = contrasts.loc[np.isfinite(contrasts.q_BY) & contrasts.q_BY.lt(.05)]
        c_frames = set(sig.loc[sig.panel.eq("c"), ["program", "region", "domain"]].itertuples(index=False, name=None))
    spatial_cases = [("VLPFC", "B02111F5", "DonorG"),
                     ("V1", "B02222E1", "DonorG"),
                     ("V1", "D00865B3", "DonorF")]
    subclass_order = ["L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB",
                      "L6 IT", "L6 CT", "L6B", "L6 CAR3", "ET", "NP",
                      "PVALB", "SST", "VIP", "LAMP5", "NDNF", "PAX6",
                      "CHANDELIER", "AST", "MICRO", "OLIGO", "OPC", "ENDO", "VLMC"]
    subclass_short = {"L2-L3 IT LINC00507": "L2/3 IT", "L3-L4 IT RORB": "L3/4 IT",
                      "L4-L5 IT RORB": "L4/5 IT", "L6B": "L6b",
                      "CHANDELIER": "Chandelier", "AST": "Ast", "MICRO": "Micro",
                      "OLIGO": "Oligo", "ENDO": "Endo"}
    old_ct = ["AST", "CHANDELIER", "ENDO", "ET", "L2-L3 IT LINC00507",
              "L3-L4 IT RORB", "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT",
              "L6B", "LAMP5", "MICRO", "NDNF", "NP", "OLIGO", "OPC", "PAX6",
              "PVALB", "SST", "VIP", "VLMC"]
    ct_colors = dict(zip(old_ct, [
        "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
        "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
        "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
        "#17becf", "#9edae5", "#393b79", "#637939"]))
    order = sorted(ann.index, key=lambda p: (subclass_order.index(ann.loc[p, "dominant_subclass"]),
                                           int(p[1:])))
    raw = ["program_" + str(ann.loc[p, "cnmf_component"]) for p in order]
    raw_to_new = dict(zip(raw, order))
    global_data = pd.read_csv(F2 / "prog_x_layer_global.tsv", sep="\t").set_index("majorDomain")
    global_matrix = global_data.loc[display_domains, raw].to_numpy(float)
    if not local_significance:
        correlations = pd.read_csv(F2 / "panelg_reproducibility.tsv", sep="\t")
        correlations = correlations.loc[correlations.program.isin(raw)].copy()
        correlation_summary = pd.read_csv(F2 / "panelg_summary.tsv", sep="\t").set_index("program")
        correlation_summary = correlation_summary.loc[raw]
    region_data = profiles.loc[profiles.record_type.eq("region_summary")]
    donor_data = profiles.loc[profiles.record_type.eq("donor")]
    regional_lookup = region_data.set_index(["program", "region", "layer"])
    case_ids = [chip for _, chip, _ in spatial_cases]
    meta = pq.read_table(CROSS / "spatial_bin50_meta.parquet",
                         columns=["bin", "chip", "x", "y", "majorDomain"],
                         filters=[("chip", "in", case_ids)]).to_pandas()
    wanted = set(meta.bin)
    case_raw = ["program_" + str(ann.loc[p, "cnmf_component"]) for p in case_programs]
    score_parts = []
    for batch in pq.ParquetFile(CROSS / "spatial_bin50_program_score_SCT.parquet").iter_batches(
            columns=["bin"] + case_raw, batch_size=400000):
        block = batch.to_pandas()
        block = block.loc[block.bin.isin(wanted)]
        if len(block):
            score_parts.append(block)
    scores = pd.concat(score_parts, ignore_index=True)
    tissue = meta.merge(scores, on="bin", how="left", validate="one_to_one")
    if tissue[case_raw].isna().any().any():
        raise ValueError("A selected tissue bin has no retained fixed-score value")
    for area, chip, person in spatial_cases:
        identity = ids.loc[ids.chip.eq(chip)]
        if len(identity) != 1 or identity.iloc[0].donor != person or identity.iloc[0].region != area:
            raise ValueError("Selected tissue identity differs from the approved specimen and donor")

    from PIL import Image
    from matplotlib import font_manager
    W, H, MM = 170.0, 200.0, 72/25.4
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.0,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": .5})
    fig = plt.figure(figsize=(W/25.4,H/25.4), facecolor="none")
    font_path = font_manager.findfont("DejaVu Sans")
    pdf_font = fitz.Font(fontfile=font_path)
    native = fitz.open(stream=SOURCE.read_bytes(), filetype="pdf")
    native_spans = []
    for block in native[0].get_text("dict", flags=fitz.TEXTFLAGS_DICT & ~fitz.TEXT_PRESERVE_IMAGES)["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                native_spans.append(dict(span, direction=line["dir"]))
    clean = fitz.open(stream=SOURCE.read_bytes(), filetype="pdf")
    for span in native_spans:
        rect = fitz.Rect(span["bbox"]) & clean[0].rect
        if not rect.is_empty:
            clean[0].add_redact_annot(rect, fill=False, cross_out=False)
    clean[0].apply_redactions(images=0, graphics=0)

    def axmm(x,y,w,h):
        return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],facecolor="none")
    def textmm(x,y,value,**kw):
        options=dict(ha="left",va="top",fontsize=7.0,color="#20252B")
        options.update(kw)
        return fig.text(x/W,1-y/H,value,**options)
    def patchmm(x,y,w,h,color):
        fig.patches.append(Rectangle((x/W,1-(y+h)/H),w/W,h/H,
            transform=fig.transFigure,facecolor=color,edgecolor="none"))
    def plot_pdf():
        buffer=io.BytesIO()
        fig.savefig(buffer,format="pdf",dpi=300,transparent=True)
        plt.close(fig)
        return fitz.open(stream=buffer.getvalue(),filetype="pdf")
    def pdf_text(page,x,y,value,size=7.0,rotate=0):
        page.insert_font(fontname="figurefont",fontfile=font_path)
        page.insert_text((x*MM,y*MM),value,fontname="figurefont",
                         fontsize=max(7.0,size),rotate=rotate,color=(.12,.15,.18))
    def pdf_box(page,box,value,size=7.3):
        page.insert_font(fontname="figurefont",fontfile=font_path)
        page.insert_textbox(fitz.Rect(*(v*MM for v in box)),value,
                           fontname="figurefont",fontsize=size,lineheight=1.18,
                           color=(.12,.15,.18))
    def put_native(page,clip,x,y,width,with_text=True,omit=None):
        clip=fitz.Rect(clip)
        scale=width*MM/clip.width
        height=width*clip.height/clip.width
        destination = fitz.Rect(x*MM,y*MM,(x+width)*MM,(y+height)*MM)
        if local_significance:
            # Remove outside content before import, not just a display clip.
            isolated = fitz.open(stream=clean.tobytes(), filetype="pdf")
            ip = isolated[0]
            bounds = ip.rect
            outside = [
                fitz.Rect(bounds.x0,bounds.y0,bounds.x1,clip.y0),
                fitz.Rect(bounds.x0,clip.y1,bounds.x1,bounds.y1),
                fitz.Rect(bounds.x0,clip.y0,clip.x0,clip.y1),
                fitz.Rect(clip.x1,clip.y0,bounds.x1,clip.y1)]
            for rect in outside:
                if not rect.is_empty:
                    ip.add_redact_annot(rect,fill=False,cross_out=False)
            ip.apply_redactions(images=2,graphics=1,text=0)
            ip.clean_contents(sanitize=True)
            ip.set_cropbox(clip)
            compact = fitz.open(stream=isolated.tobytes(garbage=4,deflate=True),filetype="pdf")
            page.show_pdf_page(destination,compact,0,keep_proportion=True)
            compact.close();isolated.close()
        else:
            page.show_pdf_page(destination,clean,0,clip=clip,keep_proportion=True)
        if with_text:
            page.insert_font(fontname="figurefont",fontfile=font_path)
            for sp in native_spans:
                rect=fitz.Rect(sp["bbox"])
                center=(rect.tl+rect.br)/2
                value=sp["text"].strip()
                if not clip.contains(center) or not value:
                    continue
                if value in list("abcdefghi") and sp["size"]>=7:
                    continue
                if value == "median program z" or (omit is not None and omit(sp)):
                    continue
                ox,oy=sp["origin"]
                dx,dy=sp["direction"]
                rotation=90 if dy<-.5 else 270 if dy>.5 else 180 if dx<-.5 else 0
                if rotation==0 and re.search(r"\d\s{2,}",sp["text"]):
                    for word in native[0].get_text("words"):
                        wr=fitz.Rect(word[:4]);wc=(wr.tl+wr.br)/2
                        if rect.contains(wc):
                            tx=x*MM+(wc.x-clip.x0)*scale
                            page.insert_text((tx-pdf_font.text_length(word[4],fontsize=7)/2,
                                              y*MM+(oy-clip.y0)*scale),
                                             word[4],fontname="figurefont",fontsize=7,
                                             color=(.08,.08,.08))
                else:
                    page.insert_text((x*MM+(ox-clip.x0)*scale,y*MM+(oy-clip.y0)*scale),
                                     sp["text"],fontname="figurefont",fontsize=7,
                                     rotate=rotation,color=(.08,.08,.08))
        return height
    def matrix(ax,values,cmap,norm):
        nr,nc=values.shape
        ax.pcolormesh(np.arange(nc+1)-.5,np.arange(nr+1)-.5,values,
                      cmap=cmap,norm=norm,shading="flat",edgecolors="none")
        ax.set_xlim(-.5,nc-.5);ax.set_ylim(nr-.5,-.5)
        ax.tick_params(length=0,pad=1)
        for spine in ax.spines.values():spine.set_visible(False)
    score_cmap=LinearSegmentedColormap.from_list("fixed_global_z",["#2868A0","#F7F7F7","#C45F4D"])
    score_norm=Normalize(-2,2,clip=True)
    heatmap_cmap=score_cmap.copy();heatmap_cmap.set_bad("#E1E4E8")
    def scaled_columns(values):
        table=pd.DataFrame(np.where(np.isfinite(values),values,np.nan))
        means=table.mean(axis=0,skipna=True)
        sds=table.std(axis=0,skipna=True,ddof=1)
        sds=sds.where(sds.gt(0)&np.isfinite(sds))
        shown=table.sub(means,axis=1).div(sds,axis=1).to_numpy(float)
        finite=shown[np.isfinite(shown)]
        span=np.ceil(np.max(np.abs(finite))*10)/10 if len(finite) else 1.0
        return shown,TwoSlopeNorm(vmin=-span,vcenter=0,vmax=span),span
    def physical_axes(ax,frame):
        xx=frame.x.to_numpy(float)/2;yy=frame.y.to_numpy(float)/2
        lim=(xx.min()-12.5,xx.max()+12.5,yy.min()-12.5,yy.max()+12.5)
        ax.set_xlim(lim[0],lim[1]);ax.set_ylim(lim[3],lim[2])
        ax.set_aspect("equal",adjustable="box");ax.set_axis_off()
        return xx,yy,lim
    def scale_bar(ax,lim):
        xmin,xmax,ymin,ymax=lim
        bx=xmin+.05*(xmax-xmin);by=ymax-.025*(ymax-ymin)
        ax.plot([bx,bx+1000],[by,by],color="#20252B",lw=.8,solid_capstyle="butt")
        ax.text(bx,by-.018*(ymax-ymin),"1 mm",ha="left",va="bottom",fontsize=7,
                bbox=dict(facecolor="white",alpha=.72,edgecolor="none",pad=.1))
    def domain_labels(ax,frame,xx,yy,lim):
        xmin,xmax,ymin,ymax=lim
        def anchor(domain,fx,fy):
            mask=frame.majorDomain.eq(domain).to_numpy()
            xs,ys=xx[mask],yy[mask]
            if not len(xs):return None
            j=np.argmin(((xs-xmin)/(xmax-xmin)-fx)**2+((ys-ymin)/(ymax-ymin)-fy)**2)
            return xs[j],ys[j]
        for j,domain in enumerate(LAYERS):
            point=anchor(domain,.82-j*.065,.18+j*.075)
            if point is not None:
                ax.annotate(domain,xy=point,xycoords="data",xytext=(1.025,.96-j*.137),
                    textcoords="axes fraction",ha="left",va="center",fontsize=7,
                    color=DOMAIN_COLORS[domain],fontweight="bold",annotation_clip=False,
                    arrowprops=dict(arrowstyle="-",color="#535B61",lw=.3,shrinkA=1,shrinkB=0))
        point=anchor("WM",.45,.59)
        if point is not None:
            ax.annotate("WM",xy=point,xycoords="data",xytext=(.74,-.065),
                textcoords="axes fraction",ha="center",va="center",fontsize=7,
                color=DOMAIN_COLORS["WM"],fontweight="bold",annotation_clip=False,
                arrowprops=dict(arrowstyle="-",color="#535B61",lw=.3,shrinkA=1,shrinkB=0))
    def domain_map(ax,frame):
        xx,yy,lim=physical_axes(ax,frame)
        colours=frame.majorDomain.map(DOMAIN_COLORS).fillna("#ECECEC").to_numpy()
        colours[frame.majorDomain.eq("ARACHNOID").to_numpy()]="#ECECEC"
        ax.scatter(xx,yy,c=colours,s=.14,marker="s",linewidths=0,rasterized=True)
        domain_labels(ax,frame,xx,yy,lim);scale_bar(ax,lim)
    seen_donors=["DonorF","DonorG","DonorH"]
    y_ranges={}
    for pid in case_programs:
        vals=donor_data.loc[donor_data.program.eq(pid)&donor_data.region.isin(["VLPFC","V1"])&
                           donor_data.observed,"score"].to_numpy(float)
        low=np.floor(np.nanmin(vals)*2)/2;high=np.ceil(np.nanmax(vals)*2)/2
        y_ranges[pid]=(low-.08*(high-low),high+.08*(high-low))
    def case_headers(y):
        for ci,(region,chip,person) in enumerate(spatial_cases):
            textmm(12.5+28*ci,y,region+" · "+person.replace("Donor","")+"\n"+chip,
                   ha="center",fontsize=7,linespacing=1.0)
        for ri,region in enumerate(["VLPFC","V1"]):
            textmm(108.5+43*ri,y,region+" (3 donors)\nProfile: global z",
                   ha="center",fontsize=7,linespacing=1.0)
    def case_row(pid,y,height=20):
        pcol="program_"+str(ann.loc[pid,"cnmf_component"])
        for ci,(region,chip,person) in enumerate(spatial_cases):
            frame=tissue.loc[tissue.chip.eq(chip)]
            ax=axmm(1+28*ci,y,23,height)
            xx,yy,lim=physical_axes(ax,frame)
            cortex=frame.majorDomain.isin(LAYERS).to_numpy()
            ax.scatter(xx[~cortex],yy[~cortex],color="#E6E8EB",s=.14,marker="s",
                       linewidths=0,rasterized=True)
            ax.scatter(xx[cortex],yy[cortex],c=frame.loc[cortex,pcol].to_numpy(float),
                       cmap=score_cmap,norm=score_norm,s=.14,marker="s",
                       linewidths=0,rasterized=True)
            domain_labels(ax,frame,xx,yy,lim);scale_bar(ax,lim)
        for ri,region in enumerate(["VLPFC","V1"]):
            ax=axmm(93+43*ri,y+.5,32,height-1.2)
            for person in seen_donors:
                part=donor_data.loc[donor_data.program.eq(pid)&donor_data.region.eq(region)&
                                   donor_data.donor.eq(person)].set_index("layer").reindex(LAYERS)
                values=part.score.to_numpy(float)
                values[~part.observed.fillna(False).to_numpy(bool)]=np.nan
                ax.plot(range(6),values,color=DONOR_COLORS[person],marker=DONOR_MARKERS[person],
                        markersize=2,linewidth=.9 if person=="DonorG" else .6,alpha=.9)
            values=[regional_lookup.loc[(pid,region,layer),"score"] for layer in LAYERS]
            ax.plot(range(6),values,color="#20252B",lw=.8,linestyle=(0,(3,1.5)))
            ax.set_xlim(-.2,5.2);ax.set_ylim(*y_ranges[pid])
            ax.set_xticks(range(6),LAYERS)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3,steps=[1,2,5,10]))
            ax.tick_params(labelsize=7,length=1.2,pad=.5)
            for side in ["top","right"]:ax.spines[side].set_visible(False)
    def donor_key(y,all_five=False,iqr=False):
        people=DONORS if all_five else seen_donors
        handles=[Line2D([],[],color=DONOR_COLORS[d],marker=DONOR_MARKERS[d],
                        markersize=3,lw=0 if iqr else .7,
                        label=d.replace("Donor","Donor ")) for d in people]
        if iqr:
            handles += [Rectangle((0,0),1,1,facecolor="#EDF0F2",edgecolor="#8A929A",lw=.5,label="Bin IQR"),
                        Line2D([],[],color="#20252B",lw=.8,label="Median")]
        else:
            handles += [Line2D([],[],color="#20252B",lw=.8,linestyle=(0,(3,1.5)),label="Region mean")]
        return fig.legend(handles=handles,loc="upper center",
                   bbox_to_anchor=((W/2 if all_five else 108)/W,1-y/H),
                   ncol=len(handles),frameon=False,fontsize=7,
                   handletextpad=.3,columnspacing=.6,borderpad=0,borderaxespad=0)
    def spatial_key(y):
        textmm(1.5,y,"Bin global z",fontsize=7)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=score_norm,cmap=score_cmap),
             cax=axmm(1.5,y+3,28,1.1),orientation="horizontal",ticks=[-2,0,2],extend="both")
        cb.ax.tick_params(labelsize=7,length=.8,pad=.3);cb.outline.set_visible(False)

    # Main a: retained native spatial evidence, relabelled at final-size 7 pt.
    textmm(.4,.1,"a",fontsize=9,fontweight="bold")
    textmm(8,.4,"VLPFC · Donor G · B02111F5",fontsize=7)
    a_labels=[("Excitatory","Activity-dep. IEG\nP31"),
              ("Inh (VIP)†","Synaptic-vesicle\nexocytosis (P21)"),
              ("Astrocyte","Glutamate transp.\nP14"),
              ("Oligodendrocyte","Oligo/myelin\nP33"),
              ("Endothelial","Vessel morphog.\nP51"),
              ("Microglia","Immune activation\nP36")]
    for i,(cell,label) in enumerate(a_labels):
        cx=(47.73+81.87*i)*142/503.643
        pid=["P31","P21","P14","P33","P51","P36"][i]
        col=ct_colors[ann.loc[pid,"dominant_subclass"]]
        textmm(cx,4.0,cell,ha="center",fontsize=7,color=tuple(.55*v for v in to_rgba(col)[:3]))
        textmm(cx,7.0,label,ha="center",fontsize=7,linespacing=.98)
    rep=pd.read_csv(F2/"repchip_meta.tsv",sep="\t")
    textmm(154.5,6,"VLPFC",ha="center",fontsize=7)
    domain_map(axmm(144.0,13.3,19.0,23.0),rep)
    textmm(107,39.2,"Activity (per program)",ha="right",fontsize=7)
    for x,value in [(112,"low"),(125.5,"0"),(139,"high")]:
        textmm(x,36.6,value,ha="center",fontsize=7)
    # Main b: all 54 programs retained on a single annotated column axis.
    textmm(.4,43.5,"b",fontsize=9,fontweight="bold")
    x0,pitch,pw=19.0,2.62,54*2.62
    groups=[]
    for subclass in subclass_order:
        positions=[j for j,p in enumerate(order) if ann.loc[p,"dominant_subclass"]==subclass]
        if positions:groups.append((subclass,min(positions),max(positions)+1))
    for subclass,begin,end in groups:
        cx=x0+pitch*(begin+end)/2
        textmm(cx,54.7,subclass_short.get(subclass,subclass),fontsize=7,rotation=65,
               ha="left",va="bottom",rotation_mode="anchor",
               color=tuple(.55*v for v in to_rgba(ct_colors[subclass])[:3]))
        patchmm(x0+pitch*begin,56.0,pitch*(end-begin)-.2,1.3,ct_colors[subclass])
    textmm(.5,50,"Program ref.\npreference",fontsize=7,linespacing=.98)
    gz,gn,gs=scaled_columns(global_matrix)
    ga=axmm(x0,58.0,pw,12.6);matrix(ga,gz,heatmap_cmap,gn)
    ga.set_xticks([]);ga.set_yticks([])
    for i,domain in enumerate(display_domains):
        tx=17.8 if i%2==0 else 11.8
        ty=58.0+(i+.5)*12.6/7
        textmm(tx,ty,domain,fontsize=7,ha="right",va="center")
        fig.lines.append(Line2D([(tx+.35)/W,18.7/W],[1-ty/H,1-ty/H],
                               transform=fig.transFigure,color="#70777D",lw=.3))
    textmm(4.0,64.3,"Column z-score\n7 domains",rotation=90,ha="center",va="center",fontsize=7)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=gn,cmap=heatmap_cmap),
                   cax=axmm(168,58,1,12.6),ticks=[-gs,0,gs])
    cb.ax.yaxis.set_ticks_position("left");cb.ax.tick_params(labelsize=7,length=.8,pad=1)
    cb.outline.set_visible(False)
    # Owner-selected descriptive display: no cell outlines in panel b.
    for j,pid in enumerate(order):
        suffix="*" if str(ann.loc[pid,"confidence"]).startswith("Lower") else ""
        textmm(x0+pitch*(j+.5),72.8,pid+suffix,fontsize=7,rotation=90,
               ha="right",va="center",rotation_mode="anchor")

    # Independent regional panel: all 84 area-layer rows, column z over the
    # complete matrix, then independent Euclidean/complete-linkage ordering.
    # Colour saturation is display-only and never enters either linkage.
    selected_layers=LAYERS
    row_keys=[(r,l) for r in REGIONS for l in selected_layers]
    rv=np.array([[regional_lookup.loc[(p,r,l),"score"] for p in order]
                 for r,l in row_keys],dtype=float)
    rz,unused_regional_norm,raw_display_span=scaled_columns(rv)
    if not np.isfinite(rz).all():
        raise ValueError("Regional clustering requires finite 84-by-54 displayed z values; no missing values were imputed")
    row_link=linkage(rz,method="complete",metric="euclidean",optimal_ordering=False)
    col_link=linkage(rz.T,method="complete",metric="euclidean",optimal_ordering=False)
    row_order=leaves_list(row_link);col_order=leaves_list(col_link)
    clustered=rz[row_order][:,col_order]
    rn=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=3)
    textmm(.4,83.0,"c",fontsize=9,fontweight="bold")
    textmm(77,83.3,"Column z-score",fontsize=7)
    textmm(120,83.3,"* Low annot. conf.",fontsize=7)
    ca=axmm(x0,87, pw,4.2)
    dendrogram(col_link,ax=ca,no_labels=True,color_threshold=0,
               above_threshold_color="#4E575E",link_color_func=lambda k:"#4E575E")
    ca.set_xlim(0,54*10);ca.set_axis_off()
    for collection in ca.collections:collection.set_linewidth(.35)
    for j,source_j in enumerate(col_order):
        pid=order[int(source_j)]
        patchmm(x0+j*pitch,92.0,pitch-.05,1.2,ct_colors[ann.loc[pid,"dominant_subclass"]])
    textmm(10.8,91.8,"Ref.",fontsize=7)
    area_colors=dict(zip(REGIONS,[plt.get_cmap("tab20").colors[i]
                          for i in [0,2,4,6,8,10,12,14,16,18,1,3,5,7]]))
    regional_y,regional_h=94.5,39.2
    rd=axmm(.5,regional_y,8.5,regional_h)
    dendrogram(row_link,ax=rd,orientation="left",no_labels=True,color_threshold=0,
               above_threshold_color="#4E575E",link_color_func=lambda k:"#4E575E")
    rd.set_ylim(84*10,0);rd.set_axis_off()
    for collection in rd.collections:collection.set_linewidth(.35)
    identities=[row_keys[int(i)] for i in row_order]
    for sx,width,colors in [
        (11.0,2.0,[to_rgba(area_colors[r]) for r,l in identities]),
        (14.0,1.3,[to_rgba(DOMAIN_COLORS[l]) for r,l in identities])]:
        strip=axmm(sx,regional_y,width,regional_h)
        strip.imshow(np.array(colors)[:,None,:],aspect="auto",interpolation="none")
        strip.set_axis_off()
    # Annotation types are identified by their matching full legends below.
    ra=axmm(x0,regional_y,pw,regional_h);matrix(ra,clustered,heatmap_cmap,rn)
    for i,(region,layer) in enumerate(identities):
        for j,source_j in enumerate(col_order):
            if (order[int(source_j)],region,layer) in c_frames:
                ra.add_patch(Rectangle((j-.48,i-.48),.96,.96,facecolor="none",
                                       edgecolor="black",lw=.28,zorder=4))
    ra.set_xticks([]);ra.set_yticks([])
    for j,source_j in enumerate(col_order):
        pid=order[int(source_j)]
        suffix="*" if str(ann.loc[pid,"confidence"]).startswith("Lower") else ""
        textmm(x0+pitch*(j+.5),134.6,pid+suffix,fontsize=7,rotation=90,
               ha="right",va="center",rotation_mode="anchor")
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=rn,cmap=heatmap_cmap),
                   cax=axmm(168,97.0,1,34.2),ticks=[-3,0,3],extend="both")
    cb.ax.yaxis.set_ticks_position("left");cb.ax.tick_params(labelsize=7,length=.8,pad=1)
    cb.outline.set_visible(False)
    counts={r:int(regional_lookup.loc[("P16",r,"L2"),"n_donors_observed"]) for r in REGIONS}
    textmm(.7,142.7,"Area (donors)",fontsize=7)
    textmm(112,142.7,"Layer",fontsize=7)
    for i,region in enumerate(REGIONS):
        sx=1.0+(i%7)*15.0;sy=145.8+(i//7)*3.2
        patchmm(sx,sy+.2,1.3,1.3,area_colors[region])
        textmm(sx+1.7,sy,region+" ("+str(counts[region])+")",fontsize=7)
    for i,layer in enumerate(LAYERS):
        sx=112+(i%3)*18;sy=145.8+(i//3)*3.2
        patchmm(sx,sy+.2,1.3,1.3,DOMAIN_COLORS[layer])
        textmm(sx+2,sy,layer,fontsize=7)

    if local_significance:
        textmm(.7,153.3,"1-donor areas: no regional inference; effect signs and unavailable contrasts: companion table.",fontsize=7)

    # The unchanged correlations are rendered only in the S10 continuation.
    if not local_significance:
        correlations["donor"]=correlations.chip.map(ids.set_index("chip").donor)
        if correlations.donor.isna().any():raise ValueError("Unmapped original correlation donor")
    def correlation_panel(y,height=11.5):
        ax=axmm(x0,y,pw,height)
        for j,name in enumerate(raw):
            block=correlations.loc[correlations.program.eq(name)]
            for di,person in enumerate(DONORS):
                part=block.loc[block.donor.eq(person)].sort_values("chip")
                if not len(part):continue
                jitter=np.linspace(-.04,.04,len(part)) if len(part)>1 else np.array([0.0])
                ax.scatter(j+(di-2)*.115+jitter,part.corr_to_mean,s=2,
                           marker=DONOR_MARKERS[person],color=DONOR_COLORS[person],
                           alpha=.72,edgecolors="none")
            median=float(correlation_summary.loc[name,"median"])
            ax.plot([j-.34,j+.34],[median,median],color="#29323A",lw=.8)
        ax.set_xlim(-.5,53.5);ax.set_ylim(-1.05,1.05)
        ax.set_xticks([]);ax.set_yticks([-1,0,1])
        ax.tick_params(labelsize=7,length=1,pad=.5)
        ax.axhline(0,color="#DADEE1",lw=.3)
        for side in ["top","right","bottom"]:ax.spines[side].set_visible(False)
        for j,pid in enumerate(order):
            suffix="*" if str(ann.loc[pid,"confidence"]).startswith("Lower") else ""
            textmm(x0+pitch*(j+.5),y+height+.8,pid+suffix,fontsize=7,rotation=90,
                   ha="right",va="center",rotation_mode="anchor")
    # Main d retains the principal P16 example unchanged.
    textmm(.4,158.8,"d",fontsize=9,fontweight="bold")
    textmm(6,159.2,"P16 Axon guidance / synaptic adhesion",fontsize=7.5,fontweight="bold")
    textmm(105,159.3,"Ref. L3/4 IT",fontsize=7,color="#486E43")
    case_headers(162.8);case_row("P16",168.9,19.1)
    spatial_key(191.2)
    main_handles=[Line2D([],[],color=DONOR_COLORS[d],marker=DONOR_MARKERS[d],
                         markersize=3,lw=.7,label=d.replace("Donor","")) for d in DONORS]
    main_handles.append(Line2D([],[],color="#20252B",lw=.8,linestyle=(0,(3,1.5)),
                               label="Region mean"))
    fig.legend(handles=main_handles,loc="upper center",bbox_to_anchor=(111/W,1-191.0/H),
               title="Donor",title_fontsize=7,ncol=6,frameon=False,fontsize=7,
               handletextpad=.3,columnspacing=.7,borderpad=0,borderaxespad=0)
    main_overlay=plot_pdf()
    main=fitz.open();page=main.new_page(width=170*MM,height=200*MM)
    put_native(page,(0,43,503.643,124),0,13.3,142,with_text=False)
    put_native(page,(401.999,131.189,459.273,137.99),112,39.5,27,with_text=False)
    page.show_pdf_page(page.rect,main_overlay,0,overlay=True)
    main.save(OUT/"Fig3.pdf",garbage=4,deflate=True)
    page.get_pixmap(dpi=300,alpha=False).save(OUT/"Fig3.png")
    print("RETAINED",OUT/"Fig3.pdf",OUT/"Fig3.png","mm",170,200,
          "main_caption_allowance_mm",25,"font_min_pt",7,flush=True)
    if local_significance:
        print("DISPLAYED_CELL_FRAMES","b",0,"c",len(c_frames),
              "panel b outlines removed by owner; statistics retained unchanged",flush=True)
        print("PRESERVED","a and d; seven-domain values; 84-row values; column-z; complete-linkage ordering; native vectors",flush=True)
        return

    # FigS10 page 1 is rebuilt from fixed native c/d/e/i, never from the output.
    supplement=fitz.open();first=supplement.new_page(width=170*MM,height=225*MM)
    pdf_text(first,2,5,"Figure S10 · Supplementary material",8.5)
    for letter,x,y in [("a",1,11),("b",87,11),("c",1,101),("d",87,101)]:
        pdf_text(first,x,y,letter,9)
    put_native(first,(121,138,253,267),2,13,80)
    put_native(first,(250,138,367,267),88,13,79)
    # Reflow the three native facets independently so the original r labels
    # and all axis labels remain legible at 7 pt without overlapping neighbours.
    facet_info=[
        ((368,166,414,216),2,"P33 Oligo","Pearson r = 0.73","Spearman rs = 0.66",False),
        ((413,166,465.8,216),30,"P51 Vessel","Pearson r = 0.21","Spearman rs = 0.21",True),
        ((461,166,503.643,216),58,"P5 Lipid","Pearson r = 0.33","Spearman rs = 0.41",False)]
    for clip,x,title,pearson,spearman,middle in facet_info:
        pdf_text(first,x+.5,106,title,7.5)
        pdf_text(first,x+.5,110,pearson,7)
        pdf_text(first,x+.5,113,spearman,7)
        put_native(first,clip,x,116,26,omit=lambda sp,m=middle:
                   "Pearson" in sp["text"] or "Spearman" in sp["text"] or
                   (m and sp["bbox"][0]>460 and sp["bbox"][1]<201))
    put_native(first,(367,219,390.7,267),1.0,152,10.0)
    for i,name in enumerate(old_ct):
        x=14.0+(i//11)*36.0;y=154.0+(i%11)*2.8
        colour=tuple(to_rgba(ct_colors[name])[:3])
        first.draw_rect(fitz.Rect(x*MM,(y-2)*MM,(x+1.5)*MM,(y-.5)*MM),
                        color=None,fill=colour,overlay=True)
        pdf_text(first,x+2,y,name,7)
    put_native(first,(330,491,497,573),90,103,76)
    for i,name in enumerate(DOMAINS):
        x=89+(i%4)*20;y=144+(i//4)*3.5
        colour=tuple(to_rgba(DOMAIN_COLORS[name])[:3])
        first.draw_rect(fitz.Rect(x*MM,(y-2.0)*MM,(x+1.5)*MM,(y-.5)*MM),
                        color=None,fill=colour,overlay=True)
        pdf_text(first,x+2,y,name,7)
    for left,clip,title in [
        (88,(334,606,387,657),"P33\nOligodendrocyte/\nmyelin"),
        (115,(385,606,438,657),"P31\nActivity-dependent\nIEG"),
        (142,(435.5,606,490,657),"P7\nCation channel\n(interneuron)")]:
        pdf_box(first,(left,149,left+25.5,159),title,7)
        put_native(first,clip,left,158,25.5,with_text=False)
    # Existing feature colour bars are retained separately; only keys are reset.
    for left,clip in [(88,(347,657,386,661)),(115,(397.5,657,436.5,661)),
                       (142,(448,657,487,661))]:
        put_native(first,clip,left+5,187.0,19.5,with_text=False)
        pdf_text(first,left,185.7,"Median z",7)
        pdf_text(first,left+4.5,193,"low",7)
        pdf_text(first,left+13.5,193,"0",7)
        pdf_text(first,left+20.8,193,"high",7)
    caption1=("Figure S10. Anatomical and cellular contexts of spatial program expression. "
        "a, Retained eight-domain profiles of P5, P48, P33, P51 and P26. "
        "b, RCTD subclass mixtures across the representative tissue. "
        "c, Retained score–weight relationships for P33–OLIGO, P51–ENDO and P5–upper-layer IT; "
        "the original Pearson and Spearman correlations are shown. "
        "d, The original expression-space embedding coloured by domain and by P33, P31 and P7 scores. "
        "All a–d scientific marks are reused from the fixed original figure. "
        "Additional distributions and astrocyte-preferred program examples continue on the next page.")
    pdf_box(first,(2,195,168,224),caption1,7.4)

    # FigS10 page 2 contains e (four distributions), f and g (two AST rows).
    H=225.0;fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor="none")
    textmm(2,1,"Figure S10 · continued",fontsize=8.5)
    textmm(.5,8,"e",fontsize=9,fontweight="bold")
    pool=distributions.loc[distributions.record_type.eq("pooled_layer_bins")]
    sections=distributions.loc[distributions.record_type.eq("chip_layer_bins")]
    restored=["P5","P13","P28","P8"]
    for pid,cx in zip(restored,[7,48,89,130]):
        title=ann.loc[pid,"functional_name"].strip()+" ("+pid+")"
        textmm(cx+17,10,"\n".join(textwrap.wrap(title,24)),fontsize=7,ha="center",
               fontweight="bold",linespacing=1.0)
        subclass=ann.loc[pid,"dominant_subclass"]
        textmm(cx+17,22.5,"Ref. "+subclass_short.get(subclass,subclass),fontsize=7,ha="center",
               color=tuple(.55*v for v in to_rgba(ct_colors[subclass])[:3]))
        ax=axmm(cx,27,34,23)
        sub=pool.loc[pool.program.eq(pid)].set_index("layer").reindex(LAYERS)
        source=sections.loc[sections.program.eq(pid)]
        for j,layer in enumerate(LAYERS):
            v=sub.loc[layer]
            ax.add_patch(Rectangle((j-.24,v.q25),.48,v.q75-v.q25,
                         facecolor=to_rgba(DOMAIN_COLORS[layer],.22),
                         edgecolor=DOMAIN_COLORS[layer],lw=.5,zorder=1))
            ax.plot([j-.24,j+.24],[v["median"]]*2,color="#20252B",lw=.8,zorder=3)
            for di,person in enumerate(DONORS):
                part=source.loc[source.layer.eq(layer)&source.donor.eq(person)].sort_values("chip")
                if not len(part):continue
                jitter=np.linspace(-.038,.038,len(part)) if len(part)>1 else np.array([0.0])
                ax.scatter(j+(di-2)*.078+jitter,part["median"],s=5,
                           marker=DONOR_MARKERS[person],color=DONOR_COLORS[person],
                           alpha=.72,edgecolors="none",zorder=4)
        ax.set_xlim(-.6,5.6);ax.set_xticks(range(6),LAYERS)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3,steps=[1,2,5,10]))
        ax.tick_params(labelsize=7,length=1.2,pad=.5)
        if pid=="P5":ax.set_ylabel("GEP score (global z)",fontsize=7,labelpad=1)
        for side in ["top","right"]:ax.spines[side].set_visible(False)
    donor_key(56,all_five=True,iqr=True)
    case_headers(68)
    for pid,letter,title_y,map_y in [("P40","f",78.0,84.0),("P48","g",118.0,124.0)]:
        textmm(.5,title_y,letter,fontsize=9,fontweight="bold")
        textmm(6,title_y+.3,pid+" "+ann.loc[pid,"functional_name"].strip(),fontsize=7.5,fontweight="bold")
        textmm(109,title_y+.3,"Ref. Ast",fontsize=7,color="#145A7C")
        case_row(pid,map_y,23)
    spatial_key(151);donor_key(154.5)
    textmm(.5,160,"h",fontsize=9,fontweight="bold")
    textmm(7,160.2,"8-domain Pearson r · 34 sections",fontsize=7)
    textmm(112,160.2,"* Low annot. conf.",fontsize=7)
    correlation_panel(164,11.5)
    second_plot=plot_pdf()
    supplement.insert_pdf(second_plot)
    second=supplement[1]
    caption2=("Figure S10 continued. e, Original pooled-bin interquartile ranges and medians across "
        "L1–L6 for P5, P13, P28 and P8; points are all retained section medians, coloured and shaped "
        "by the five true donors. f,g, P40 L-glutamate transport and P48 astrocyte metal homeostasis, "
        "both astrocyte-preferred in the reference. Each row shows VLPFC donor G (B02111F5), "
        "V1 donor G (B02222E1) and V1 donor F (D00865B3), with unchanged scores and a common colour "
        "scale. Grey denotes non-cortical-layer context; leaders indicate the frozen domains. "
        "All scale bars are 1 mm. Right: complete observed F/G/H donor profiles in VLPFC and V1; "
        "the dashed curve is the retained equal-donor regional mean. These are descriptive "
        "functional-profile examples, not evidence of significant regional interactions or strict "
        "layer specificity. h, Original eight-domain Pearson r for all 54 programs in 34 complete "
        "sections, against their inclusive mean; points use the donor key in e and black bars "
        "are the retained medians. Asterisks indicate lower annotation confidence, not significance. "
        "No scores or correlations were recomputed.")
    pdf_box(second,(2,190,168,223),caption2,7.4)
    supplement.save(OUT/"FigS10.pdf",garbage=4,deflate=True)
    previews=[]
    for page in supplement:
        pix=page.get_pixmap(dpi=300,alpha=False)
        previews.append(Image.frombytes("RGB",(pix.width,pix.height),pix.samples))
    gutter=12
    preview=Image.new("RGB",(max(p.width for p in previews),sum(p.height for p in previews)+gutter),"white")
    offset=0
    for im in previews:
        preview.paste(im,(0,offset));offset+=im.height+gutter
    preview.save(OUT/"FigS10.png",dpi=(300,300))
    print("RETAINED",OUT/"FigS10.pdf",OUT/"FigS10.png",
          "pdf_pages",len(supplement),"page_mm",[(170,225),(170,225)],
          "png_role","two-page complete preview","font_min_pt",7,flush=True)
    print("PANEL_MAPPING","Fig3 a unchanged spatial/context; b 7-domain maxima and P labels; c independently clustered 84-area-layer by 54-program matrix; d unchanged P16. FigS10 a-g retained; h original 54-program by 34-section correlation panel",flush=True)
    print("HEATMAP_DISPLAY","independent column z, ddof=1; global 7 rows; regional 84 rows; data not clipped; fixed regional colour range with extensions",
          "global_limits",(-gs,gs),"regional_limits",(-3,3),"unclipped_regional_extent",raw_display_span,flush=True)
    print("CLUSTERING","84 rows and 54 columns independently; untruncated column-z; Euclidean; complete linkage; no optimal ordering",flush=True)
    print("GLOBAL_FRAMES","local contrast frames only; no column-maximum rule",flush=True)
    print("REGIONAL_ROW_ORDER",[row_keys[int(i)] for i in row_order],flush=True)
    print("REGIONAL_COLUMN_ORDER",[order[int(i)] for i in col_order],flush=True)
    print("SPATIAL_CASES",spatial_cases,"scores/curves unchanged; no sampling or smoothing; DNB/2 and 1mm bars",flush=True)
    print("UNCHANGED",DIST_PATH,PROFILE_PATH,"S11-S20; no S21 or fixed-source copy",flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--local-significance", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    ann, ids = inputs()
    if args.local_significance:
        profiles = pd.read_csv(PROFILE_PATH, sep="\t", keep_default_na=True)
        render(ann, ids, None, profiles, local_significance=True)
        return
    if args.render_only:
        distributions = pd.read_csv(DIST_PATH, sep="\t", keep_default_na=True)
        profiles = pd.read_csv(PROFILE_PATH, sep="\t", keep_default_na=True)
    else:
        distributions = summarize_bins(ann, ids)
        profiles = regional_profiles(ann, ids)
        selected, peaks, _ = examples(ann, distributions)
        distributions["display_example"] = distributions.program.isin(selected)
        profiles["display_example"] = profiles.program.isin(selected)
        profiles["illustrated_layer"] = profiles.program.map(peaks).fillna("")
        distributions.to_csv(DIST_PATH, sep="\t", index=False, na_rep="")
        profiles.to_csv(PROFILE_PATH, sep="\t", index=False, na_rep="")
        print("RETAINED", DIST_PATH, "rows", len(distributions), flush=True)
        print("RETAINED", PROFILE_PATH, "rows", len(profiles), flush=True)
        print("OBSERVATION_UNITS", "44 sections / 5 true donors / 14 regions; six frozen cortical layers; 54 retained programs", flush=True)
        print("NONFINITE_BIN_SCORES", int(distributions.loc[distributions.record_type.eq("chip_layer_bins"), "n_nonfinite"].sum()), flush=True)
    render(ann, ids, distributions, profiles)


if __name__ == "__main__":
    main()
