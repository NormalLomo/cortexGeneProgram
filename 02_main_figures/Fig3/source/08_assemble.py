"""Assemble fig2_full.pdf from the 8 per-panel vector PDFs.
Rasterize each panel via pdftoppm (high DPI), compose into a Nature-style
multi-panel layout with bold panel letters. Per-panel PDFs remain vector masters."""
import subprocess, os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import matplotlib.gridspec as gridspec
from PIL import Image

outd="CORTEX_PROGRAM_ROOT/figures/fig2/"
tmp="CORTEX_PROGRAM_ROOT/scripts/fig2/_raster/"
os.makedirs(tmp, exist_ok=True)
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","DejaVu Sans"]})

panels="abcdefgh"
DPI=400
imgs={}
for p in panels:
    src=outd+f"fig2_{p}.pdf"
    base=tmp+f"p{p}"
    subprocess.run(["pdftoppm","-r",str(DPI),"-png","-singlefile",src,base],check=True)
    im=Image.open(base+".png")
    imgs[p]=np.asarray(im)
    print(p, imgs[p].shape)

def place(ax, p, letter):
    ax.imshow(imgs[p]); ax.axis("off")
    h,w=imgs[p].shape[:2]
    ax.set_xlim(0,w); ax.set_ylim(h,0)
    ax.text(-0.02, 1.01, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="right")

# Layout matched to native panel aspects. e is tall -> give it a full-height left column.
# Panel e is now a wide 2-col split (half height) -> gets its own wide row.
# Left column (a small + c) rebalanced against the taller b on the right of row1-2.
fig=plt.figure(figsize=(14.0, 18.0))
gs=gridspec.GridSpec(6, 6, figure=fig, hspace=0.10, wspace=0.06,
    height_ratios=[1.30,1.05,1.05,0.95,1.00,0.95])

# row1: a (cols 0-2, enlarged so layer key + chip ID legible) | b (cols 2-6)
ax_a=fig.add_subplot(gs[0,0:2]); place(ax_a,"a","a")
ax_b=fig.add_subplot(gs[0,2:6]); place(ax_b,"b","b")
# row2: c (cols 0-3, slightly wider) | d (cols 3-6)
ax_c=fig.add_subplot(gs[1,0:3]); place(ax_c,"c","c")
ax_d=fig.add_subplot(gs[1,3:6]); place(ax_d,"d","d")
# row3: e wide 2-col heatmap, spans most of the width
ax_e=fig.add_subplot(gs[2,0:5]); place(ax_e,"e","e")
# row4: f wide
ax_f=fig.add_subplot(gs[3,0:6]); place(ax_f,"f","f")
# row5: g wide
ax_g=fig.add_subplot(gs[4,0:6]); place(ax_g,"g","g")
# row6: h (centered)
ax_h=fig.add_subplot(gs[5,1:5]); place(ax_h,"h","h")

fig.suptitle("Figure 2 | Spatial program backbone of human cortex (bin50, 50 µm)",
             fontsize=17, y=0.997, x=0.012, ha="left", fontweight="bold")
fig.savefig(outd+"fig2_full.pdf", dpi=300, bbox_inches="tight", pad_inches=0.15)
fig.savefig(outd+"fig2_full.png", dpi=200, bbox_inches="tight", pad_inches=0.15)
print("ASSEMBLED")
