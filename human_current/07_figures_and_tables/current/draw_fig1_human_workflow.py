#!/usr/bin/env python3
"""Approved human-only Fig1a. Retained vector PDF and publication PNG."""
from pathlib import Path
import fitz
ROOT=Path((__import__("os").environ["NMF_WORK_ROOT"] + ""))
OUT=ROOT/'figures/human_revision/panels/fig1_human_workflow'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    d=fitz.open();p=d.new_page(width=590,height=110)
    blocks=[('Human nuclei','1,036,039 nuclei\n14 areas; 10 donors',(.16,.38,.58)),('Program discovery','Consensus NMF\n54 retained programs',(.16,.45,.43)),('Human tissue','Fixed-basis projection\n44 sections; 25-µm bins',(.43,.32,.60)),('External annotation','Disease and aging\ngene-set overlap',(.58,.34,.23))]
    for i,(title,body,col) in enumerate(blocks):
        x=8+i*147;p.draw_rect(fitz.Rect(x,8,x+135,101),color=col,fill=(.97,.98,.98),width=.8)
        p.insert_text((x+7,31),title,fontsize=10,fontname='hebo',color=col)
        p.insert_text((x+7,56),body,fontsize=8.8,lineheight=1.8)
        if i<3:
            p.draw_line((x+136,56),(x+145,56),width=1,color=(.4,.4,.4))
            p.draw_line((x+142,53),(x+145,56),width=1,color=(.4,.4,.4))
            p.draw_line((x+142,59),(x+145,56),width=1,color=(.4,.4,.4))
    d.save(OUT/'Fig1_human_workflow.pdf',garbage=1,deflate=True)
    p.get_pixmap(dpi=300,alpha=False).save(OUT/'Fig1_human_workflow.png')
    print('Produced',OUT/'Fig1_human_workflow.pdf',flush=True)
if __name__=='__main__':main()
