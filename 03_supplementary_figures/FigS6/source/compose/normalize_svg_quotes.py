#!/usr/bin/env python
"""Normalize svglite single-quote attrs to double-quotes so ink_crop's regex parses them.
Only rewrites attribute quotes on the <svg ...> opening tag (enough for the cropper),
but to be safe we convert all single-quoted attr values to double on the whole file
(SVG content uses no literal single-quote text in svglite output)."""
import sys, glob, os, re
d = sys.argv[1]
# ONLY convert single->double quotes inside the opening <svg ...> tag (no nested quotes there).
# Leave style='...font-family: "X"...' bodies untouched (nested-quote conflict otherwise).
def fix_svg_tag(t):
    m = re.search(r"<svg\b[^>]*>", t)
    if not m: return t
    tag = m.group(0)
    tag2 = re.sub(r"=\'([^\']*)\'", r'="\1"', tag)
    return t[:m.start()] + tag2 + t[m.end():]
for f in glob.glob(os.path.join(d, "*.svg")):
    t = open(f, encoding="utf-8").read()
    open(f, "w", encoding="utf-8").write(fix_svg_tag(t))
    print("normalized", os.path.basename(f))
