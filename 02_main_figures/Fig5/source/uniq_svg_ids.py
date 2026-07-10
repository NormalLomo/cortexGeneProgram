#!/usr/bin/env python3
# Namespace every id="..."/url(#...) in each figB_*.svg by a per-panel prefix so
# that when svgutils merges the panels into one composite, svglite's deterministic
# (geometry-hashed) clipPath IDs cannot COLLIDE across panels. Without this, two
# panels that share an identical clip-rect emit the same base64 id; in the merged
# SVG the later definition wins and clips the earlier panel's content to a sliver
# (symptom: the title band rendered as a single stray "(" glyph). Pure-text rewrite
# of id tokens only; geometry/positions untouched.
import os, re, glob

SVGD = "CORTEX_PROGRAM_ROOT/scripts/figmarkcorr_B/svg_panels"

ID_DEF = re.compile(r'\bid="([^"]+)"')
URL_REF = re.compile(r'url\(#([^)]+)\)')
HREF_REF = re.compile(r'(xlink:href|href)="#([^"]+)"')


def uniq_one(path, prefix):
    s = open(path, encoding="utf-8").read()
    ids = set(ID_DEF.findall(s))
    for old in sorted(ids, key=len, reverse=True):
        new = f"{prefix}_{old}"
        s = s.replace(f'id="{old}"', f'id="{new}"')
        s = s.replace(f'url(#{old})', f'url(#{new})')
        s = s.replace(f'href="#{old}"', f'href="#{new}"')
    open(path, "w", encoding="utf-8").write(s)
    return len(ids)


def main():
    for p in sorted(glob.glob(os.path.join(SVGD, "figB_*.svg"))):
        pid = os.path.basename(p).replace("figB_", "").replace(".svg", "")
        n = uniq_one(p, f"p{pid}")
        print(f"{os.path.basename(p)}: namespaced {n} ids -> prefix p{pid}_")


if __name__ == "__main__":
    main()
