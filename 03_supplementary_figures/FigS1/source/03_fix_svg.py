#!/usr/bin/env python
# Edit FigS1.svg title tspan: rebuild text + per-char x array for new shorter title.
# old text (literal in XML, & is &#x0026;):
#   "Extended Data Fig. 3  |  cNMF K-robustness of cross-region variability &#x0026; cell drivers"
# new title (single spaces around |):
#   "Figure S1 | cNMF K-robustness of cross-region variability & cell drivers"
import re, fitz

SVG = "CORTEX_PROGRAM_ROOT/figure_release/SUBMISSION_v8/supplementary_figures/FigS1.svg"
TTF = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
FONT_SIZE = 9.0

# new title as plain unicode (what should render); '&' will be XML-escaped to &#x0026; in output
new_plain = "Figure S1 | cNMF K-robustness of cross-region variability & cell drivers"

# measure per-char advance using fitz.Font (same engine as PDF render)
font = fitz.Font(fontfile=TTF)
xs = [0.0]
acc = 0.0
for ch in new_plain[:-1]:  # cumulative start positions; last char needs no following start
    w = font.char_lengths(ch, fontsize=FONT_SIZE)[0]
    acc += w
    xs.append(acc)
# format x array like original (space separated, trimmed floats)
def fmt(v):
    s = ("%.7f" % v).rstrip("0").rstrip(".")
    return s if s else "0"
x_attr = " ".join(fmt(v) for v in xs)

# XML-escape the new text for the tspan body (& -> &#x0026; to match original style)
new_xml_text = new_plain.replace("&", "&#x0026;")

with open(SVG, "r", encoding="utf-8") as f:
    svg = f.read()

# locate the unique title <text ...><tspan ... x="...">OLD</tspan></text> at L2107
# We match the tspan that contains the old title body.
old_body = "Extended Data Fig. 3  |  cNMF K-robustness of cross-region variability &#x0026; cell drivers"
assert svg.count(old_body) == 1, "old title body not unique: count=%d" % svg.count(old_body)

# Find the enclosing tspan opening tag to replace its x="..." attribute.
# Pattern: <tspan y="0" x="OLD_X">old_body</tspan>
m = re.search(r'(<tspan\b[^>]*\bx=")([^"]*)("[^>]*>)' + re.escape(old_body) + r'(</tspan>)', svg)
assert m is not None, "tspan with old title not found by regex"
new_tspan = m.group(1) + x_attr + m.group(3) + new_xml_text + m.group(4)
svg2 = svg[:m.start()] + new_tspan + svg[m.end():]

assert svg2 != svg, "no change made"
assert old_body not in svg2, "old body still present"
assert new_xml_text in svg2, "new body not present"

with open(SVG, "w", encoding="utf-8") as f:
    f.write(svg2)

print("OK FigS1.svg updated")
print("new text:", new_plain)
print("n chars:", len(new_plain), "n x-coords:", len(xs))
print("x array head:", x_attr[:80], "...")
print("x array tail:", x_attr[-60:])
