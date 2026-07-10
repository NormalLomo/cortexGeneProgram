#!/usr/bin/env python3
"""Render the public README workflow image without service or backend content."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "Fig1A_workflow_overview.png"


def box(axis: plt.Axes, x: float, title: str, detail: str, color: str) -> None:
    axis.add_patch(FancyBboxPatch((x, 0.35), 0.17, 0.32, boxstyle="round,pad=0.012,rounding_size=0.015", facecolor=color, edgecolor="#263238", linewidth=1.2))
    axis.text(x + 0.085, 0.57, title, ha="center", va="center", fontsize=11, weight="bold", color="#102027")
    axis.text(x + 0.085, 0.43, detail, ha="center", va="center", fontsize=8.5, color="#263238")


def main() -> None:
    figure, axis = plt.subplots(figsize=(14, 3.8), dpi=180)
    figure.patch.set_facecolor("#f7faf9")
    figure.patch.set_alpha(1)
    axis.set_facecolor("#f7faf9")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    boxes = [
        (0.02, "Source records", "Human, mouse, and\nprimate inputs", "#d7eef2"),
        (0.22, "Raw UMI", "Validated count matrices\nand gene identifiers", "#dfead5"),
        (0.42, "cNMF", "100 replicates per\nconfigured K", "#f6e7b8"),
        (0.62, "Derived analyses", "Regional, spatial, and\ncross-species summaries", "#ead9ee"),
        (0.82, "Final assets", "Figures, tables, and\nsupplementary data", "#f5d6ca"),
    ]
    for x, title, detail, color in boxes:
        box(axis, x, title, detail, color)
    for left, right in zip(boxes, boxes[1:]):
        axis.add_patch(FancyArrowPatch((left[0] + 0.175, 0.51), (right[0] - 0.012, 0.51), arrowstyle="-|>", mutation_scale=15, linewidth=1.3, color="#546e7a"))
    axis.text(0.03, 0.85, "Human Cortical Gene Programs", fontsize=18, weight="bold", color="#102027")
    axis.text(0.03, 0.76, "Public analysis workflow", fontsize=11, color="#546e7a")
    figure.savefig(OUTPUT, bbox_inches="tight", facecolor="#f7faf9", transparent=False)


if __name__ == "__main__":
    main()
