"""
Generate the overview diagram of MINERVA (docs/images/minerva_overview.svg), shown in the
README: the main stages of one Monte Carlo simulation. The detailed diagrams are made by
scripts/make_flowchart.py.

Usage (from the repository root):
    python scripts/make_overview.py docs/images/minerva_overview.svg
"""
import html
import sys

W, H = 1100, 400
FONT = "Helvetica, Arial, 'DejaVu Sans', sans-serif"
MONO = "'DejaVu Sans Mono', Menlo, Consolas, monospace"
out = []


def text(x, y, s, size=14, weight="normal", fill="#222", anchor="middle", style="normal", family=FONT):
    out.append(f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" font-weight="{weight}" '
               f'font-style="{style}" fill="{fill}" text-anchor="{anchor}">{html.escape(s, quote=False)}</text>')


def box(x, y, w, h, title, lines, fill, stroke, module=None):
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}" stroke="{stroke}" stroke-width="2"/>')
    cx = x + w / 2
    mod = module or []
    n = len(lines)
    block = 24 + 21 * n + (10 + 17 * len(mod) if mod else 0)
    top = y + (h - block) / 2 + 17
    text(cx, top, title, size=17.5, weight="bold", fill="#1B2631")
    for i, s in enumerate(lines):
        text(cx, top + 27 + 21 * i, s, size=14.5, fill="#34495E")
    for i, m in enumerate(mod):
        text(cx, top + 27 + 21 * n + 10 + 17 * i, m, size=12.5, fill=stroke, family=MONO)
    return dict(x=x, y=y, w=w, h=h, cx=cx, cy=y + h / 2, right=x + w, bottom=y + h)


def arrow(points, color="#555", width=2.2):
    d = " ".join(f"{x},{y}" for x, y in points)
    out.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="{width}" marker-end="url(#arrow)"/>')


out.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
out.append('<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
           'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#555"/></marker></defs>')
out.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# Monte Carlo frame
out.append(f'<rect x="12" y="12" width="{W - 24}" height="{H - 24}" rx="18" fill="#FBFCFC" stroke="#AAB7B8" '
           f'stroke-width="1.6" stroke-dasharray="9,6"/>')
text(32, 42, "Monte Carlo: n_iterations simulations per target", size=16, weight="bold", fill="#566573", anchor="start")
text(W - 30, 42, "runner.py · simulation.py", size=13, fill="#566573", anchor="end", family=MONO)

Y, BH = 70, 196
b0 = box(30, Y, 150, BH, "Target", ["catalog row", "or single target", "(config.yaml)"], "#F2F3F4", "#7F8C8D")
b1 = box(212, Y, 200, BH, "1 · Synthetic star", ["age, rotation", "turnover time τ", "log R'HK, K_act"],
         "#FDEDEC", "#B03A2E", ["star.py"])
b2 = box(444, Y, 200, BH, "2 · RV time series", ["real cadence, noise", "activity signal", "planet, error bars"],
         "#E9F7EF", "#1E8449", ["signals.py"])
b3 = box(676, Y, 228, BH, "3 · Sequential recovery", ["GLS, pre-whitening", "planet or planet +", "activity fit, criteria"],
         "#E8F6F3", "#117A65", ["periodogram.py", "fitting.py · detection.py"])
b4 = box(936, Y, 134, BH, "N_RV", ["for each SNR", "threshold"], "#FDEBD0", "#CA6F1E")

for a, b in [(b0, b1), (b1, b2), (b2, b3), (b3, b4)]:
    arrow([(a["right"] + 2, a["cy"]), (b["x"] - 2, b["cy"])])

# loop of the sequential recovery
ly = b3["bottom"] + 34
arrow([(b3["right"] - 30, b3["bottom"] + 2), (b3["right"] - 30, ly), (b3["x"] + 30, ly), (b3["x"] + 30, b3["bottom"] + 2)],
      color="#117A65")
text(b3["cx"], ly + 22, "+ step observations, until the planet is", size=13.5, style="italic", fill="#117A65")
text(b3["cx"], ly + 40, "detected or the observations run out", size=13.5, style="italic", fill="#117A65")

text(36, H - 34, "→ per target: median and 16th–84th percentiles of N_RV over the simulations, and the failures",
     size=14.5, fill="#34495E", anchor="start")

out.append("</svg>")
target = sys.argv[1] if len(sys.argv) > 1 else "docs/images/minerva_overview.svg"
open(target, "w", encoding="utf-8").write("\n".join(out))
print("written", target)
