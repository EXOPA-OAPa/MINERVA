"""
Generate the flowchart of MINERVA as a vector image (SVG), used in README.md and docs/.

Usage (from the repository root):
    python scripts/make_flowchart.py docs/images/minerva_flowchart.svg               # complete diagram
    python scripts/make_flowchart.py docs/images/minerva_simulation.svg simulation   # left panel only
    python scripts/make_flowchart.py docs/images/minerva_recovery.svg recovery       # right panel only

Edit the boxes below when the code changes; the layout is computed from them.
"""
import html
import sys
from pathlib import Path

LOGO = Path(__file__).resolve().parents[1] / "docs" / "images" / "minerva_logo.svg"
MODE = sys.argv[2] if len(sys.argv) > 2 else "full"     # "full", "simulation" or "recovery"
assert MODE in ("full", "simulation", "recovery"), MODE
TOP = 62 if MODE == "full" else 0   # vertical space added for the logo

W, H = 1500, 1475 + 62
FONT = "Helvetica, Arial, 'DejaVu Sans', sans-serif"
out = []


def esc(s):
    return html.escape(s, quote=False)


def text(x, y, s, size=13, weight="normal", fill="#222", anchor="middle", style="normal", family=FONT):
    out.append(f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" font-weight="{weight}" '
               f'font-style="{style}" fill="{fill}" text-anchor="{anchor}">{esc(s)}</text>')


def box(cx, y, w, lines, fill, stroke, h=None, rx=9, title_size=14.5, sub_size=12):
    """Rounded box centred at cx, top at y. First line bold, the others smaller."""
    h = h or (22 + 17 * len(lines))
    x = cx - w / 2
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="1.6"/>')
    n = len(lines)
    total = 18 + 16 * (n - 1)
    ty = y + (h - total) / 2 + 13
    text(cx, ty, lines[0], size=title_size, weight="bold")
    for i, s in enumerate(lines[1:]):
        text(cx, ty + 18 + 16 * i, s, size=sub_size, fill="#3a3a3a")
    return dict(cx=cx, x=x, y=y, w=w, h=h, top=y, bottom=y + h, left=x, right=x + w, cy=y + h / 2)


def diamond(cx, cy, w, h, lines, fill="#FEF9E7", stroke="#B7950B"):
    pts = f"{cx},{cy - h/2} {cx + w/2},{cy} {cx},{cy + h/2} {cx - w/2},{cy}"
    out.append(f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="1.6"/>')
    ty = cy - 8 * (len(lines) - 1) + 4.5
    for i, s in enumerate(lines):
        text(cx, ty + 16 * i, s, size=13, weight="bold")
    return dict(cx=cx, cy=cy, top=cy - h / 2, bottom=cy + h / 2, left=cx - w / 2, right=cx + w / 2)


def arrow(points, color="#444", dashed=False, width=1.6):
    d = " ".join(f"{x},{y}" for x, y in points)
    dash = ' stroke-dasharray="6,4"' if dashed else ""
    out.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="{width}"{dash} marker-end="url(#arrow)"/>')


def label(x, y, s, anchor="middle", color="#555", size=12, weight="bold"):
    text(x, y, s, size=size, weight=weight, fill=color, anchor=anchor)


def frame(x, y, w, h, stroke, title, sub, dash=True):
    d = ' stroke-dasharray="8,5"' if dash else ""
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="none" stroke="{stroke}" stroke-width="1.8"{d}/>')
    text(x + 14, y + 22, title, size=14, weight="bold", fill=stroke, anchor="start")
    text(x + w - 12, y + 22, sub, size=12, fill=stroke, anchor="end", family="'DejaVu Sans Mono', Menlo, monospace")


# ----------------------------------------------------------------------------------------
out.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
out.append('<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
           'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#444"/></marker></defs>')
out.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# Vector logo embedded in the diagram (an SVG shown with <img> cannot load external files)
if MODE == "full":
  logo = LOGO.read_text()
  inner = logo[logo.index('>', logo.index('<svg')) + 1:logo.rindex('</svg>')]
  inner = inner.replace('<title>MINERVA</title>', '')
  logo_h = 100
  k = logo_h / 263
  out.append(f'<g transform="translate({(W - 522 * k) / 2:.1f},10) scale({k:.4f})">{inner}</g>')
  text(W / 2, 128, "one Monte Carlo simulation of a target: synthetic star and RV time series, then sequential recovery",
       size=14, fill="#555")
out.append(f'<g transform="translate(0,{TOP})">')

# ---------------- left panel: simulation ----------------
LP = dict(x=20, y=90, w=715, h=1365)
out.append(f'<rect x="{LP["x"]}" y="{LP["y"]}" width="{LP["w"]}" height="{LP["h"]}" rx="16" fill="#FBFCFC" stroke="#D5D8DC" stroke-width="1.5"/>')
text(LP["x"] + LP["w"] / 2, LP["y"] + 32, "1 · Synthetic star and RV time series", size=19, weight="bold", fill="#1B2631")

# star chain
SX = 200
frame(36, 140, 330, 690, "#B03A2E", "Synthetic star", "star.py")
STAR = dict(fill="#FDEDEC", stroke="#B03A2E")
tp = box(SX, 180, 290, ["Target parameters", "catalog row or single_target_params", "optional bootstrap within the errors"], h=72, **STAR)
text(SX, tp["bottom"] + 22, "each step is skipped when its value is given", size=12, style="italic", fill="#7B241C")
chain = [
    ["Age", "uniform, up to the end of the main sequence"],
    ["Rotation regime", "slow 25% · median 50% · fast 25%"],
    ["Rotation period  P_rot", "Gallet & Bouvier (2015)"],
    ["Turnover time τ,  Ro = P_rot / τ", "YaPSI tracks (Spada et al. 2017)"],
    ["Activity level  log R'HK", "Mamajek & Hillenbrand (2008)", "M dwarfs: Astudillo-Defru et al. (2017)"],
    ["Activity amplitude  K_act", "Santos et al. (2000)", "K dwarfs: GAPS/ArMS relation"],
]
y = tp["bottom"] + 38
nodes = []
for c in chain:
    h = 56 if len(c) == 2 else 72
    n = box(SX, y, 290, c, h=h, **STAR)
    if nodes:
        arrow([(SX, nodes[-1]["bottom"]), (SX, n["top"])])
    else:
        arrow([(SX, tp["bottom"] + 30), (SX, n["top"])])
    nodes.append(n)
    y = n["bottom"] + 30
arrow([(SX, tp["bottom"]), (SX, tp["bottom"] + 8)])
kact = nodes[-1]

# signals
GX = 553
frame(384, 140, 338, 690, "#1E8449", "RV signals", "signals.py")
SIG = dict(fill="#E9F7EF", stroke="#1E8449")
cad = box(GX, 180, 306, ["Observing cadence", "real epochs (weather, seasons)", "USP sampling for P < 2 d"], h=72, **SIG)
wn = box(GX, cad["bottom"] + 46, 306, ["White noise", "σ = √(σ_instr² + σ_photon² + σ_planets²)"], h=56, **SIG)
act = box(GX, wn["bottom"] + 46, 306, ["Activity signal", "P_rot and P_rot (1 + α): differential rotation",
                                       "new phases and amplitudes every season", "harmonic P_rot / 2 for high activity"], h=90, **SIG)
pl = box(GX, act["bottom"] + 46, 306, ["Planet (Keplerian)", "K_pl from the mass, P and T0 from the transit,", "eccentricity e, periastron ω"],
         h=72, fill="#EBF5FB", stroke="#2471A3")
eb = box(GX, pl["bottom"] + 46, 306, ["Error bars", "√(σ_instr² + σ_photon² + σ_activity²)", "activity term from the GAPS database"], h=72, **SIG)

# star -> activity signal and error bars
arrow([(kact["right"], kact["cy"]), (372, kact["cy"]), (372, act["cy"]), (act["left"], act["cy"])], color="#B03A2E")
logr = nodes[-2]
arrow([(logr["right"], logr["cy"] + 10), (378, logr["cy"] + 10), (378, eb["cy"] + 12), (eb["left"], eb["cy"] + 12)], color="#B03A2E", dashed=True)

# RV time series
rv = box(LP["x"] + LP["w"] / 2, 875, 470, ["Synthetic RV time series",
                                           "RV(t) = white noise + activity + planet, with error bars"],
         h=62, fill="#FCF3CF", stroke="#B7950B")
arrow([(GX, 830), (GX, rv["top"])])
arrow([(SX, 830), (SX, rv["top"])])

mc = box(LP["x"] + LP["w"] / 2, 1000, 640, ["Monte Carlo",
                                           "n_iterations simulations per target, one random generator per target (random_seed + i)",
                                           "→ median and 16th–84th percentiles of N_RV, failures, for every SNR threshold"],
         h=78, fill="#F2F3F4", stroke="#7F8C8D")
text(LP["x"] + LP["w"] / 2, mc["bottom"] + 26, "catalog targets can run in parallel (n_workers)", size=12.5, style="italic", fill="#566573")
refd = box(LP["x"] + LP["w"] / 2, mc["bottom"] + 52, 640, ["Reference data  (data/)", "HARPS-N photon noise (NSTS) · maximum stellar ages · GAPS activity database",
           "YaPSI evolutionary tracks · default cadences (HD 164922, TOI-1807)"], h=78, fill="#F2F3F4", stroke="#7F8C8D")
outb = box(LP["x"] + LP["w"] / 2, refd["bottom"] + 24,  640, ["Outputs", "results_<target>.txt: one line per simulation",
           "results_<catalog>.txt: statistics per target   ·   Log_*.txt   ·   Plots/*.pdf"], h=78, fill="#F2F3F4", stroke="#7F8C8D")
text(LP["x"] + LP["w"] / 2, outb["bottom"] + 30, "runner.py · simulation.py · summary.py", size=12, fill="#566573",
     family="'DejaVu Sans Mono', Menlo, monospace")

# ---------------- right panel: recovery ----------------
RP = dict(x=755, y=90, w=725, h=1365)
out.append(f'<rect x="{RP["x"]}" y="{RP["y"]}" width="{RP["w"]}" height="{RP["h"]}" rx="16" fill="#FBFCFC" stroke="#D5D8DC" stroke-width="1.5"/>')
text(RP["x"] + RP["w"] / 2, RP["y"] + 32, "2 · Sequential recovery", size=19, weight="bold", fill="#1B2631")
text(RP["x"] + RP["w"] / 2, RP["y"] + 54, "periodogram.py · fitting.py · detection.py", size=12, fill="#566573",
     family="'DejaVu Sans Mono', Menlo, monospace")

CX, SX2 = 975, 1300
REC = dict(fill="#E8F6F3", stroke="#117A65")
r1 = box(CX, 162, 330, ["First n = n_start observations", "in time order (sampling_strategy)"], h=56, **REC)
r2 = box(CX, r1["bottom"] + 30, 350, ["GLS periodogram", "peaks above the FAP thresholds"], h=56, **REC)
arrow([(CX, r1["bottom"]), (CX, r2["top"])])
d1 = diamond(CX, r2["bottom"] + 70, 220, 92, ["Intruder", "peaks?"])
arrow([(CX, r2["bottom"]), (CX, d1["top"])])

r3 = box(SX2, d1["cy"] - 33, 300, ["Pre-whitening", "fit and subtract the intruders", "in order of decreasing power"], h=66, **REC)
arrow([(d1["right"], d1["cy"]), (r3["left"], d1["cy"])]); label(d1["right"] + 20, d1["cy"] - 8, "yes")
d2 = diamond(SX2, r3["bottom"] + 72, 230, 92, ["Planet peak", "improves?"])
arrow([(SX2, r3["bottom"]), (SX2, d2["top"])])

d3 = diamond(CX, d2["cy"] + 70, 220, 92, ["Planet is the", "highest peak?"])
arrow([(CX, d1["bottom"]), (CX, d3["top"])]); label(CX + 14, d1["bottom"] + 18, "no", anchor="start")
arrow([(d2["left"], d2["cy"]), (CX, d2["cy"]), (CX, d3["top"])])
label(d2["left"] - 14, d2["cy"] - 8, "none accepted", anchor="end")

r4 = box(CX, d3["bottom"] + 30, 330, ["Reverse pre-whitening", "subtract the planet, GLS of the residuals"], h=56, **REC)
arrow([(CX, d3["bottom"]), (CX, r4["top"])]); label(CX + 14, d3["bottom"] + 18, "yes", anchor="start")
d4 = diamond(CX, r4["bottom"] + 66, 230, 92, ["Stellar peak", "above the FAP?"])
arrow([(CX, r4["bottom"]), (CX, d4["top"])])

sim = box(SX2, d4["cy"] - 45, 300, ["Simultaneous fit", "planet + sinusoid at P_star,", "amplitude and phase per season", "(planet-only fit if it fails)"], h=90, **REC)
arrow([(SX2, d2["bottom"]), (SX2, sim["top"])]); label(SX2 + 14, d2["bottom"] + 18, "yes", anchor="start")
arrow([(d4["right"], d4["cy"]), (sim["left"], d4["cy"])]); label(d4["right"] + 20, d4["cy"] - 8, "yes")

po = box(CX, d4["bottom"] + 36, 290, ["Planet-only fit", "circular orbit, P and T0 fixed"], h=56, **REC)
arrow([(CX, d4["bottom"]), (CX, po["top"])]); label(CX + 14, d4["bottom"] + 20, "no", anchor="start")
arrow([(d3["left"], d3["cy"]), (812, d3["cy"]), (812, po["cy"]), (po["left"], po["cy"])])
label(d3["left"] - 14, d3["cy"] - 8, "no", anchor="end")

st = box(CX, po["bottom"] + 36, 350, ["Fit statistics", "χ²_red,  σ_K × max(1, √χ²_red),  SNR = K / σ_K", "ΔBIC against the null hypothesis"], h=74, **REC)
arrow([(CX, po["bottom"]), (CX, st["top"])])
arrow([(SX2, sim["bottom"]), (SX2, st["cy"]), (st["right"], st["cy"])])

crit = box(CX, st["bottom"] + 30, 350, ["Detection criteria", "ΔBIC > delta_bic_threshold,  K > 0",
          "K consistent with the injected K (3σ)", "P_star = P_rot, P_rot/2 or P_rot/3 (15%)", "K / σ_K above the SNR thresholds"],
          h=104, fill="#FEF9E7", stroke="#B7950B")
arrow([(CX, st["bottom"]), (CX, crit["top"])])
d5 = diamond(CX, crit["bottom"] + 58, 230, 84, ["All met at the", "highest SNR?"])
arrow([(CX, crit["bottom"]), (CX, d5["top"])])
end = box(SX2, d5["cy"] - 34, 290, ["N_RV for every SNR threshold", "detection (status 1)"], h=68, fill="#FDEBD0", stroke="#CA6F1E")
arrow([(d5["right"], d5["cy"]), (end["left"], d5["cy"])]); label(d5["right"] + 20, d5["cy"] - 8, "yes")

d6 = diamond(CX, d5["bottom"] + 60, 230, 84, ["Observations", "left?"])
arrow([(CX, d5["bottom"]), (CX, d6["top"])]); label(CX + 14, d5["bottom"] + 18, "no", anchor="start")
fail = box(SX2, d6["cy"] - 34, 290, ["Observations exhausted", "N_RV = N_total (status 0 or 3)"], h=68, fill="#F4ECF7", stroke="#7D3C98")
arrow([(d6["right"], d6["cy"]), (fail["left"], d6["cy"])]); label(d6["right"] + 20, d6["cy"] - 8, "no")

add_y = d6["bottom"] + 26
add = box(CX, add_y, 300, ["Add `step` observations", "(5 by default)"], h=50, fill="#E8F6F3", stroke="#117A65")
arrow([(CX, d6["bottom"]), (CX, add["top"])]); label(CX + 14, d6["bottom"] + 17, "yes", anchor="start")
arrow([(add["left"], add["cy"]), (783, add["cy"]), (783, r2["cy"]), (r2["left"], r2["cy"])], color="#117A65")

# link between the panels
if MODE == "full":
    arrow([(rv["right"], rv["cy"]), (745, rv["cy"]), (745, r1["cy"]), (r1["left"], r1["cy"])], color="#B7950B", width=2.2)
if MODE == "simulation":
    arrow([(rv["right"], rv["cy"]), (rv["right"] + 50, rv["cy"])], color="#B7950B", width=2.2)
    text(rv["right"] + 2, rv["bottom"] + 20, "→ sequential recovery", size=12.5, style="italic", fill="#9A7D0A", anchor="start")
if MODE == "recovery":
    arrow([(772, r1["cy"]), (r1["left"], r1["cy"])], color="#B7950B", width=2.2)
    text(770, r1["cy"] - 12, "RV series", size=12, style="italic", fill="#9A7D0A", anchor="start")

out.append("</g>")
out.append("</svg>")
if MODE != "full":
    vx = 12 if MODE == "simulation" else 747
    vw, vy, vh = 731, 82, 1381
    out[0] = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vx} {vy} {vw} {vh}" width="{vw}" height="{vh}">'
    out[2] = f'<rect x="{vx}" y="{vy}" width="{vw}" height="{vh}" fill="white"/>'
assert add["bottom"] < H - 10, add["bottom"]
target = sys.argv[1] if len(sys.argv) > 1 else "docs/images/minerva_flowchart.svg"
open(target, "w", encoding="utf-8").write("\n".join(out))
print("written", target)
