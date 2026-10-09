"""
Vector MINERVA logo (SVG), rebuilt from the measured geometry of the original raster logo.

Letters: Montserrat (SIL Open Font License 1.1), weight 775, converted to paths, so the
logo does not depend on installed fonts. The M is drawn as one outline: its second diagonal
is the thick start of the RV curve. Coordinates are those of the original 522 x 263 px image.

Usage (from the repository root; needs fontTools and the Montserrat variable font,
https://github.com/google/fonts/tree/main/ofl/montserrat):
    python scripts/make_logo.py docs/images/minerva_logo.svg Montserrat[wght].ttf

The font file may be the variable font (the weight is set here) or a static instance.
"""
import sys

import numpy as np
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

OUT = sys.argv[1] if len(sys.argv) > 1 else 'docs/images/minerva_logo.svg'
FONT = sys.argv[2] if len(sys.argv) > 2 else 'Montserrat[wght].ttf'
WEIGHT = 775

BASE, CAP = 167.0, 68.0                       # baseline and cap height [px]
LEFT = {'I': 102, 'N': 134, 'E': 212, 'R': 278, 'V': 340}
A_RIGHT = 479                                 # right edge of the A (V and A overlap)
TEXT_STOPS = [(9, '#2252A6'), (50, '#2673B7'), (80, '#1E8DC9'), (140, '#1EB3D2'), (195, '#48AF9D'),
              (250, '#67A972'), (290, '#9AAA5C'), (320, '#D6AC3D'), (350, '#E59C39'), (410, '#D66041'),
              (479, '#D43141')]
WAVE_STOPS = [(47, '#206FB9'), (79, '#2B82B1'), (100, '#4191A9'), (130, '#509B93'), (180, '#4A9479'),
              (300, '#559A7A')]
GREEN, DOT = '#4F9578', '#4B9278'
LINE = 7.2                                    # width of the curve and of the orbit
CX, CY, RX, RY = 391.0, 130.0, 117.6, 120.9    # orbit (slightly elliptical, as in the original)
DOTS = [(117, 8), (182, 85), (244, 8), (290, 64)]
DOT_R = 7.3
PLANET = (475.5, 46.5, 10.6)                  # x, y, radius

font = TTFont(FONT)
if 'fvar' in font:                                  # variable font: instance at WEIGHT
    from fontTools.varLib import instancer
    font = instancer.instantiateVariableFont(font, {'wght': WEIGHT})
gs, cmap = font.getGlyphSet(), font.getBestCmap()
scale = CAP / font['OS/2'].sCapHeight


def glyph_path(ch, left=None, right=None):
    from fontTools.pens.boundsPen import BoundsPen
    g = gs[cmap[ord(ch)]]
    bp = BoundsPen(gs); g.draw(bp); x0, _, x1, _ = bp.bounds
    dx = (left - x0 * scale) if left is not None else (right - x1 * scale)
    pen = SVGPathPen(gs, ntos=lambda v: f'{v:.2f}'.rstrip('0').rstrip('.'))
    g.draw(TransformPen(pen, (scale, 0, 0, -scale, dx, BASE)))
    return pen.getCommands()


def stops(lst, x0, x1):
    return ''.join(f'<stop offset="{(x - x0) / (x1 - x0):.3f}" stop-color="{c}"/>' for x, c in lst)


# --- wave: four cubic Bezier arcs through the measured points of the original
#     (horizontal tangents at the two maxima and at the minimum; start along the diagonal of
#     the M; arrival in the last dot with the measured direction)
A = np.array([50.2, 159.5])                       # vertex of the M (centre of its base)
d_start = np.array([0.40, -0.92])                 # direction of the second diagonal of the M
d_end = np.array([0.47, 0.88])                    # direction of the curve entering the last dot
P1, T, P2, E = map(np.array, DOTS[:1] + DOTS[1:2] + DOTS[2:3] + DOTS[3:4])
P1, T, P2, E = np.array(DOTS[0], float), np.array(DOTS[1], float), np.array(DOTS[2], float), np.array(DOTS[3], float)
h = 27.0                                          # length of the horizontal tangents
segs = [(A, A + 95 * d_start, P1 - (17, 0), P1),
        (P1, P1 + (h, 0), T - (h, 0), T),
        (T, T + (h, 0), P2 - (h, 0), P2),
        (P2, P2 + (13, 0), E - 26 * d_end, E)]
f = lambda v: f'{v[0]:.1f},{v[1]:.1f}'
wave_d = f'M{f(A)}' + ''.join(f' C{f(c1)} {f(c2)} {f(p3)}' for _, c1, c2, p3 in segs)

# Thick start of the curve (second diagonal of the M): the first arc, widened from 15.5 px at
# the vertex to the line width at the top of the letters, offset along its normals.
def bez(seg, t):
    p0, c1, c2, p3 = seg
    return (1 - t)**3 * p0 + 3 * (1 - t)**2 * t * c1 + 3 * (1 - t) * t**2 * c2 + t**3 * p3
ts = np.linspace(0, 1, 400)
pts = np.array([bez(segs[0], t) for t in ts])
top = np.argmax(pts[:, 1] <= 100.0)
pts = pts[:top + 1]
tang = np.gradient(pts, axis=0); tang /= np.linalg.norm(tang, axis=1)[:, None]
normal = np.c_[-tang[:, 1], tang[:, 0]]
frac = (pts[0, 1] - pts[:, 1]) / (pts[0, 1] - pts[-1, 1])
half = 7.75 + (LINE / 2 - 7.75) * frac
left, right = pts + normal * half[:, None], pts - normal * half[:, None]
left[:, 1] = np.minimum(left[:, 1], A[1]); right[:, 1] = np.minimum(right[:, 1], A[1])   # flat base at the vertex
# which offset is on the left of the curve
if left[:, 0].mean() > right[:, 0].mean():
    left, right = right, left

# --- the M as one outline: left stem, first diagonal, notch, thick start of the curve, flat vertex
diag_x = lambda y: 23.5 + 0.61 * (y - 100.0)          # inner edge of the first diagonal (measured)
k = np.argmax(left[:, 0] - diag_x(left[:, 1]) > 0)     # notch: going up, the curve passes to the right of the diagonal
notch = left[k]
outer_x = lambda y: 23.0 + (41.0 - 23.0) * (y - 127.5) / (159.5 - 127.5)   # outer edge of the first diagonal
right_up = [tuple(p) for p in right[::-1] if p[1] < A[1] - 0.05]       # right edge, above the base
base_r = (right_up[-1][0], A[1])                                       # where the right edge meets the base
m_poly = ([(9, 100), (23.5, 100), tuple(notch)] + [tuple(p) for p in left[k:]] +
          right_up + [base_r, (41.0, A[1]), (23.0, 127.5), (23, 167), (9, 167)])
m_left = 'M' + ' L'.join(f(p) for p in m_poly) + ' Z'

# right stem: its top-left corner follows the right edge of the curve at a constant gap
GAP = 6.5
edge = [(x + GAP, y) for x, y in right[::-1] if x + GAP >= 72 and y <= 140]
m_right = 'M88,100 ' + ' '.join(f'L{x:.2f},{y:.2f}' for x, y in edge) + ' L72,167 L88,167 Z'

# --- orbit: from the last dot clockwise to where it disappears behind the R
t0 = np.degrees(np.arctan2((DOTS[-1][1] - CY) / RY, (DOTS[-1][0] - CX) / RX)) % 360
t1 = 153.0
x0, y0 = CX + RX * np.cos(np.radians(t0)), CY + RY * np.sin(np.radians(t0))
x1, y1 = CX + RX * np.cos(np.radians(t1)), CY + RY * np.sin(np.radians(t1))
sweep = (t1 - t0) % 360
orbit_d = f'M{x0:.1f},{y0:.1f} A{RX},{RY} 0 {1 if sweep > 180 else 0} 1 {x1:.1f},{y1:.1f}'

# --- the M: left stem + first diagonal, right stem, and the thick start of the curve

px, py, pr = PLANET
svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 522 263" width="522" height="263">
  <title>MINERVA</title>
  <defs>
    <linearGradient id="text" gradientUnits="userSpaceOnUse" x1="9" y1="0" x2="479" y2="0">{stops(TEXT_STOPS, 9, 479)}</linearGradient>
    <linearGradient id="wave" gradientUnits="userSpaceOnUse" x1="47" y1="0" x2="300" y2="0">{stops(WAVE_STOPS, 47, 300)}</linearGradient>
    <clipPath id="above-vertex"><rect x="0" y="0" width="522" height="{A[1]:.1f}"/></clipPath>
    <linearGradient id="planet" x1="0.15" y1="0.85" x2="0.85" y2="0.15">
      <stop offset="0" stop-color="#2E7576"/><stop offset="1" stop-color="#58A57E"/></linearGradient>
  </defs>
  <path d="{orbit_d}" fill="none" stroke="{GREEN}" stroke-width="{LINE}" stroke-linecap="round"/>
  <path d="{wave_d}" fill="none" stroke="url(#wave)" stroke-width="{LINE}" stroke-linecap="butt" stroke-linejoin="round" clip-path="url(#above-vertex)"/>
  <g fill="url(#text)">
    <path d="{m_left}"/><path d="{m_right}"/>
    <path d="{glyph_path('I', left=LEFT['I'])}"/>
    <path d="{glyph_path('N', left=LEFT['N'])}"/>
    <path d="{glyph_path('E', left=LEFT['E'])}"/>
    <path d="{glyph_path('R', left=LEFT['R'])}"/>
    <path d="{glyph_path('V', left=LEFT['V'])}"/>
    <path d="{glyph_path('A', right=A_RIGHT)}"/>
  </g>
  <g fill="{DOT}">{''.join(f'<circle cx="{x}" cy="{y}" r="{DOT_R}"/>' for x, y in DOTS)}</g>
  <g transform="translate({px},{py}) rotate(-34)">
    <path d="M-20,0 A20,4.5 0 0 1 20,0" fill="none" stroke="#3E8A6A" stroke-width="2.8"/>
    <circle r="{pr}" fill="url(#planet)" transform="rotate(34)"/>
    <path d="M-20,0 A20,4.5 0 0 0 20,0" fill="none" stroke="#5BA683" stroke-width="2.8"/>
  </g>
</svg>
'''
open(OUT, 'w').write(svg)
print('written', OUT, len(svg) // 1024, 'kB')
