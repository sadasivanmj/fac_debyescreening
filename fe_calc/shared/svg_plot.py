"""Minimal SVG line/marker plots (standard library only; no matplotlib needed).

Series with missing values (None / NaN) are NOT connected: the line is broken there, so failed or
absent calculations never produce a valid-looking segment.
"""

from __future__ import annotations

import math
from html import escape as _esc
from pathlib import Path

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]   # colour-blind safe


def _nice_ticks(lo: float, hi: float, n: int = 6) -> list[float]:
    if hi == lo:
        hi = lo + 1.0
    raw = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=raw)
    start = math.ceil(lo / step) * step
    ticks, t = [], start
    while t <= hi + 1e-12 * abs(hi):
        ticks.append(round(t, 12))
        t += step
    return ticks


def _fmt(v: float) -> str:
    if v == 0:
        return "0"
    if abs(v) >= 1e4 or abs(v) < 1e-3:
        return f"{v:.0e}".replace("e+0", "e").replace("e-0", "e-")
    return f"{v:g}"


def line_plot(path: Path, title: str, xlabel: str, ylabel: str, series: list[dict], logx: bool = False,
              logy: bool = False, width: int = 760, height: int = 480, hlines: list[tuple[float, str]] | None = None,
              note: str = "", ylim: tuple[float, float] | None = None) -> None:
    """ylim = (low, high) clips the y range (points outside are dropped, never drawn on the edge)."""
    """series: dicts with name, x, y, and optional color, marker ('o','s','^','x', None), line (bool), dashed (bool), hollow."""
    left, right, top, bottom = 78, 190, 44, 66
    pw, ph = width - left - right, height - top - bottom
    if ylim:
        series = [{**s, "y": [y if (y is not None and math.isfinite(y) and ylim[0] <= y <= ylim[1]) else None for y in s["y"]]} for s in series]
    pts = [(x, y) for s in series for x, y in zip(s["x"], s["y"])
           if y is not None and x is not None and math.isfinite(y) and math.isfinite(x) and (not logy or y > 0) and (not logx or x > 0)]
    allhl = [v for v, _ in (hlines or [])]
    if not pts:
        pts = [(0.0, 0.0), (1.0, 1.0)]
    fx = (lambda v: math.log10(v)) if logx else (lambda v: v)
    fy = (lambda v: math.log10(v)) if logy else (lambda v: v)
    xs, ys = [fx(p[0]) for p in pts], [fy(p[1]) for p in pts] + [fy(v) for v in allhl if (not logy or v > 0)]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    padx, pady = (x1 - x0) * 0.04 or 0.5, (y1 - y0) * 0.06 or 0.5
    x0, x1, y0, y1 = x0 - padx, x1 + padx, y0 - pady, y1 + pady
    sx = lambda v: left + (fx(v) - x0) / (x1 - x0) * pw
    sy = lambda v: top + ph - (fy(v) - y0) / (y1 - y0) * ph
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" font-family="sans-serif" font-size="12">',
         f'<rect width="{width}" height="{height}" fill="white"/>',
         f'<text x="{left}" y="24" font-size="15" font-weight="bold">{_esc(title)}</text>',
         f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>']
    # ticks
    def decades(lo, hi):
        ks = list(range(math.ceil(lo), math.floor(hi) + 1))
        step = max(1, math.ceil(len(ks) / 7))          # at most ~7 labelled decades
        return [10.0 ** k for k in ks[::step]]
    xt = decades(x0, x1) if logx else _nice_ticks(x0, x1)
    yt = decades(y0, y1) if logy else _nice_ticks(y0, y1)
    for t in xt:
        X = sx(t)
        o += [f'<line x1="{X:.1f}" y1="{top}" x2="{X:.1f}" y2="{top + ph}" stroke="#ddd"/>',
              f'<text x="{X:.1f}" y="{top + ph + 16}" text-anchor="middle">{_fmt(t)}</text>']
    for t in yt:
        Y = sy(t)
        o += [f'<line x1="{left}" y1="{Y:.1f}" x2="{left + pw}" y2="{Y:.1f}" stroke="#ddd"/>',
              f'<text x="{left - 6}" y="{Y + 4:.1f}" text-anchor="end">{_fmt(t)}</text>']
    o += [f'<text x="{left + pw / 2}" y="{height - 28}" text-anchor="middle">{_esc(xlabel)}</text>',
          f'<text transform="translate(16 {top + ph / 2}) rotate(-90)" text-anchor="middle">{_esc(ylabel)}</text>']
    for v, label in hlines or []:
        if (not logy) or v > 0:
            Y = sy(v)
            o += [f'<line x1="{left}" y1="{Y:.1f}" x2="{left + pw}" y2="{Y:.1f}" stroke="#777" stroke-dasharray="2 4"/>',
                  f'<text x="{left + pw - 4}" y="{Y - 4:.1f}" text-anchor="end" fill="#555">{_esc(label)}</text>']
    for k, s in enumerate(series):
        col = s.get("color", PALETTE[k % len(PALETTE)])
        marker, line, dashed = s.get("marker", "o"), s.get("line", True), s.get("dashed", False)
        seg: list[tuple[float, float]] = []
        segs = []
        for x, y in zip(s["x"], s["y"]):
            ok = y is not None and x is not None and math.isfinite(y) and math.isfinite(x) and (not logy or y > 0) and (not logx or x > 0)
            if ok:
                seg.append((sx(x), sy(y)))
            else:
                segs.append(seg)
                seg = []
        segs.append(seg)
        if line:
            for sg in segs:
                if len(sg) > 1:
                    o.append(f'<polyline fill="none" stroke="{col}" stroke-width="1.8"' + (' stroke-dasharray="6 4"' if dashed else '')
                             + f' points="{" ".join(f"{a:.1f},{b:.1f}" for a, b in sg)}"/>')
        for sg in segs:
            for a, b in sg:
                if marker == "o":
                    o.append(f'<circle cx="{a:.1f}" cy="{b:.1f}" r="3.6" fill="{"white" if s.get("hollow") else col}" stroke="{col}" stroke-width="1.5"/>')
                elif marker == "s":
                    o.append(f'<rect x="{a - 3.5:.1f}" y="{b - 3.5:.1f}" width="7" height="7" fill="{"white" if s.get("hollow") else col}" stroke="{col}" stroke-width="1.5"/>')
                elif marker == "^":
                    o.append(f'<polygon points="{a:.1f},{b - 4.5:.1f} {a - 4:.1f},{b + 3.5:.1f} {a + 4:.1f},{b + 3.5:.1f}" fill="{"white" if s.get("hollow") else col}" stroke="{col}" stroke-width="1.5"/>')
                elif marker == "x":
                    o.append(f'<path d="M{a - 4:.1f} {b - 4:.1f} L{a + 4:.1f} {b + 4:.1f} M{a - 4:.1f} {b + 4:.1f} L{a + 4:.1f} {b - 4:.1f}" stroke="{col}" stroke-width="2"/>')
        ly = top + 14 + 17 * k
        o.append(f'<line x1="{left + pw + 14}" y1="{ly - 4}" x2="{left + pw + 36}" y2="{ly - 4}" stroke="{col}" stroke-width="2"' + (' stroke-dasharray="6 4"' if dashed else '') + '/>')
        o.append(f'<text x="{left + pw + 42}" y="{ly}" font-size="11">{_esc(s["name"])}</text>')
    if note:
        o.append(f'<text x="{left}" y="{height - 8}" font-size="11" fill="#555">{_esc(note)}</text>')
    o.append("</svg>")
    Path(path).write_text("\n".join(o))
