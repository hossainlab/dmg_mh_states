"""Publication output for figures: vector PDF + 300 dpi PNG, with the geometric legibility check run on every save.

Journals take the PDF (vector text, Type-42 fonts, set by dmg_atlas.figstyle) and use the PNG for review copies.
Figures are authored at the final print size - 183 mm (7.2 in) double-column - so a point is a point on the page and
nothing is rescaled later; a figure drawn at 300 mm and shrunk to 183 mm would land its 7 pt text at 4.3 pt, under
every journal's floor.
"""
from __future__ import annotations

from pathlib import Path

PAGE_W = 7.2          # 183 mm, double column
COL_W = 3.42          # 87 mm, single column
MIN_PT = 5.0          # smallest text any journal accepts at final size

# In-plot cluster labels are drawn in neutral dark, not in their cluster's colour. Colouring the text by series
# seems helpful but fights the background, because these labels sit ON their own points: "OC-like" in sky blue
# over the red OPC-like territory was unreadable however much the text was darkened. Position already ties a
# label to its cluster, so the colour channel is redundant here and only costs contrast.
LABEL_INK = "#1A1A1A"


def _label_box(col):
    """White plate behind an in-plot label, outlined in the population's own colour.

    These figures carry no legends - populations are named in plot - so the label is the ONLY place a reader can
    learn which colour means which cell type. Two earlier attempts both failed at that. Colouring the text itself
    was unreadable, because a label sits on its own points and competes with them. Making the text neutral dark
    fixed the contrast but severed the colour-to-name link entirely: "AC-like" printed over a patch holding red,
    blue and rose points tells a reader nothing about which of the three is AC-like.

    A white plate with a coloured border carries both: the plate lifts dark text clear of whatever it covers, and
    the border is an unambiguous colour key attached to the name. It also measures as part of the text's window
    extent, so the existing collision logic keeps working unchanged.
    """
    return dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=col, linewidth=1.0, alpha=0.93)


def check_figure(fig, min_pt: float = MIN_PT):
    """Geometric check: text boxes that overlap each other or a spine, text outside the canvas, text below min_pt."""
    import matplotlib as mpl
    fig.canvas.draw()  # text has no position until the figure is drawn; measuring before this reports phantom clashes
    r = fig.canvas.get_renderer()
    # matplotlib keeps Text objects for ticks outside the view limits and never draws them; they would
    # otherwise be reported as clashes and as text off the canvas
    stale = set()
    for ax in fig.axes:
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = sorted(lim)
            for tick, loc in zip(axis.get_major_ticks(), axis.get_majorticklocs()):
                if not (lo <= loc <= hi):
                    stale.add(id(tick.label1)); stale.add(id(tick.label2))
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and id(t) not in stale]
    ticklabels = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both")) for ax in fig.axes}
    # a colourbar's label legitimately sits against its own axes, so those spines are not clash candidates
    cbar_axes = {getattr(ax, "_colorbar", None) and ax for ax in fig.axes if hasattr(ax, "_colorbar")}
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes if ax not in cbar_axes
              for s in ax.spines.values() if s.get_visible()]
    all_ticks = set().union(*ticklabels.values()) if ticklabels else set()

    def rotated_tick(t):
        """Angled tick labels have overlapping boxes while their glyphs stay clear, so a box test
        over-reports on them; only upright text is compared pairwise."""
        return t in all_ticks and abs(t.get_rotation() % 180) > 1

    issues = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb) and not (rotated_tick(a) and rotated_tick(b)):
                issues.append(f"text overlap: {a.get_text()[:26]!r} / {b.get_text()[:26]!r}")
    for t, bt in texts:
        for s, bs in spines:
            if bt.overlaps(bs) and t not in ticklabels.get(s.axes, ()):
                issues.append(f"text on spine: {t.get_text()[:26]!r}")
    fb = fig.bbox
    for t, bt in texts:
        if bt.x0 < fb.x0 - 1 or bt.y0 < fb.y0 - 1 or bt.x1 > fb.x1 + 1 or bt.y1 > fb.y1 + 1:
            issues.append(f"text outside canvas: {t.get_text()[:26]!r}")
        if t.get_fontsize() < min_pt:
            issues.append(f"text below {min_pt} pt: {t.get_text()[:26]!r} at {t.get_fontsize():.1f} pt")
    w, h = fig.get_size_inches()
    if w > PAGE_W + 0.05:
        issues.append(f"figure {w:.2f} in wide ({25.4 * w:.0f} mm) exceeds the {25.4 * PAGE_W:.0f} mm page width")
    return sorted(set(issues))


def save_figure(fig, png_path, log=None, dpi: int = 300, close: bool = True):
    """Save PDF (vector, for the journal) and PNG (300 dpi, for review), reporting any legibility issues."""
    png = Path(png_path)
    issues = check_figure(fig)
    # The saved page must be exactly the authored size, or the figure does not land at 183 mm in the journal's
    # template. savefig(bbox_inches=None) does NOT do that: None falls back to rcParams["savefig.bbox"], which the
    # house style sets to "tight", and a tight box crops to the drawn content. The rcParam itself has to be cleared.
    import matplotlib as mpl
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(png, dpi=dpi)
        fig.savefig(png.with_suffix(".pdf"))
    msg = (f"saved {png.name} + {png.with_suffix('.pdf').name} "
           f"({fig.get_size_inches()[0] * 25.4:.0f} x {fig.get_size_inches()[1] * 25.4:.0f} mm)"
           + (f"; {len(issues)} legibility issue(s): " + "; ".join(issues[:6]) if issues else "; no legibility issues"))
    (log.write(msg) if log else print(msg))
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    return issues


# Nomenclature follows the report: hyphenated subpopulation and programme names, "T/NK" for the compartment.
_SPECIAL = [("T_NK_normal", "T/NK (normal)"), ("T_NK", "T/NK"), ("uncertainty_knn", "kNN uncertainty"),
            ("T_cell", "T-cell"), ("_like_normal", "-like (normal)"), ("_like", "-like"), ("_normal", " (normal)")]


def pretty(s):
    """Underscored machine labels -> reading text for axes and legends."""
    s = str(s)
    for k, v in _SPECIAL:
        s = s.replace(k, v)
    return s.replace("_", " ")


def panel(ax, s, title=None, dx=-0.06, dy=1.032, pad=3.0, letter_pt=8.0, title_pt=7.0):
    """Bold uppercase panel letter clear of the axes, with a left-aligned short title.

    The letter is right-aligned just outside the axes and the title grows rightward from the axes edge, so a long
    title cannot run into the letter however narrow the panel gets (centred titles did, once panels were
    re-authored at 183 mm).
    """
    ax.text(dx, dy, str(s).upper(), transform=ax.transAxes, fontsize=letter_pt,
            fontweight="bold", va="bottom", ha="right")
    if title:
        ax.set_title(title, loc="left", fontsize=title_pt, pad=pad)


def _readable_text(col, max_luminance: float = 0.22):
    """Darken a colour until a bold label in it is legible on white.

    In-plot labels are drawn in their cluster's own colour so the reader can tie label to points without a legend.
    A fixed 0.62 multiplier was used for this, which is fine for dark palettes and fails for light ones: 0.62 of
    viridis yellow is still yellow, and the label disappeared. Darkening now continues until the text clears a
    luminance ceiling, so the tie to the cluster colour survives for dark palettes and legibility survives for
    light ones.
    """
    import numpy as np
    import matplotlib.colors as mc
    rgb = np.array(mc.to_rgb(col))
    f = lambda c: np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    for _ in range(24):
        if float((np.array([0.2126, 0.7152, 0.0722]) * f(rgb)).sum()) <= max_luminance:
            break
        rgb = rgb * 0.88
    return tuple(rgb)


def cell_grid(ax, nrows: int, ncols: int, color: str = "0.75", lw: float = 0.5):
    """Draw cell borders on an imshow heatmap.

    A heatmap whose low end is white loses its cells against a white page - a significance matrix where most
    entries are zero stops looking like a matrix at all. Borders keep the grid readable whatever the colormap does.
    """
    import numpy as np
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", color=color, linewidth=lw)
    ax.tick_params(which="minor", length=0)
    ax.set_axisbelow(False)


def species_axis(ax, columns, fontsize: float = 6.0, group_fontsize: float = 6.0, rule: bool = True):
    """Label a heatmap x-axis whose columns are 'human:<name>' and 'mouse:<name>'.

    The prefix belongs above the axis as a group label, not on every tick as a one-letter code: 'H AC-like' and
    'M astrocyte-like (normal)' are cryptic and long, and the length is what forces these panels to reserve a third
    of their height for rotated text.
    """
    cols = [str(c) for c in columns]
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([pretty(c.split(":", 1)[1] if ":" in c else c) for c in cols],
                       rotation=90, fontsize=fontsize)
    nh = sum(c.startswith("human:") for c in cols)
    if not (0 < nh < len(cols)):
        return
    if rule:
        ax.axvline(nh - 0.5, color="k", lw=1)
    for x0, x1, name in ((0, nh, "Human"), (nh, len(cols), "Mouse")):
        ax.text((x0 + x1 - 1) / 2, 1.02, name, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=group_fontsize, style="italic")


def annotate_clusters(ax, coords, labels, palette, *, fontsize=6.0, min_label=1,
                      max_spread=0.35, background=(), names=None):
    """Name populations inside the plot instead of in a legend; returns the labels it could not place.

    Anchors on the densest bin of a 24x24 histogram rather than the median, which lands outside the population
    when it is scattered; labels are then repelled from one another and carry a white stroke so they stay
    readable over points.
    """
    import numpy as np
    import matplotlib.colors as mc
    import matplotlib.patheffects as mpe
    names = names or {}
    xy = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    xr = float(np.ptp(ax.get_xlim())) or 1.0
    yr = float(np.ptp(ax.get_ylim())) or 1.0
    items, skipped = [], []
    for lab in np.unique(labels):
        g = xy[labels == lab]
        spread = float(np.median(np.abs(g[:, 0] - np.median(g[:, 0])))) / xr
        if lab in background or len(g) < min_label or spread > max_spread:
            skipped.append(f"{names.get(lab, pretty(lab))} ({len(g)})")
            continue
        if len(g) < 400:
            x, y = float(np.median(g[:, 0])), float(np.median(g[:, 1]))
        else:
            H, xe, ye = np.histogram2d(g[:, 0], g[:, 1], bins=24)
            i, j = np.unravel_index(np.argmax(H), H.shape)
            x, y = float((xe[i] + xe[i + 1]) / 2), float((ye[j] + ye[j + 1]) / 2)
        items.append((x, y, names.get(lab, pretty(lab)), palette.get(lab, "0.35"), len(g)))
    x0, x1 = sorted(ax.get_xlim())
    y0, y1 = sorted(ax.get_ylim())
    mx, my = 0.11 * xr, 0.05 * yr  # keep the text box off the spines

    def clamp(v, lo, hi):
        return min(max(v, lo), hi)

    placed, drawn = [], []
    for x, y, lab, col, n in sorted(items, key=lambda t: -t[4]):
        dx, dy = 0.150 * xr, 0.095 * yr
        x = clamp(x, x0 + mx, x1 - mx)
        y = clamp(y, y0 + my, y1 - my)
        for _ in range(80):
            clash = [(px, py) for px, py in placed if abs(x - px) < dx and abs(y - py) < dy]
            if not clash:
                break
            y = clamp(y + (0.015 * yr) * (1 if y >= clash[0][1] else -1), y0 + my, y1 - my)
        placed.append((x, y))
        t = ax.text(x, y, lab, fontsize=fontsize, ha="center", va="center",
                    color=LABEL_INK, zorder=6, bbox=_label_box(col))
        drawn.append(t)
    _nudge_inside(ax, drawn)
    return skipped


def _nudge_inside(ax, texts, pad=2.0, rounds=12):
    """Push labels off the spines and off each other, measuring the rendered text box.

    Reserving a fixed fraction of the panel is not enough: a centred label is half its own width wide, so a long
    name for an edge cluster still crosses the spine. Measuring is exact, and needs no redraw because these
    labels carry explicit coordinates.
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    abox = ax.get_window_extent()
    for _ in range(rounds):
        moved = False
        boxes = [(t, t.get_window_extent(r)) for t in texts]
        for i, (t, b) in enumerate(boxes):
            dx = max(0.0, abox.x0 + pad - b.x0) - max(0.0, b.x1 - (abox.x1 - pad))
            dy = max(0.0, abox.y0 + pad - b.y0) - max(0.0, b.y1 - (abox.y1 - pad))
            for _t2, b2 in boxes[:i]:
                if b.overlaps(b2):
                    dy += (b2.y1 - b.y0 + pad) if b.y0 >= b2.y0 else -(b.y1 - b2.y0 + pad)
                    break
            if abs(dx) > 0.5 or abs(dy) > 0.5:
                X, Y = ax.transData.transform(t.get_position())
                t.set_position(ax.transData.inverted().transform((X + dx, Y + dy)))
                moved = True
        if not moved:
            break
