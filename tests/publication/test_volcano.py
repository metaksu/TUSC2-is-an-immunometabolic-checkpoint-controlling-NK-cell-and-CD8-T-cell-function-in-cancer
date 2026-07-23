import itertools

from tusc2_deg.publication import volcano, pub_style
from tests.publication.panel_asserts import assert_no_clutter


def test_cd8_volcano_labels_include_hspa6_pdcd1_gzmk():
    # The rendered label set = curated provenance + figure anchors. HSPA6 is
    # absent from the curated provenance but is the lone CD8 gene past 1.5x and is
    # called out in the text, so it is force-labelled via volcano.ANCHORS.
    labels = volcano.label_genes("cd8")
    assert "HSPA6" in labels
    assert "PDCD1" in labels
    assert "GZMK" in labels
    # and it must be a genuinely significant point, not a fabricated label
    from tusc2_deg.publication import data
    cd8 = data.load_pseudobulk("cd8").set_index("gene")
    assert cd8.loc["HSPA6", "padj"] < 0.05 and abs(cd8.loc["HSPA6", "log2FoldChange"]) > 0.585


def test_cd8_render_clean():
    pub_style.apply()
    fig = volcano.render("cd8")
    assert_no_clutter(fig)


def test_nk_render_clean():
    pub_style.apply()
    fig = volcano.render("nk")
    assert_no_clutter(fig)


def test_axis_label_states_direction_not_formula():
    pub_style.apply()
    fig = volcano.render("cd8")
    xlabels = [a.get_xlabel() for a in fig.axes]
    assert any("log" in x.lower() for x in xlabels)
    assert all("1.07x" not in x for x in xlabels)   # no FC-decode formula


def _label_boxes(fig):
    """Display-space bounding boxes of every curated gene label in the panel."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax = fig.axes[0]
    return [(t.get_text(), t.get_window_extent(renderer=renderer)) for t in ax.texts]


def _min_horizontal_gap(fig):
    """Smallest horizontal gap (px) between any two y-overlapping labels.

    Negative means the two label boxes overlap in x; ~0 means they abut and
    read as a single run-on token.
    """
    boxes = _label_boxes(fig)
    worst = float("inf")
    for (ta, ba), (tb, bb) in itertools.combinations(boxes, 2):
        if ba.y1 < bb.y0 or bb.y1 < ba.y0:      # no vertical overlap -> can't run on
            continue
        if ba.x1 <= bb.x0:
            gap = bb.x0 - ba.x1
        elif bb.x1 <= ba.x0:
            gap = ba.x0 - bb.x1
        else:
            gap = -1.0                          # overlapping in x
        worst = min(worst, gap)
    return worst


def test_cd8_curated_labels_do_not_abut():
    # Regression: GZMK and CD69 (and PRF1/IFNG) must not collide into a run-on
    # token. Require a clear horizontal gap between any two y-overlapping labels.
    pub_style.apply()
    fig = volcano.render("cd8")
    assert _min_horizontal_gap(fig) >= 4.0, "CD8 curated labels abut (run-on token)"


def test_nk_signal_is_small_zero_past_1p5x():
    from tusc2_deg.publication import data
    nk = data.load_pseudobulk("nk")
    past = (nk["padj"] < 0.05) & (nk["baseMean"] > 5) & (nk["log2FoldChange"].abs() > 0.585)
    assert past.sum() == 0          # NK: no gene past 1.5-fold (the 'selective' message)
