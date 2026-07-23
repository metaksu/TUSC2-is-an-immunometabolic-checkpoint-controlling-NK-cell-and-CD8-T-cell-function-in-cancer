from pathlib import Path
import matplotlib.pyplot as plt
from tusc2_deg.publication import pub_style

def test_apply_sets_vector_safe_fonttype():
    pub_style.apply()
    import matplotlib as mpl
    assert mpl.rcParams["pdf.fonttype"] == 42      # embed TrueType, not type-3
    assert mpl.rcParams["svg.fonttype"] == "none"
    assert mpl.rcParams["font.family"][0] in ("Arial", "Helvetica", "DejaVu Sans")

def test_save_writes_pdf_and_png(tmp_out):
    pub_style.apply()
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1], [0, 1])
    base = Path(tmp_out) / "probe"
    pub_style.save(fig, base)
    assert base.with_suffix(".pdf").stat().st_size > 0
    assert base.with_suffix(".png").stat().st_size > 0

def test_pdf_is_vector_text_not_rasterized(tmp_out):
    pub_style.apply()
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_xlabel("log2 fold change")
    base = Path(tmp_out) / "vec"
    pub_style.save(fig, base)
    raw = base.with_suffix(".pdf").read_bytes()
    # A vector PDF with embedded text contains font objects.
    assert b"/Font" in raw

def test_diverging_colors_present():
    assert pub_style.UP_COLOR != pub_style.DOWN_COLOR
    assert pub_style.NS_COLOR
