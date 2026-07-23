"""Shared assertions: panels must be clean (no clutter text) and vector."""
from pathlib import Path

# Forbidden on-panel substrings (these belong in the caption, never the panel).
FORBIDDEN = ["UP=", "DOWN=", "1.07x", "1.19x", "1.5x", "cap=",
             "permutation floor", "Benjamini", "padj<0.05", "n=110", "n=99"]

def collect_text(fig) -> str:
    return " ".join(t.get_text() for t in fig.findobj(match=_is_text))

def _is_text(o):
    import matplotlib.text as mtext
    return isinstance(o, mtext.Text)

def assert_no_clutter(fig):
    txt = collect_text(fig)
    bad = [s for s in FORBIDDEN if s in txt]
    assert not bad, f"panel carries caption-only clutter: {bad}"

def assert_vector_pdf(pdf_path):
    raw = Path(pdf_path).read_bytes()
    assert b"/Font" in raw, "PDF text is rasterized/outlined, not vector"
