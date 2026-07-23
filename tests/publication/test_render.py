from pathlib import Path
from tusc2_deg.publication import render
from tests.publication.panel_asserts import assert_vector_pdf

def test_render_one_panel_writes_pdf_and_png(tmp_out):
    out = render.render_panel("a", out_dir=tmp_out)
    assert Path(out + ".pdf").stat().st_size > 0
    assert Path(out + ".png").stat().st_size > 0
    assert_vector_pdf(out + ".pdf")
