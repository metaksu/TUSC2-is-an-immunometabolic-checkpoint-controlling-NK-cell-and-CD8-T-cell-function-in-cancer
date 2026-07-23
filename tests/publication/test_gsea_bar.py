from tusc2_deg.publication import gsea_bar, pub_style
from tests.publication.panel_asserts import assert_no_clutter

def test_selected_sets_include_spine_programs():
    sel = gsea_bar.select_hallmark_sets()
    terms = " ".join(sel["nk"]["term"].tolist() + sel["cd8"]["term"].tolist()).lower()
    for prog in ("oxidative phosphorylation", "e2f", "tnf", "interferon"):
        assert prog in terms

def test_signs_match_direction_convention():
    sel = gsea_bar.select_hallmark_sets()
    cd8 = sel["cd8"].set_index("term")
    # OXPHOS higher in TUSC2-positive => negative NES
    ox = cd8[cd8.index.str.contains("oxidative phosphorylation", case=False)]
    assert (ox["NES"] < 0).all()

def test_render_returns_clean_figure():
    pub_style.apply()
    fig = gsea_bar.render()
    assert_no_clutter(fig)
