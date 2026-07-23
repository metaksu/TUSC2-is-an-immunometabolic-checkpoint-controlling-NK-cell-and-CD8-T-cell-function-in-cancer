import math
import numpy as np
from tusc2_deg.publication import running_es, pub_style
from tests.publication.panel_asserts import assert_no_clutter

def test_running_es_min_is_negative_for_oxphos():
    res = running_es.compute("cd8", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    # OXPHOS enriched at the bottom (TUSC2-positive) => most-extreme ES negative
    assert res["es_curve"].min() < 0
    assert abs(res["es_curve"].min()) >= abs(res["es_curve"].max())

def test_running_es_curve_returns_to_zero():
    res = running_es.compute("cd8", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    assert abs(res["es_curve"][-1]) < 1e-6

def test_render_clean():
    pub_style.apply()
    fig = running_es.render("cd8", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    assert_no_clutter(fig)

def test_nk_oxphos_nes_is_hallmark_not_go_row():
    # NK OXPHOS must annotate the Hallmark NES (-1.94), not the higher-|NES| GO
    # row (-1.97) — so the panel matches the GSEA bar.
    res = running_es.compute("nk", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    assert math.isclose(res["nes"], -1.94, abs_tol=0.02)

def test_paired_render_clean_and_shares_es_axis():
    pub_style.apply()
    fig = running_es.render_pair("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    assert_no_clutter(fig)
    # both ES sub-axes must share the same y-limits so the depth contrast is honest
    es_axes = [ax for ax in fig.findobj(lambda o: hasattr(o, "get_ylabel"))
               if getattr(ax, "get_ylabel", lambda: "")() == "Enrichment score (ES)"]
    assert len(es_axes) == 2
    assert es_axes[0].get_ylim() == es_axes[1].get_ylim()
