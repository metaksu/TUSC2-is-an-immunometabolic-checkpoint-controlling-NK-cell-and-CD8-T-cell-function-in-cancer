"""Single entry point for the TUSC2 DEG + GSEA pipeline.

Stages (in order):
    1. pseudobulk        (data.build_pseudobulk per ctx)
    2. degs              (degs.run_deseq2 + run_wilcoxon per ctx)
    3. gsea              (gsea.run_prerank per ctx)

Alongside the stages, two audit tables are written from the live DEG frames:
the effect-size significance ladder (deg_sig_ladder.tsv) and the curated-label
provenance table (curated_label_provenance.tsv).

Caching: each stage checks its expected output files on disk. --force busts
the cache for that stage and downstream. --force-from <stage> reruns from
the named stage onward.
"""
from __future__ import annotations
import argparse
import json
import datetime
from pathlib import Path
from typing import Callable
import pandas as pd

from tusc2_deg.cd8 import config, data, degs, gsea


STAGES = ["pseudobulk", "degs", "gsea"]


def _cached_or_run(out_path: Path, runner: Callable, force: bool) -> str:
    """Return 'cached' if out_path exists and not force; else run runner and return 'ran'."""
    if out_path.exists() and not force:
        print(f"  [cached] {out_path.name}", flush=True)
        return "cached"
    runner()
    return "ran"


def _write_manifest(stage_status: dict[str, str]) -> Path:
    config.AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = config.AUDIT_DIR / "pipeline_manifest.json"
    payload = {
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "stages": stage_status,
        "config_snapshot": {
            "CONTEXTS":            config.CONTEXTS,
            "LFC_THRESH":          config.LFC_THRESH,
            "PADJ_THRESH":         config.PADJ_THRESH,
            "BASEMEAN_THRESH":     config.BASEMEAN_THRESH,
            "MIN_CELLS_PER_GROUP": config.MIN_CELLS_PER_GROUP,
            "GSEA_DIR":            str(config.GSEA_DIR),
            "LFC_NULL":            config.LFC_NULL,
        },
    }
    manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return manifest


def _write_sig_ladder(deg_dict: dict[str, pd.DataFrame], contexts: list[str]) -> Path:
    """Effect-size sensitivity ladder: UP/DOWN DEG counts at each |log2FC|
    cutoff in config.SIG_LFC_LADDER, with padj+baseMean held to the canonical
    gate. Reports counts across cutoffs (0.1/0.25/0.585) rather than a single
    permissive-gate number."""
    config.AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for ctx in contexts:
        deg = deg_dict[ctx]
        base = (deg["padj"] < config.PADJ_THRESH) & (deg["baseMean"] > config.BASEMEAN_THRESH)
        # Calibrated lfcThreshold-test count: how many genes survive the
        # calibrated gate (padj_calibrated-sig & baseMean-sig), reported so Methods
        # can state how many genes exceed the calibrated 1.1-fold null. Repeated
        # per lfc_cutoff row for the context (guarded on the column existing for
        # cached frames / fixtures that lack it).
        if "padj_calibrated" in deg.columns:
            cal = ((deg["padj_calibrated"] < config.PADJ_THRESH) &
                   (deg["baseMean"] > config.BASEMEAN_THRESH))
            n_cal = int(cal.sum())
            # ...and within the |log2FC| discovery floor — the count Methods reports
            n_cal_floored = int((cal & (deg["log2FoldChange"].abs() > config.LFC_THRESH)).sum())
        else:
            n_cal = n_cal_floored = -1  # column absent — sentinel
        for L in config.SIG_LFC_LADDER:
            up = int((base & (deg["log2FoldChange"] >  L)).sum())
            dn = int((base & (deg["log2FoldChange"] < -L)).sum())
            rows.append({"context": ctx, "lfc_cutoff": L, "UP": up, "DOWN": dn,
                         "total": up + dn, "n_calibrated_sig": n_cal,
                         "n_calibrated_sig_floored": n_cal_floored})
    out = config.AUDIT_DIR / "deg_sig_ladder.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"  [ladder] wrote {out.name}", flush=True)
    return out


def _write_curated_label_provenance(deg_dict: dict[str, pd.DataFrame],
                            contexts: list[str]) -> Path:
    """Per curated MUST_LABEL gene, emit its live LFC/padj/baseMean + a
    significance tag so the volcano labelling is auditable (the panel labels a
    literature-curated marker set, not top-N by effect size). Tags:
      [sig]              — passes the full 3-part gate (config.is_significant)
      [relaxed]          — padj-sig & baseMean-sig but |LFC| < LFC_THRESH
      [apriori-disclose] — curated but fails even the relaxed gate (or absent)
    """
    config.AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for ctx in contexts:
        must = (getattr(config, "MUST_LABEL_BY_CONTEXT", {}) or {}).get(ctx)
        if not must:
            continue
        deg = deg_dict[ctx]
        by_gene = deg.set_index("gene")
        for g in sorted(must):
            if g in by_gene.index:
                r = by_gene.loc[g]
                # guard against duplicate-gene index (take first)
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                lfc, padj, bm = (float(r["log2FoldChange"]),
                                 float(r["padj"]), float(r["baseMean"]))
                one = pd.DataFrame([{"log2FoldChange": lfc, "padj": padj,
                                     "baseMean": bm}])
                if bool(config.is_significant(one).iloc[0]):
                    tag = "[sig]"
                elif (padj < config.PADJ_THRESH) and (bm > config.BASEMEAN_THRESH):
                    tag = "[relaxed]"
                else:
                    tag = "[apriori-disclose]"
            else:
                lfc = padj = bm = float("nan")
                tag = "[apriori-disclose]"
            rows.append({"context": ctx, "gene": g, "family": config.gene_family(g),
                         "log2FoldChange": lfc, "padj": padj, "baseMean": bm, "tag": tag})
    out = config.AUDIT_DIR / "curated_label_provenance.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"  [provenance] wrote {out.name}", flush=True)
    return out


def _force_from(force_from: str | None, stage: str, force_all: bool) -> bool:
    if force_all:
        return True
    if force_from is None:
        return False
    if force_from not in STAGES:
        raise ValueError(f"--force-from {force_from} not in {STAGES}")
    return STAGES.index(stage) >= STAGES.index(force_from)


def run(force: bool = False, force_from: str | None = None,
        only_context: str | None = None) -> dict[str, str]:
    """Run all stages cache-aware. Returns stage->status dict."""
    contexts = [only_context] if only_context else config.CONTEXTS
    for ctx in contexts:
        if ctx not in config.CONTEXTS:
            raise ValueError(f"unknown context {ctx!r}")
    stage_status: dict[str, str] = {}

    # ------- 1. Load AnnDatas (always, since they're memory-only) -------
    print("[stage 1/3] load contexts", flush=True)
    adatas = data.load_all_contexts()

    # ------- 2. Pseudobulk + degs -------
    print("[stage 2/3] pseudobulk + degs", flush=True)
    f_deg = _force_from(force_from, "degs", force)
    deg_dict = {}
    for ctx in contexts:
        out = config.DEG_DIR / f"deg_pseudobulk_{ctx}.csv"
        def _runner(c=ctx, o=out):
            pb = data.build_pseudobulk(adatas[c], c)
            degs.run_deseq2(pb, c, out_path=o)
        stage_status[f"degs_pb_{ctx}"] = _cached_or_run(out, _runner, f_deg)
        deg_dict[ctx] = pd.read_csv(out)
        wout = config.DEG_DIR / f"deg_wilcoxon_{ctx}.csv"
        def _wrunner(c=ctx, o=wout):
            degs.run_wilcoxon(adatas[c], c, out_path=o)
        stage_status[f"degs_wc_{ctx}"] = _cached_or_run(wout, _wrunner, f_deg)
        # concordance audit (always re-append, no cache — additive)
        wilc = pd.read_csv(wout)
        degs.write_concordance_audit(deg_dict[ctx], wilc, ctx,
                                     audit_path=config.AUDIT_DIR / "concordance_pb_vs_wilcoxon.tsv")
    _write_sig_ladder(deg_dict, contexts)

    # ------- 3. GSEA prerank (threshold-free enrichment) -------
    print("[stage 3/3] gsea prerank", flush=True)
    f_gsea = _force_from(force_from, "gsea", force)
    gsea_dict = {}
    for ctx in contexts:
        gout = config.GSEA_DIR / f"gsea_{ctx}.csv"
        def _grunner(c=ctx, o=gout):
            gsea.run_prerank(deg_dict[c], c, out_path=o)
        stage_status[f"gsea_{ctx}"] = _cached_or_run(gout, _grunner, f_gsea)
        if gout.exists():
            gsea_dict[ctx] = pd.read_csv(gout)

    # Curated-label provenance audit — cheap, always refresh against live DEGs.
    _write_curated_label_provenance(deg_dict, contexts)

    _write_manifest(stage_status)
    return stage_status


def main():
    p = argparse.ArgumentParser(description="TUSC2 DEG + GSEA pipeline")
    p.add_argument("--force", action="store_true", help="rerun every stage")
    p.add_argument("--force-from", type=str, default=None,
                   help=f"rerun from this stage onward; one of: {STAGES}")
    p.add_argument("--context", type=str, default=None,
                   help="run only this context (else all configured contexts)")
    args = p.parse_args()
    status = run(force=args.force, force_from=args.force_from, only_context=args.context)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
