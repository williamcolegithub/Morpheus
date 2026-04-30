"""PROBE-22: paired stats across folds. foundation-model vs each baseline."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from morpheus.paths import RESULTS

MODELS = ["geneformer", "geneformer_random_init", "bag_of_genes"]
PROBES = ["layer1", "layer2", "layer3_hub"]


def _load_all() -> pd.DataFrame:
    rows = []
    for f in (RESULTS / "probes").glob("*.parquet"):
        try:
            rows.append(pd.read_parquet(f))
        except Exception as exc:  # pragma: no cover
            print(f"  skip {f}: {exc}")
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def _bootstrap_ci(values: np.ndarray, n_boot: int = 5000, alpha: float = 0.05, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    if len(values) == 0:
        return (np.nan, np.nan)
    boots = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


def main() -> None:
    df = _load_all()
    if df.empty:
        print("no probe results to summarize")
        return
    print(f"loaded {len(df)} rows from results/probes/")

    # Per (probe, model[, permuted]) summary.
    if "permuted" not in df.columns:
        df["permuted"] = False
    summaries = []
    for probe in PROBES:
        sub = df[df["probe_name"] == probe]
        if sub.empty:
            continue
        for model in sub["model"].unique():
            ms = sub[sub["model"] == model]
            vals = ms["value"].to_numpy()
            lo, hi = _bootstrap_ci(vals)
            summaries.append({
                "probe": probe,
                "model": model,
                "n": int(len(vals)),
                "median": float(np.median(vals)),
                "mean": float(vals.mean()),
                "ci95_lo": lo,
                "ci95_hi": hi,
                "metric": ms["metric"].iloc[0],
            })
    summary_df = pd.DataFrame(summaries)
    print("\n=== per-probe-per-model summary ===")
    print(summary_df.to_string(index=False))

    # Paired comparisons: geneformer vs each baseline, on the SAME probe targets/folds.
    paired = []
    for probe in PROBES:
        sub = df[df["probe_name"] == probe]
        if sub.empty:
            continue
        gf = sub[sub["model"] == "geneformer"][["target", "fold", "value"]].rename(columns={"value": "v_gf"})
        for baseline in [m for m in sub["model"].unique() if m != "geneformer"]:
            bl = sub[sub["model"] == baseline][["target", "fold", "value"]].rename(columns={"value": "v_bl"})
            merged = gf.merge(bl, on=["target", "fold"])
            if len(merged) < 3:
                continue
            diffs = (merged["v_gf"] - merged["v_bl"]).to_numpy()
            try:
                stat, p = wilcoxon(diffs, alternative="two-sided")
            except ValueError:
                stat, p = (np.nan, np.nan)
            lo, hi = _bootstrap_ci(diffs)
            paired.append({
                "probe": probe,
                "baseline": baseline,
                "n_pairs": int(len(merged)),
                "median_delta": float(np.median(diffs)),
                "ci95_lo": lo,
                "ci95_hi": hi,
                "wilcoxon_stat": float(stat) if stat == stat else np.nan,
                "p_value": float(p) if p == p else np.nan,
            })
    paired_df = pd.DataFrame(paired)
    if not paired_df.empty:
        print("\n=== paired: geneformer minus baseline (Wilcoxon two-sided + bootstrap CI on the diff) ===")
        print(paired_df.to_string(index=False))

    out = RESULTS / "stats_summary.json"
    out.write_text(json.dumps({"per_model": summaries, "paired": paired}, indent=2))
    print(f"\nwrote → {out}")


if __name__ == "__main__":
    main()
