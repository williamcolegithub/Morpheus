"""PROBE-23: figures.

Fig 1: probe AUC bars + controls (per probe layer).
Fig 2: per-TF AUC distribution (Layer 2).
Fig 3: RSA values across models.
Fig 4: attention case study — average attention to GATA1 across heads/layers, on cells of erythroid lineage.
"""
from __future__ import annotations

import json
import pickle

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from morpheus.paths import CELLS_H5AD, EMBEDDINGS, GENEFORMER_DIR, RESULTS

FIG_DIR = RESULTS / "figures"


def _save(fig: plt.Figure, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  wrote → {FIG_DIR / name}.{{pdf,png}}")


def _load_all_probes() -> pd.DataFrame:
    rows = []
    for f in (RESULTS / "probes").glob("*.parquet"):
        rows.append(pd.read_parquet(f))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def fig1_probe_bars(df: pd.DataFrame) -> None:
    if df.empty:
        return
    probes = [p for p in ["layer1", "layer2", "layer3_hub"] if p in df["probe_name"].unique()]
    fig, axes = plt.subplots(1, len(probes), figsize=(5 * len(probes), 4), squeeze=False)
    for ax, probe in zip(axes[0], probes):
        sub = df[df["probe_name"] == probe]
        order = ["geneformer", "geneformer_random_init", "bag_of_genes"]
        present = [m for m in order if m in sub["model"].unique()]
        # Also include permuted variants if present.
        permuted = [m for m in sub["model"].unique() if m.endswith("__permuted")]
        groups = present + permuted
        means = [sub[sub["model"] == m]["value"].mean() for m in groups]
        sems = [sub[sub["model"] == m]["value"].sem() for m in groups]
        ax.bar(range(len(groups)), means, yerr=sems, capsize=4, color=["#1f77b4", "#888", "#ccc"][: len(groups)])
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=20, ha="right")
        ax.set_title(probe)
        metric = sub["metric"].iloc[0] if not sub.empty else ""
        ax.set_ylabel(metric)
        ax.axhline(0.5, color="r", linestyle=":", linewidth=0.8, alpha=0.5)
    fig.suptitle("Probe performance vs controls")
    fig.tight_layout()
    _save(fig, "fig1_probe_bars")


def fig2_layer2_distribution(df: pd.DataFrame) -> None:
    sub = df[df["probe_name"] == "layer2"].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    models = sorted(sub["model"].unique())
    data = [sub[sub["model"] == m]["value"].to_numpy() for m in models]
    ax.boxplot(data, labels=models, showfliers=False)
    ax.axhline(0.5, color="r", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_ylabel("AUC")
    ax.set_title("Layer 2: per-TF AUC distribution")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    _save(fig, "fig2_layer2_distribution")


def fig3_rsa() -> None:
    rsa_files = list((RESULTS / "rsa").glob("*.parquet"))
    if not rsa_files:
        return
    rows = pd.concat([pd.read_parquet(f) for f in rsa_files], ignore_index=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(len(rows)), rows["value"], color="#1f77b4")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(rows["lhs"].str.replace("gene_embeddings_", ""), rotation=20, ha="right")
    ax.set_ylabel("Spearman ρ")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_title("RSA: gene-embedding cosine vs DoRothEA TF-TF adjacency")
    fig.tight_layout()
    _save(fig, "fig3_rsa")


def fig4_attention_case_study() -> None:
    """Average attention heatmap for cells of erythroid lineage, focused on GATA1 token position."""
    npz = EMBEDDINGS / "geneformer" / "attention_sample.npz"
    if not npz.exists():
        return
    with open(GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "gene_name_id_dict_gc30M.pkl", "rb") as f:
        name2ensembl = pickle.load(f)
    with open(GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl", "rb") as f:
        ensembl2token = {k: int(v) for k, v in pickle.load(f).items()}
    gata1_ens = name2ensembl.get("GATA1")
    gata1_tok = ensembl2token.get(gata1_ens) if gata1_ens else None
    if gata1_tok is None:
        return

    npz_data = np.load(npz, allow_pickle=True)
    cell_ids = [str(x) for x in npz_data["cell_ids"]]
    attention_obj = npz_data["attention"]  # object array of (L, H, T, T) per cell, T variable
    print(f"  attention sample: {len(attention_obj)} cells, first shape {attention_obj[0].shape if len(attention_obj) else '?'}")

    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    obs = adata.obs[["cell_id", "cell_type_name"]].copy()
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs["cell_type_name"] = obs["cell_type_name"].astype(str)
    name_map = dict(zip(obs["cell_id"], obs["cell_type_name"]))
    types = [name_map.get(cid, "?").lower() for cid in cell_ids]

    erythroid_mask = np.array(["erythro" in t or "rbc" in t for t in types])
    other_mask = ~erythroid_mask
    print(f"  erythroid cells in sample: {erythroid_mask.sum()} / {len(cell_ids)}")

    # Per-layer mean attention magnitude across heads/tokens, comparing erythroid vs. other.
    n_layers = attention_obj[0].shape[0] if len(attention_obj) else 0
    def _layer_means(mask):
        if not mask.any() or n_layers == 0:
            return np.zeros(n_layers)
        means = np.zeros(n_layers)
        c = 0
        for i in np.flatnonzero(mask):
            a = attention_obj[i].astype(np.float32)  # (L, H, T, T)
            means += a.mean(axis=(1, 2, 3))
            c += 1
        return means / max(c, 1)
    layer_means_ery = _layer_means(erythroid_mask)
    layer_means_other = _layer_means(other_mask)

    fig, ax = plt.subplots(figsize=(6, 4))
    layers = np.arange(n_layers)
    width = 0.35
    ax.bar(layers - width / 2, layer_means_ery, width, label=f"erythroid (n={int(erythroid_mask.sum())})")
    ax.bar(layers + width / 2, layer_means_other, width, label=f"other (n={int(other_mask.sum())})")
    ax.set_xlabel("layer")
    ax.set_ylabel("mean attention weight")
    ax.set_title("Layer-wise attention: erythroid vs. other")
    ax.legend()
    fig.tight_layout()
    _save(fig, "fig4_attention_layer_means")


def main() -> None:
    df = _load_all_probes()
    print(f"loaded {len(df)} probe rows")
    fig1_probe_bars(df)
    fig2_layer2_distribution(df)
    fig3_rsa()
    fig4_attention_case_study()


if __name__ == "__main__":
    main()
