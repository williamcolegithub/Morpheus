"""Refreshed three-scale figures for the v2 manuscript.

Generates Figures 1-4 referenced from manuscript/sections/results.md.
Outputs land in results/figures/v2/ as both PDF and PNG.

Run:
    uv run python -m morpheus.analysis.figures_v2
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROBES = ROOT / "results" / "probes"
OUT = ROOT / "results" / "figures" / "v2"
OUT.mkdir(parents=True, exist_ok=True)

# Nature-ish palette
C_V1 = "#2E5BA8"        # blue
C_V2_104 = "#E69F00"    # orange
C_V2_316 = "#C0392B"    # red
C_BAG = "#7F7F7F"       # mid-grey
C_RND = "#BDBDBD"       # light-grey
C_PERM = "#D5D5D5"

plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def _save(fig, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)


def _load(probe: str, model: str) -> pd.DataFrame:
    df = pd.read_parquet(PROBES / f"{probe}__{model}.parquet")
    metric = "macro_f1" if probe == "layer1" else "auc"
    return df[df["metric"] == metric].copy()


def _agg_per_fold(probe: str, model: str) -> np.ndarray:
    """Mean per fold across targets (layer2/3 multi-target collapsed to one number per fold)."""
    df = _load(probe, model)
    return df.groupby("fold")["value"].mean().values


def _ci95(values: np.ndarray) -> tuple[float, float, float]:
    mean = float(np.mean(values))
    if len(values) < 2:
        return mean, mean, mean
    se = float(np.std(values, ddof=1) / np.sqrt(len(values)))
    half = 1.96 * se
    return mean, mean - half, mean + half


# ---------------------------------------------------------------------------
# Figure 1: three-scale headline (L1 / L2 / L3)
# ---------------------------------------------------------------------------
def figure_1() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))

    # Panel A: L1 macro-F1
    ax = axes[0]
    cats = [
        ("V1", "geneformer", C_V1),
        ("V2-104M", "Geneformer-V2-104M", C_V2_104),
        ("V2-316M", "Geneformer-V2-316M", C_V2_316),
        ("bag-of-genes", "bag_of_genes", C_BAG),
        ("V1 rand", "geneformer_random_init", C_RND),
        ("V2-104M rand", "Geneformer-V2-104M_random_init", C_RND),
        ("V2-316M rand", "Geneformer-V2-316M_random_init", C_RND),
    ]
    means, los, his, colors, labels = [], [], [], [], []
    for label, model, color in cats:
        vals = _agg_per_fold("layer1", model)
        m, lo, hi = _ci95(vals)
        means.append(m); los.append(m - lo); his.append(hi - m)
        colors.append(color); labels.append(label)
    xs = np.arange(len(labels))
    ax.bar(xs, means, yerr=[los, his], color=colors, edgecolor="black", linewidth=0.6, capsize=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("L1 cell-type macro-F1")
    ax.set_title("A  Layer 1: cell-type identity")
    ax.set_ylim(0, 0.6)
    ax.axhline(0.404, ls="--", color=C_BAG, lw=0.7, zorder=0)

    # Panel B: L2 AUC distribution per scale (boxplot)
    ax = axes[1]
    box_models = [
        ("V1", "geneformer", C_V1),
        ("V2-104M", "Geneformer-V2-104M", C_V2_104),
        ("V2-316M", "Geneformer-V2-316M", C_V2_316),
        ("bag", "bag_of_genes", C_BAG),
        ("permuted", "geneformer__permuted", C_PERM),
    ]
    data, labels, colors = [], [], []
    for label, model, color in box_models:
        df = _load("layer2", model)
        per_target = df.groupby("target")["value"].mean()
        data.append(per_target.values)
        labels.append(label)
        colors.append(color)
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("L2 TF→target AUC (per-TF mean)")
    ax.set_title("B  Layer 2: TF→target edges")
    ax.axhline(0.5, ls=":", color="black", lw=0.6)

    # Panel C: L3 hub-AUC
    ax = axes[2]
    cats = [
        ("V1", "geneformer", C_V1),
        ("V2-104M", "Geneformer-V2-104M", C_V2_104),
        ("V2-316M", "Geneformer-V2-316M", C_V2_316),
        ("bag-of-genes", "bag_of_genes", C_BAG),
        ("V1 rand", "geneformer_random_init", C_RND),
        ("V2-104M rand", "Geneformer-V2-104M_random_init", C_RND),
        ("V2-316M rand", "Geneformer-V2-316M_random_init", C_RND),
    ]
    means, los, his, colors, labels = [], [], [], [], []
    for label, model, color in cats:
        vals = _agg_per_fold("layer3", model)
        m, lo, hi = _ci95(vals)
        means.append(m); los.append(m - lo); his.append(hi - m)
        colors.append(color); labels.append(label)
    xs = np.arange(len(labels))
    ax.bar(xs, means, yerr=[los, his], color=colors, edgecolor="black", linewidth=0.6, capsize=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("L3 hub-TF AUC")
    ax.set_title("C  Layer 3: hub identity")
    ax.set_ylim(0, 0.85)
    ax.axhline(0.5, ls=":", color="black", lw=0.6)

    fig.tight_layout()
    _save(fig, "fig1")


# ---------------------------------------------------------------------------
# Figure 2: forest plot of paired contrasts
# ---------------------------------------------------------------------------
def _paired_delta(probe: str, m_a: str, m_b: str, key: str = "target", agg: str = "target") -> tuple[float, float, float, int]:
    """Compute mean Δ and 95% CI for paired comparison m_a - m_b.

    For L1/L3 (target is a singleton: cell_type / is_hub) we pair on fold (n=5).
    For L2 we average folds within target then bootstrap across targets.
    """
    a = _load(probe, m_a)
    b = _load(probe, m_b)
    merged = a.merge(b, on=["target", "fold"], suffixes=("_a", "_b"))
    merged["delta"] = merged["value_a"] - merged["value_b"]
    if agg == "fold":
        per = merged.groupby("fold")["delta"].mean().values
    else:
        per = merged.groupby("target")["delta"].mean().values
    mean = float(np.mean(per))
    rng = np.random.default_rng(0)
    boots = [rng.choice(per, size=len(per), replace=True).mean() for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return mean, float(lo), float(hi), len(per)


def figure_2() -> None:
    rows = [
        ("L1  V2-104M − V1",          "layer1", "Geneformer-V2-104M", "geneformer",         "fold"),
        ("L1  V2-316M − V1",          "layer1", "Geneformer-V2-316M", "geneformer",         "fold"),
        ("L1  V2-316M − V2-104M",     "layer1", "Geneformer-V2-316M", "Geneformer-V2-104M", "fold"),
        ("L1  V2-316M − bag-of-genes","layer1", "Geneformer-V2-316M", "bag_of_genes",       "fold"),
        ("L2  V2-104M − V1",          "layer2", "Geneformer-V2-104M", "geneformer",         "target"),
        ("L2  V2-316M − V1",          "layer2", "Geneformer-V2-316M", "geneformer",         "target"),
        ("L2  V2-316M − V2-104M",     "layer2", "Geneformer-V2-316M", "Geneformer-V2-104M", "target"),
        ("L2  V2-316M − bag-of-genes","layer2", "Geneformer-V2-316M", "bag_of_genes",       "target"),
        ("L3  V2-104M − V1",          "layer3", "Geneformer-V2-104M", "geneformer",         "fold"),
        ("L3  V2-316M − V1",          "layer3", "Geneformer-V2-316M", "geneformer",         "fold"),
        ("L3  V2-316M − V2-104M",     "layer3", "Geneformer-V2-316M", "Geneformer-V2-104M", "fold"),
        ("L3  V2-316M − bag-of-genes","layer3", "Geneformer-V2-316M", "bag_of_genes",       "fold"),
    ]

    means, los, his, labels, ns = [], [], [], [], []
    for label, probe, a, b, agg in rows:
        m, lo, hi, n = _paired_delta(probe, a, b, agg=agg)
        means.append(m); los.append(lo); his.append(hi); labels.append(label); ns.append(n)

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    y = np.arange(len(labels))[::-1]
    for yi, m, lo, hi, lbl in zip(y, means, los, his, labels):
        if "L1" in lbl:
            color = C_V1 if "V1" in lbl.split("−")[1] else C_V2_104
        if lbl.startswith("L1"): color = "#1f4e8a"
        elif lbl.startswith("L2"): color = "#a85b00"
        else: color = "#882e22"
        ax.errorbar(m, yi, xerr=[[m - lo], [hi - m]], fmt="o", color=color, capsize=3, lw=1.4, ms=5)
    ax.axvline(0.0, color="black", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}  (n={n})" for l, n in zip(labels, ns)])
    ax.set_xlabel("paired Δ (mean, 95% bootstrap CI)")
    ax.set_title("Paired contrasts: within-V2 vs cross-version")
    fig.tight_layout()
    _save(fig, "fig2")


# ---------------------------------------------------------------------------
# Figure 3: random-init geometry control (L3 hub)
# ---------------------------------------------------------------------------
def figure_3() -> None:
    scales = [("V1\n(256d)", "geneformer", "geneformer_random_init", C_V1),
              ("V2-104M\n(768d)", "Geneformer-V2-104M", "Geneformer-V2-104M_random_init", C_V2_104),
              ("V2-316M\n(1152d)", "Geneformer-V2-316M", "Geneformer-V2-316M_random_init", C_V2_316)]

    trained_means, trained_err = [], []
    rand_means, rand_err = [], []
    labels, colors = [], []
    for label, mt, mr, color in scales:
        tv = _agg_per_fold("layer3", mt); rv = _agg_per_fold("layer3", mr)
        m1, lo1, hi1 = _ci95(tv); m2, lo2, hi2 = _ci95(rv)
        trained_means.append(m1); trained_err.append([(m1 - lo1), (hi1 - m1)])
        rand_means.append(m2); rand_err.append([(m2 - lo2), (hi2 - m2)])
        labels.append(label); colors.append(color)

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    x = np.arange(len(labels))
    w = 0.35
    t_err = np.array(trained_err).T
    r_err = np.array(rand_err).T
    ax.bar(x - w / 2, trained_means, w, yerr=t_err, color=colors, edgecolor="black", lw=0.6, capsize=3, label="trained")
    ax.bar(x + w / 2, rand_means, w, yerr=r_err, color=C_RND, edgecolor="black", lw=0.6, capsize=3, label="random init")
    ax.axhline(0.5, ls=":", color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("L3 hub-TF AUC")
    ax.set_title("Trained vs random-init at matched dimensionality")
    ax.set_ylim(0, 0.85)
    ax.legend(frameon=False, loc="upper right")

    # delta annotations
    for xi, t, r in zip(x, trained_means, rand_means):
        ax.annotate(f"Δ={t - r:+.2f}", xy=(xi, max(t, r) + 0.04), ha="center", fontsize=8)

    fig.tight_layout()
    _save(fig, "fig3")


# ---------------------------------------------------------------------------
# Figure 4: layer-wise V2-104M L1 macro-F1
# ---------------------------------------------------------------------------
def figure_4() -> None:
    df = pd.read_parquet(PROBES / "layer1_per_layer__Geneformer-V2-104M_layerwise.parquet")
    df = df[df.metric == "macro_f1"]
    g = df.groupby("layer")["value"].agg(["mean", "std", "count"]).reset_index()
    g["se"] = g["std"] / np.sqrt(g["count"])

    bag_vals = _agg_per_fold("layer1", "bag_of_genes")
    bag_mean = float(np.mean(bag_vals))
    final_mean = _agg_per_fold("layer1", "Geneformer-V2-104M").mean()

    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.errorbar(g["layer"], g["mean"], yerr=g["se"], fmt="o-", color=C_V2_104, lw=1.6, ms=5, capsize=3, label="V2-104M layer-wise")
    ax.axhline(bag_mean, ls="--", color=C_BAG, lw=1.0, label=f"bag-of-genes ({bag_mean:.3f})")
    ax.axhline(final_mean, ls=":", color=C_V2_104, lw=1.0, label=f"V2-104M final extraction ({final_mean:.3f})")
    ax.set_xlabel("Transformer layer (0 = embedding, 12 = final)")
    ax.set_ylabel("L1 cell-type macro-F1")
    ax.set_title("Layer-wise probe of Geneformer V2-104M (L1)")
    ax.set_xticks(range(0, 13))
    ax.legend(frameon=False, loc="lower right")
    ax.set_ylim(0.30, 0.55)
    fig.tight_layout()
    _save(fig, "fig4")


def main() -> None:
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    print(f"wrote figures to {OUT}")


if __name__ == "__main__":
    main()
