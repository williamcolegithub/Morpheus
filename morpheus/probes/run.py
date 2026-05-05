"""Probe runner CLI: PROBE-13 / PROBE-14 / PROBE-15 + control variants.

Usage:
  python -m morpheus.probes.run --probe layer1 --model geneformer
  python -m morpheus.probes.run --probe layer2 --model geneformer
  python -m morpheus.probes.run --probe layer3 --model geneformer
  python -m morpheus.probes.run --probe layer2 --model geneformer --permute-labels   # PROBE-19

Models: geneformer | geneformer_random_init | bag_of_genes

Reads embeddings from data/embeddings/<model>/, splits from data/processed/splits.parquet,
labels from cells.h5ad (Layer 1) or dorothea_edges.parquet (Layer 2/3).
Writes per-fold rows to results/probes/<probe>__<model>[__permuted].parquet.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from tqdm import tqdm

from morpheus.paths import CELLS_H5AD, DOROTHEA_EDGES, EMBEDDINGS, RESULTS, SPLITS

CELL_TYPE_MIN_CELLS = 200  # drop ultra-rare types from Layer 1; they wreck CV folds
TF_MIN_TARGETS = 30
HUB_TOP_DECILE = 0.10
MLP_HIDDEN_LAYER = (64,)


def _make_classifier(probe_type: str, multiclass: bool, seed: int):
    if probe_type == "mlp":
        return MLPClassifier(
            hidden_layer_sizes=MLP_HIDDEN_LAYER,
            max_iter=300,
            early_stopping=True,
            random_state=seed,
        )
    # linear, default
    if multiclass:
        return LogisticRegression(max_iter=1000, n_jobs=-1, C=1.0, random_state=seed)
    return LogisticRegression(max_iter=2000, C=1.0, random_state=seed)


def _binary_score(clf, X):
    """Return real-valued scores for binary AUC, working for both LogisticRegression and MLPClassifier."""
    if hasattr(clf, "decision_function"):
        try:
            return clf.decision_function(X)
        except (AttributeError, NotImplementedError):
            pass
    proba = clf.predict_proba(X)
    return proba[:, 1] if proba.ndim == 2 and proba.shape[1] == 2 else proba


def _load_cell_embeddings(model: str) -> tuple[np.ndarray | "csr_matrix", pd.DataFrame, dict]:
    d = EMBEDDINGS / model
    meta = json.loads((d / "META.json").read_text())
    cell_idx = pd.read_parquet(d / "cell_index.parquet")
    if (d / "cells.npz").exists():
        X = load_npz(d / "cells.npz")
    else:
        X = np.load(d / "cells.npy", mmap_mode="r")
    return X, cell_idx, meta


def _load_gene_embeddings(model: str) -> tuple[np.ndarray, pd.DataFrame, dict]:
    d = EMBEDDINGS / model
    meta = json.loads((d / "META.json").read_text())
    gene_idx = pd.read_parquet(d / "gene_index.parquet")
    G = np.load(d / "genes.npy")
    return G, gene_idx, meta


def _result_path(probe: str, model: str, permuted: bool) -> Path:
    suffix = "__permuted" if permuted else ""
    return RESULTS / "probes" / f"{probe}__{model}{suffix}.parquet"


def layer1(model: str, seed: int, probe_type: str = "linear") -> pd.DataFrame:
    """Cell-type identity. Macro-F1 per fold, donor-stratified."""
    X, cell_idx, meta = _load_cell_embeddings(model)
    splits = pd.read_parquet(SPLITS).merge(cell_idx, on="cell_id")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    obs = adata.obs[["cell_id", "cell_type"]].copy()
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    df = splits.merge(obs, on="cell_id").sort_values("row").reset_index(drop=True)
    counts = df["cell_type"].value_counts()
    keep_types = counts[counts >= CELL_TYPE_MIN_CELLS].index
    df = df[df["cell_type"].isin(keep_types)].reset_index(drop=True)
    print(f"  layer1 ({probe_type}): kept {len(df):,} cells across {df['cell_type'].nunique()} types (min={CELL_TYPE_MIN_CELLS})")

    rows = []
    for fold in sorted(df["fold"].unique()):
        train = df[df["fold"] != fold]
        test = df[df["fold"] == fold]
        Xtr = X[train["row"].to_numpy()]
        Xte = X[test["row"].to_numpy()]
        if hasattr(Xtr, "toarray"):
            Xtr = Xtr.toarray()
            Xte = Xte.toarray()
        scaler = StandardScaler(with_mean=False)
        Xtr = scaler.fit_transform(Xtr)
        Xte = scaler.transform(Xte)
        clf = _make_classifier(probe_type, multiclass=True, seed=seed)
        # MLPClassifier in sklearn 1.8 chokes on string class labels under early_stopping
        # (calls np.isnan on string predictions). Encode to ints up front; decode after predict.
        le = LabelEncoder()
        ytr = le.fit_transform(train["cell_type"].values)
        yte_true = test["cell_type"].values
        clf.fit(Xtr, ytr)
        pred_int = clf.predict(Xte)
        pred = le.inverse_transform(pred_int)
        f1 = f1_score(yte_true, pred, average="macro", zero_division=0)
        rows.append({"probe_name": "layer1", "target": "cell_type", "model": model, "fold": int(fold), "metric": "macro_f1", "value": float(f1), "n_train": int(len(train)), "n_test": int(len(test)), "seed": seed})
        print(f"    fold {fold}: macro_f1={f1:.4f}")
    return pd.DataFrame(rows)


def _gene_row_lookup(gene_idx: pd.DataFrame) -> dict[str, int]:
    return dict(zip(gene_idx["ensembl_id"].astype(str), gene_idx["row"].astype(int)))


def _expression_bins(adata: ad.AnnData, gene_idx: pd.DataFrame, n_bins: int = 10) -> dict[str, int]:
    """Bin genes by mean log-expression for matched-negative sampling. Cheap streaming computation."""
    var_to_idx = {e: i for i, e in enumerate(adata.var_names.values)}
    means: dict[str, float] = {}
    # Compute mean expression per gene across all cells, in chunks.
    n_cells = adata.n_obs
    sums = np.zeros(adata.n_vars, dtype=np.float64)
    chunk = 5000
    for s in range(0, n_cells, chunk):
        e = min(s + chunk, n_cells)
        X = adata.X[s:e]
        if hasattr(X, "toarray"):
            X = X.toarray()
        sums += np.asarray(X, dtype=np.float64).sum(axis=0)
    mean_per_gene = sums / max(n_cells, 1)
    # Map ensembl → mean.
    for ensembl in gene_idx["ensembl_id"].astype(str):
        i = var_to_idx.get(ensembl)
        if i is not None:
            means[ensembl] = float(mean_per_gene[i])
    # Bin among gene_idx genes.
    arr = np.array([means.get(e, 0.0) for e in gene_idx["ensembl_id"].astype(str)])
    quantiles = np.quantile(arr, np.linspace(0, 1, n_bins + 1))
    bins = np.clip(np.searchsorted(quantiles, arr, side="right") - 1, 0, n_bins - 1)
    return dict(zip(gene_idx["ensembl_id"].astype(str), bins.tolist()))


def layer2(
    model: str,
    seed: int,
    permute: bool,
    probe_type: str = "linear",
    neg_mode: str = "matched",
    regulon_source: str = "all",
) -> pd.DataFrame:
    """TF→target probe per TF. AUC distribution across TFs.

    Sensitivity-analysis flags (audit hardening):
      - probe_type: linear (LogReg) | mlp (1-hidden 64-unit MLP)
      - neg_mode: matched (expression-bin matched) | random (uniform from non-targets)
      - regulon_source: all (DoRothEA ∪ CollecTRI) | dorothea | collectri
    """
    G, gene_idx, _ = _load_gene_embeddings(model)
    g2row = _gene_row_lookup(gene_idx)
    edges = pd.read_parquet(DOROTHEA_EDGES)
    edges = edges[edges["tf_ensembl"].isin(g2row) & edges["target_ensembl"].isin(g2row)].copy()
    if regulon_source != "all":
        if "source" not in edges.columns:
            raise RuntimeError("dorothea_edges.parquet has no 'source' column; cannot filter by regulon_source")
        edges = edges[edges["source"] == regulon_source].copy()
        if edges.empty:
            raise RuntimeError(f"no edges with source={regulon_source}")

    if neg_mode == "matched":
        print(f"  layer2 ({probe_type}, neg={neg_mode}, regulon={regulon_source}): computing gene-expression bins for matched negatives ...")
        adata = ad.read_h5ad(CELLS_H5AD, backed="r")
        bins_map = _expression_bins(adata, gene_idx)
    else:
        print(f"  layer2 ({probe_type}, neg={neg_mode}, regulon={regulon_source}): random-pair negatives")
        bins_map = None

    rng = np.random.default_rng(seed)
    rows = []
    tf_groups = edges.groupby("tf_ensembl")["target_ensembl"].apply(lambda s: set(s)).to_dict()

    candidate_genes = gene_idx["ensembl_id"].astype(str).to_numpy()
    candidate_bins = np.array([bins_map[g] for g in candidate_genes]) if bins_map is not None else None

    for tf, targets in tqdm(tf_groups.items(), desc="layer2"):
        if len(targets) < TF_MIN_TARGETS:
            continue
        pos_genes = list(targets)
        if permute:
            mask = candidate_genes != tf
            pos_genes = list(rng.choice(candidate_genes[mask], size=len(pos_genes), replace=False))

        pos_set = set(pos_genes) | {tf}
        if neg_mode == "matched":
            neg_pool: list[str] = []
            target_bin_counts = pd.Series([bins_map.get(g, 0) for g in pos_genes]).value_counts().to_dict()
            for b, n in target_bin_counts.items():
                mask = (candidate_bins == b) & np.array([g not in pos_set for g in candidate_genes])
                pool = candidate_genes[mask]
                if len(pool) == 0:
                    continue
                chosen = rng.choice(pool, size=min(n, len(pool)), replace=False)
                neg_pool.extend(chosen.tolist())
        else:  # random uniform negatives
            mask = np.array([g not in pos_set for g in candidate_genes])
            pool = candidate_genes[mask]
            n_neg = min(len(pos_genes), len(pool))
            neg_pool = rng.choice(pool, size=n_neg, replace=False).tolist()
        if not neg_pool:
            continue
        all_genes = pos_genes + neg_pool
        y = np.array([1] * len(pos_genes) + [0] * len(neg_pool))
        Xrows = np.array([g2row[g] for g in all_genes])
        Xg = G[Xrows]

        n = len(y)
        idx = rng.permutation(n)
        split = int(0.8 * n)
        tr, te = idx[:split], idx[split:]
        if y[te].sum() == 0 or y[te].sum() == len(te):
            continue
        clf = _make_classifier(probe_type, multiclass=False, seed=seed)
        clf.fit(Xg[tr], y[tr])
        scores = _binary_score(clf, Xg[te])
        auc = float(roc_auc_score(y[te], scores))
        rows.append({"probe_name": "layer2", "target": tf, "model": model, "fold": 0, "metric": "auc", "value": auc, "n_train": int(len(tr)), "n_test": int(len(te)), "seed": seed})
    if rows:
        df = pd.DataFrame(rows)
        print(f"  layer2: {len(df)} TFs probed; AUC median={df['value'].median():.3f}, IQR=[{df['value'].quantile(0.25):.3f}, {df['value'].quantile(0.75):.3f}]")
        return df
    return pd.DataFrame(rows)


def layer3(model: str, seed: int, probe_type: str = "linear", regulon_source: str = "all") -> pd.DataFrame:
    """Hub identity probe."""
    G, gene_idx, _ = _load_gene_embeddings(model)
    g2row = _gene_row_lookup(gene_idx)
    edges = pd.read_parquet(DOROTHEA_EDGES)
    edges = edges[edges["tf_ensembl"].isin(g2row) & edges["target_ensembl"].isin(g2row)].copy()
    if regulon_source != "all":
        edges = edges[edges["source"] == regulon_source].copy()

    out_degree = edges.groupby("tf_ensembl").size()
    threshold = out_degree.quantile(1 - HUB_TOP_DECILE)
    tfs = out_degree.index.to_numpy()
    is_hub = (out_degree >= threshold).astype(int).to_numpy()
    print(f"  layer3 ({probe_type}, regulon={regulon_source}): {len(tfs)} TFs; {is_hub.sum()} hubs (out-degree ≥ {threshold:.0f}, top {int(HUB_TOP_DECILE*100)}%)")

    Xrows = np.array([g2row[t] for t in tfs])
    Xg = G[Xrows]

    rng = np.random.default_rng(seed)
    rows = []
    K = 5
    fold_idx = rng.integers(0, K, size=len(tfs))
    for k in range(K):
        tr = fold_idx != k
        te = fold_idx == k
        if is_hub[te].sum() == 0 or is_hub[te].sum() == te.sum():
            continue
        clf = _make_classifier(probe_type, multiclass=False, seed=seed)
        clf.fit(Xg[tr], is_hub[tr])
        scores = _binary_score(clf, Xg[te])
        auc = float(roc_auc_score(is_hub[te], scores))
        rows.append({"probe_name": "layer3_hub", "target": "is_hub", "model": model, "fold": int(k), "metric": "auc", "value": auc, "n_train": int(tr.sum()), "n_test": int(te.sum()), "seed": seed})
        print(f"    fold {k}: hub_auc={auc:.4f}")

    return pd.DataFrame(rows)


def rsa(model: str, seed: int) -> pd.DataFrame:
    """RSA: Spearman between gene-embedding cosine similarity and DoRothEA adjacency, restricted to TFs."""
    G, gene_idx, _ = _load_gene_embeddings(model)
    g2row = _gene_row_lookup(gene_idx)
    edges = pd.read_parquet(DOROTHEA_EDGES)
    edges = edges[edges["tf_ensembl"].isin(g2row) & edges["target_ensembl"].isin(g2row)].copy()
    tfs = sorted(edges["tf_ensembl"].unique())
    if len(tfs) < 10:
        return pd.DataFrame()
    rows_idx = np.array([g2row[t] for t in tfs])
    M = G[rows_idx]
    Mn = M / np.clip(np.linalg.norm(M, axis=1, keepdims=True), 1e-9, None)
    cos = Mn @ Mn.T

    # Adjacency on TFs (binary directed: TF i regulates TF j).
    tf_to_i = {t: i for i, t in enumerate(tfs)}
    A = np.zeros((len(tfs), len(tfs)), dtype=np.int8)
    for _, row in edges.iterrows():
        i = tf_to_i.get(row["tf_ensembl"])
        j = tf_to_i.get(row["target_ensembl"])
        if i is not None and j is not None:
            A[i, j] = 1

    iu = np.triu_indices(len(tfs), k=1)
    sim = cos[iu]
    adj_sym = (A | A.T)[iu]  # symmetric: regulates-or-regulated-by
    rho, _ = spearmanr(sim, adj_sym)
    print(f"  rsa: n_tfs={len(tfs)}, spearman_rho={rho:.4f}")
    return pd.DataFrame([{
        "lhs": f"gene_embeddings_{model}",
        "rhs": "dorothea_tf_tf_adjacency_symmetric",
        "correlation_metric": "spearman",
        "value": float(rho),
        "n": int(len(sim)),
        "seed": seed,
    }])


def _variant_tag(probe_type: str, neg_mode: str, regulon_source: str, permute: bool) -> str:
    parts: list[str] = []
    if probe_type != "linear":
        parts.append(probe_type)
    if neg_mode != "matched":
        parts.append(f"neg-{neg_mode}")
    if regulon_source != "all":
        parts.append(f"reg-{regulon_source}")
    if permute:
        parts.append("permuted")
    return "__" + "__".join(parts) if parts else ""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--probe", required=True, choices=["layer1", "layer2", "layer3", "rsa", "all"])
    p.add_argument("--model", required=True)
    p.add_argument("--permute-labels", action="store_true")
    p.add_argument("--probe-type", default="linear", choices=["linear", "mlp"], help="probe model class (default linear LR; mlp = 1-hidden-layer 64-unit MLP)")
    p.add_argument("--neg-mode", default="matched", choices=["matched", "random"], help="layer2 negative sampling: matched (expression-bin) or random uniform")
    p.add_argument("--regulon-source", default="all", choices=["all", "dorothea", "collectri"], help="which regulon subset to use for layer2/layer3")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    tag = _variant_tag(args.probe_type, args.neg_mode, args.regulon_source, args.permute_labels)

    probes = ["layer1", "layer2", "layer3", "rsa"] if args.probe == "all" else [args.probe]
    for probe in probes:
        if probe == "layer1":
            df = layer1(args.model, args.seed, probe_type=args.probe_type)
        elif probe == "layer2":
            df = layer2(args.model, args.seed, args.permute_labels, probe_type=args.probe_type, neg_mode=args.neg_mode, regulon_source=args.regulon_source)
        elif probe == "layer3":
            df = layer3(args.model, args.seed, probe_type=args.probe_type, regulon_source=args.regulon_source)
        elif probe == "rsa":
            df = rsa(args.model, args.seed)
        else:
            raise ValueError(probe)

        if probe == "rsa":
            out = RESULTS / "rsa" / f"rsa__{args.model}{tag}.parquet"
        else:
            out = RESULTS / "probes" / f"{probe}__{args.model}{tag}.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        if tag and len(df) > 0 and "model" in df.columns:
            df = df.copy()
            df["model"] = df["model"].astype(str) + tag
        if len(df) > 0:
            df.to_parquet(out, index=False)
            print(f"wrote → {out} ({len(df)} rows)")
        else:
            print(f"  {probe}: no rows produced")


if __name__ == "__main__":
    main()
