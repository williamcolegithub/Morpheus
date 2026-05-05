"""Layer-wise Layer-1 probe (cell-type identity) on a model with cells_layers.npy saved.

For each transformer layer (including embedding layer 0), runs the standard donor-stratified
5-fold logistic-regression probe on the mean-pooled hidden state at that layer, and reports
macro-F1 per layer.

Usage:
  python -m morpheus.probes.layerwise --model Geneformer-V2-104M_layerwise
"""
from __future__ import annotations

import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler

from morpheus.paths import CELLS_H5AD, EMBEDDINGS, RESULTS, SPLITS

CELL_TYPE_MIN_CELLS = 200


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    d = EMBEDDINGS / args.model
    cell_idx = pd.read_parquet(d / "cell_index.parquet")
    layers = np.load(d / "cells_layers.npy", mmap_mode="r")
    n_cells, n_layers, embed_dim = layers.shape
    print(f"loaded {args.model}/cells_layers.npy: {layers.shape}")

    splits = pd.read_parquet(SPLITS).merge(cell_idx, on="cell_id")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    obs = adata.obs[["cell_id", "cell_type"]].copy()
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    df = splits.merge(obs, on="cell_id").sort_values("row").reset_index(drop=True)
    counts = df["cell_type"].value_counts()
    keep = counts[counts >= CELL_TYPE_MIN_CELLS].index
    df = df[df["cell_type"].isin(keep)].reset_index(drop=True)
    print(f"  kept {len(df):,} cells, {df['cell_type'].nunique()} cell types")

    rows = []
    for layer in range(n_layers):
        for fold in sorted(df["fold"].unique()):
            train = df[df["fold"] != fold]
            test = df[df["fold"] == fold]
            Xtr = layers[train["row"].to_numpy(), layer, :]
            Xte = layers[test["row"].to_numpy(), layer, :]
            scaler = StandardScaler(with_mean=False)
            Xtr = scaler.fit_transform(Xtr)
            Xte = scaler.transform(Xte)
            clf = LogisticRegression(max_iter=1000, C=1.0, random_state=args.seed)
            clf.fit(Xtr, train["cell_type"].values)
            pred = clf.predict(Xte)
            f1 = f1_score(test["cell_type"].values, pred, average="macro", zero_division=0)
            rows.append({
                "probe_name": "layer1_per_layer",
                "model": args.model,
                "layer": int(layer),
                "fold": int(fold),
                "metric": "macro_f1",
                "value": float(f1),
                "n_train": int(len(train)),
                "n_test": int(len(test)),
                "seed": args.seed,
            })
        layer_f1 = pd.DataFrame(rows)[lambda x: x["layer"] == layer]["value"]
        print(f"  layer {layer:>2}: macro_f1 mean={layer_f1.mean():.4f} median={layer_f1.median():.4f}")

    out = RESULTS / "probes" / f"layer1_per_layer__{args.model}.parquet"
    pd.DataFrame(rows).to_parquet(out, index=False)
    print(f"wrote → {out}")


if __name__ == "__main__":
    main()
