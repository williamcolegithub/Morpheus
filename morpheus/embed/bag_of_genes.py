"""Bag-of-genes baseline (PROBE-18).

- Cell embedding: vocab-restricted log1p-CP10K, stored sparse on disk (.npz).
  Probes load this and densify per fold-batch.
- Gene embedding: each gene's expression pattern across cells, then PCA-projected
  to 256 dims for fair comparison with Geneformer's 256-dim space.
"""
from __future__ import annotations

import json
import pickle
from datetime import UTC, datetime

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack, save_npz
from sklearn.decomposition import TruncatedSVD
from tqdm import tqdm

from morpheus.paths import CELLS_H5AD, EMBEDDINGS, GENEFORMER_DIR

GENE_EMB_DIM = 256


def main() -> None:
    out_dir = EMBEDDINGS / "bag_of_genes"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl", "rb") as f:
        ensembl2token = {k: int(v) for k, v in pickle.load(f).items() if k.startswith("ENSG")}

    print("loading cells.h5ad header ...")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    n_cells = adata.n_obs
    var_ensembl = adata.var_names.values
    var_in_vocab = np.array([e in ensembl2token for e in var_ensembl])
    vocab_var_idx = np.flatnonzero(var_in_vocab)
    vocab_ensembl = var_ensembl[vocab_var_idx]
    n_vocab = int(vocab_var_idx.size)
    print(f"  n_cells={n_cells:,}, vocab genes={n_vocab:,}")

    chunk = 5000
    sparse_chunks = []
    cell_ids: list[str] = []
    for s in tqdm(range(0, n_cells, chunk), desc="bag-of-genes"):
        e = min(s + chunk, n_cells)
        X = adata.X[s:e]
        if not isinstance(X, csr_matrix):
            X = csr_matrix(X)
        X = X[:, vocab_var_idx]  # already log1p-CP10K
        sparse_chunks.append(X)
        cell_ids.extend(adata.obs["cell_id"].iloc[s:e].astype(str).tolist())
    cells_sparse = vstack(sparse_chunks).tocsr()
    print(f"  cells_sparse: {cells_sparse.shape}, nnz={cells_sparse.nnz:,}")

    save_npz(out_dir / "cells.npz", cells_sparse)
    pd.DataFrame({"row": np.arange(len(cell_ids), dtype=np.int64), "cell_id": cell_ids}).to_parquet(
        out_dir / "cell_index.parquet", index=False
    )

    # Gene embeddings via TruncatedSVD on (n_genes, n_cells) — i.e. each row is a gene's pattern across cells.
    # Equivalent to PCA on transposed expression matrix (mean-centered).
    print("computing gene embeddings via TruncatedSVD ...")
    gene_matrix = cells_sparse.T  # (n_genes, n_cells), CSR after .T → CSC; convert.
    gene_matrix = gene_matrix.tocsr()
    svd = TruncatedSVD(n_components=GENE_EMB_DIM, random_state=0)
    gene_emb = svd.fit_transform(gene_matrix).astype(np.float32)
    print(f"  gene_emb: {gene_emb.shape}, explained_variance_ratio sum={svd.explained_variance_ratio_.sum():.3f}")
    np.save(out_dir / "genes.npy", gene_emb)
    pd.DataFrame({"row": np.arange(n_vocab, dtype=np.int64), "ensembl_id": vocab_ensembl}).to_parquet(
        out_dir / "gene_index.parquet", index=False
    )

    meta = {
        "model": "bag_of_genes",
        "weights_source": "log1p_cp10k_per_cell + truncated_svd_gene_pattern",
        "embedding_dim": GENE_EMB_DIM,
        "cell_embedding_format": "sparse_csr",
        "n_cells": int(n_cells),
        "vocab_genes": n_vocab,
        "explained_variance_ratio_sum": float(svd.explained_variance_ratio_.sum()),
        "extracted_at": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    (out_dir / "META.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
