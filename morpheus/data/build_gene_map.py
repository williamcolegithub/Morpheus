"""PROBE-7: build gene_id_map.parquet, apply Geneformer-vocab filter to dorothea_edges.

Reads:
  - cells.h5ad  (var indexed by Ensembl ID; var['hgnc_symbol'])
  - Geneformer token_dictionary_gc30M.pkl  (Ensembl ID → token int)
  - dorothea_edges.prevocab.parquet  (TF-target edges by HGNC symbol)
Writes:
  - gene_id_map.parquet
  - dorothea_edges.parquet (vocab-filtered, with both Ensembl IDs resolved)
"""
from __future__ import annotations

import json
import pickle

import anndata as ad
import pandas as pd

from morpheus.paths import CELLS_H5AD, DOROTHEA_EDGES, GENE_ID_MAP, GENEFORMER_DIR


def _load_geneformer_vocab() -> tuple[dict[str, int], dict[str, str]]:
    vocab_path = GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl"
    name2id_path = GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "gene_name_id_dict_gc30M.pkl"
    with open(vocab_path, "rb") as f:
        ensembl2token = {k: int(v) for k, v in pickle.load(f).items() if k.startswith("ENSG")}
    with open(name2id_path, "rb") as f:
        name2ensembl = pickle.load(f)
    return ensembl2token, name2ensembl


def main() -> None:
    print("loading cells.h5ad header ...")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    cells_genes = adata.var.copy()
    if "ensembl_id" not in cells_genes.columns:
        cells_genes.index.name = "ensembl_id"
        cells_genes = cells_genes.reset_index()
    else:
        cells_genes = cells_genes.reset_index(drop=True)
    print(f"  cells.h5ad has {len(cells_genes):,} genes")

    print("loading Geneformer vocab ...")
    ensembl2token, name2ensembl = _load_geneformer_vocab()
    print(f"  vocab: {len(ensembl2token):,} ensembl IDs; name→id table: {len(name2ensembl):,}")

    edges_pre = pd.read_parquet(DOROTHEA_EDGES.with_suffix(".prevocab.parquet"))
    print(f"  pre-filter edges: {len(edges_pre):,}")

    # Build gene_id_map.
    cells_genes["in_geneformer_vocab"] = cells_genes["ensembl_id"].isin(ensembl2token)
    # Mark genes that appear as TF or target in the pre-vocab edge list (resolve their ensembl via name2ensembl).
    edge_symbols = set(edges_pre["tf_symbol"]) | set(edges_pre["target_symbol"])
    edge_ensembls = {name2ensembl[s] for s in edge_symbols if s in name2ensembl}
    cells_genes["in_dorothea"] = cells_genes["ensembl_id"].isin(edge_ensembls)
    cells_genes["hgnc_symbol"] = cells_genes["hgnc_symbol"].astype(str).fillna("")

    print(
        f"  in vocab: {cells_genes['in_geneformer_vocab'].sum():,}; "
        f"in dorothea: {cells_genes['in_dorothea'].sum():,}; "
        f"both: {(cells_genes['in_geneformer_vocab'] & cells_genes['in_dorothea']).sum():,}"
    )

    map_df = cells_genes[["ensembl_id", "hgnc_symbol", "in_geneformer_vocab", "in_dorothea"]].copy()
    map_df.to_parquet(GENE_ID_MAP, index=False)
    print(f"wrote → {GENE_ID_MAP}")

    # Vocab-filter edges. Symbol → Ensembl via name2ensembl, then drop edges where either side isn't in Geneformer vocab.
    edges = edges_pre.copy()
    edges["tf_ensembl"] = edges["tf_symbol"].map(name2ensembl)
    edges["target_ensembl"] = edges["target_symbol"].map(name2ensembl)
    pre_n = len(edges)
    edges = edges.dropna(subset=["tf_ensembl", "target_ensembl"])
    edges = edges[edges["tf_ensembl"].isin(ensembl2token) & edges["target_ensembl"].isin(ensembl2token)].copy()
    edges = edges[["tf_ensembl", "tf_symbol", "target_ensembl", "target_symbol", "confidence", "mor", "source"]].drop_duplicates()
    edges.to_parquet(DOROTHEA_EDGES, index=False)
    print(f"  edges after vocab filter: {len(edges):,} (was {pre_n:,})")
    print(f"wrote → {DOROTHEA_EDGES}")
    print(f"  source dist:\n{edges['source'].value_counts().to_string()}")
    print(f"  unique TFs in filtered edges: {edges['tf_ensembl'].nunique()}")

    # Mapping rate against cells.h5ad genes that have a non-empty hgnc_symbol.
    has_symbol = cells_genes[cells_genes["hgnc_symbol"] != ""]
    mappable = has_symbol["ensembl_id"].isin({name2ensembl[s] for s in name2ensembl if s in set(has_symbol["hgnc_symbol"])})
    rate = mappable.mean() if len(has_symbol) else 0.0
    print(f"  HGNC↔Ensembl mapping rate (cells.h5ad symbols → name2ensembl): {rate:.3f}")

    summary = {
        "n_genes_cells": int(len(cells_genes)),
        "n_in_vocab": int(cells_genes["in_geneformer_vocab"].sum()),
        "n_in_dorothea": int(cells_genes["in_dorothea"].sum()),
        "n_edges_post_vocab": int(len(edges)),
        "n_unique_tfs": int(edges["tf_ensembl"].nunique()),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
