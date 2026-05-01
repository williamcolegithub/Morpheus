"""PROBE-7: build per-variant gene_id_map and resolve DoRothEA edge endpoints to Ensembl IDs.

Reads:
  - cells.h5ad  (var indexed by Ensembl ID; var['hgnc_symbol'])
  - Geneformer vocab dictionaries for the chosen --variant (gc30M for V1, gc104M for V2)
  - dorothea_edges.prevocab.parquet  (TF-target edges by HGNC symbol)
Writes:
  - gene_id_map_<variant>.parquet — vocab-membership flags for that variant (provenance)
  - dorothea_edges.parquet — all edges with both Ensembl IDs resolved; NO vocab filter applied.
    Per-model vocab restriction happens at probe time via each model's gene_index.parquet.
    This makes V1/V2 probes draw from the same canonical edge set.
"""
from __future__ import annotations

import argparse
import json
import pickle

import anndata as ad
import pandas as pd

from morpheus.paths import CELLS_H5AD, DATA_PROCESSED, DOROTHEA_EDGES, GENEFORMER_DIR


def _vocab_token_path(variant: str):
    if variant.startswith("Geneformer-V2"):
        return GENEFORMER_DIR / "geneformer" / "token_dictionary_gc104M.pkl"
    return GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl"


def _load_variant_token_map(variant: str) -> dict[str, int]:
    with open(_vocab_token_path(variant), "rb") as f:
        return {k: int(v) for k, v in pickle.load(f).items() if k.startswith("ENSG")}


def _load_unioned_name2ensembl() -> dict[str, str]:
    """Symbol → Ensembl, unioning every Geneformer name dictionary present on disk.

    The two gc30M and gc104M maps disagree slightly because they were built against
    different GENCODE releases. Taking the union (newer wins on collision) gives the
    most permissive resolution and keeps V1/V2 probes reading from one canonical set
    of edges.
    """
    base = GENEFORMER_DIR / "geneformer"
    candidates = [
        base / "gene_dictionaries_30m" / "gene_name_id_dict_gc30M.pkl",
        base / "gene_name_id_dict_gc104M.pkl",  # later in list → wins on collision
    ]
    merged: dict[str, str] = {}
    for path in candidates:
        if not path.exists():
            continue
        with open(path, "rb") as f:
            merged.update(pickle.load(f))
    if not merged:
        raise FileNotFoundError(
            "No Geneformer name→ensembl dict found. "
            "Pull at least one of gene_name_id_dict_gc{30M,104M}.pkl into data/raw/geneformer_weights/geneformer/."
        )
    return merged


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--variant", default="Geneformer-V1-10M", help="model variant whose vocab to flag membership against")
    args = p.parse_args()

    print(f"variant: {args.variant}")
    print("loading cells.h5ad header ...")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    cells_genes = adata.var.copy()
    if "ensembl_id" not in cells_genes.columns:
        cells_genes.index.name = "ensembl_id"
        cells_genes = cells_genes.reset_index()
    else:
        cells_genes = cells_genes.reset_index(drop=True)
    print(f"  cells.h5ad has {len(cells_genes):,} genes")

    ensembl2token = _load_variant_token_map(args.variant)
    name2ensembl = _load_unioned_name2ensembl()
    print(f"  variant vocab: {len(ensembl2token):,} ensembl IDs; unioned name→id table: {len(name2ensembl):,}")

    edges_pre = pd.read_parquet(DOROTHEA_EDGES.with_suffix(".prevocab.parquet"))
    print(f"  prevocab edges: {len(edges_pre):,}")

    cells_genes["in_geneformer_vocab"] = cells_genes["ensembl_id"].isin(ensembl2token)
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
    map_path = DATA_PROCESSED / f"gene_id_map_{args.variant}.parquet"
    map_df.to_parquet(map_path, index=False)
    print(f"wrote → {map_path}")

    # Resolve edges to Ensembl IDs (for both endpoints) but DO NOT vocab-filter.
    edges = edges_pre.copy()
    edges["tf_ensembl"] = edges["tf_symbol"].map(name2ensembl)
    edges["target_ensembl"] = edges["target_symbol"].map(name2ensembl)
    n_resolved = len(edges)
    edges = edges.dropna(subset=["tf_ensembl", "target_ensembl"])
    edges = edges[["tf_ensembl", "tf_symbol", "target_ensembl", "target_symbol", "confidence", "mor", "source"]].drop_duplicates()
    edges.to_parquet(DOROTHEA_EDGES, index=False)
    print(f"  edges with both endpoints resolvable to Ensembl: {len(edges):,} (was {n_resolved:,})")
    print(f"wrote → {DOROTHEA_EDGES}")
    print(f"  source dist:\n{edges['source'].value_counts().to_string()}")
    print(f"  unique TFs in resolved edges: {edges['tf_ensembl'].nunique()}")

    summary = {
        "variant": args.variant,
        "n_genes_cells": int(len(cells_genes)),
        "n_in_variant_vocab": int(cells_genes["in_geneformer_vocab"].sum()),
        "n_in_dorothea": int(cells_genes["in_dorothea"].sum()),
        "n_edges_resolved": int(len(edges)),
        "n_unique_tfs": int(edges["tf_ensembl"].nunique()),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
