"""Canonical paths. Single source of truth — see shared/contracts.md §file-system layout."""
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

DATA_RAW = REPO / "data" / "raw"
DATA_PROCESSED = REPO / "data" / "processed"
EMBEDDINGS = REPO / "data" / "embeddings"
RESULTS = REPO / "results"

CELLS_H5AD = DATA_PROCESSED / "cells.h5ad"
GENE_ID_MAP = DATA_PROCESSED / "gene_id_map.parquet"
DOROTHEA_EDGES = DATA_PROCESSED / "dorothea_edges.parquet"
SPLITS = DATA_PROCESSED / "splits.parquet"

GENEFORMER_DIR = DATA_RAW / "geneformer_weights"


def embed_dir(model: str) -> Path:
    return EMBEDDINGS / model


for d in (DATA_RAW, DATA_PROCESSED, EMBEDDINGS, RESULTS / "probes", RESULTS / "rsa", RESULTS / "figures"):
    d.mkdir(parents=True, exist_ok=True)
