"""PROBE-4 + PROBE-5: pull blood/immune subset from CELLxGENE Census.

v1 scale: 100k cells, healthy donors only, protein-coding genes only, log-normalized,
saved as cells.h5ad (sparse). Donor-stratified — we sample by donor so full donors stay
in the dataset (needed for clean donor-stratified splits later).
"""
from __future__ import annotations

import json
from datetime import UTC, datetime

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

from morpheus.paths import CELLS_H5AD, DATA_RAW, DATA_PROCESSED

CENSUS_VERSION = "stable"
TARGET_CELLS = 100_000
SEED = 0

OBS_FILTER = (
    "tissue_general == 'blood' "
    "and disease == 'normal' "
    "and is_primary_data == True "
    "and assay in ['10x 3\\' v2', '10x 3\\' v3', '10x 5\\' v1', '10x 5\\' v2']"
)
VAR_FILTER = None  # Census var has no biotype column on this version; protein-coding restriction happens via Geneformer vocab in build_gene_map.


def main() -> None:
    print(f"opening Census version={CENSUS_VERSION} ...")
    raw_dir = DATA_RAW / "cellxgene"
    raw_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Phase 1: get just the obs metadata so we can sample donors before pulling X.
        print("phase 1: pulling obs metadata for blood/normal cells ...")
        obs_df = (
            census["census_data"]["homo_sapiens"]
            .obs.read(
                value_filter=OBS_FILTER,
                column_names=["soma_joinid", "donor_id", "cell_type", "cell_type_ontology_term_id", "dataset_id", "assay"],
            )
            .concat()
            .to_pandas()
        )
        print(f"  matched obs rows: {len(obs_df):,}; donors: {obs_df['donor_id'].nunique():,}; datasets: {obs_df['dataset_id'].nunique():,}")

        # Phase 2: donor-stratified sampling. Keep full donors where possible; trim only inside the boundary donor.
        donor_sizes = obs_df.groupby("donor_id").size().sample(frac=1.0, random_state=SEED)
        cum = donor_sizes.cumsum()
        # Donors strictly under target...
        under = donor_sizes.index[cum < TARGET_CELLS].tolist()
        # ...plus one more donor that crosses the boundary.
        crossing_idx = (cum >= TARGET_CELLS).idxmax() if (cum >= TARGET_CELLS).any() else None
        keep_donors = under + ([crossing_idx] if crossing_idx else [])
        if not keep_donors:
            keep_donors = [donor_sizes.index[0]]
        sampled_obs = obs_df[obs_df["donor_id"].isin(keep_donors)].copy()
        # Trim only inside the boundary donor — keep all earlier donors whole.
        if len(sampled_obs) > TARGET_CELLS and crossing_idx is not None:
            other = sampled_obs[sampled_obs["donor_id"] != crossing_idx]
            need = TARGET_CELLS - len(other)
            crossing_rows = sampled_obs[sampled_obs["donor_id"] == crossing_idx]
            keep_idx = rng.choice(crossing_rows.index, size=max(need, 0), replace=False)
            sampled_obs = pd.concat([other, crossing_rows.loc[keep_idx]])
        print(f"  sampled cells: {len(sampled_obs):,} across {sampled_obs['donor_id'].nunique()} donors")

        soma_ids = sampled_obs["soma_joinid"].to_numpy()

        # Phase 3: fetch X + var for sampled cells.
        print("phase 3: pulling AnnData for sampled cells ...")
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_coords=soma_ids,
            obs_column_names=["soma_joinid", "donor_id", "cell_type", "cell_type_ontology_term_id", "dataset_id", "assay", "sex", "self_reported_ethnicity"],
            var_column_names=["feature_id", "feature_name", "feature_type", "feature_length"],
        )

    print(f"  raw AnnData: {adata.shape}, X density={adata.X.nnz / (adata.shape[0] * adata.shape[1]):.4f}")

    # Restrict to protein-coding genes. CELLxGENE uses 'feature_type' or biotype is in feature_biotype which we filtered to 'gene' — but 'gene' may include lncRNA. Use feature_name + biotype if available; fallback: use Geneformer vocab match later anyway.
    # Census var has 'feature_biotype' = 'gene' for everything matching our filter; the actual biotype detail isn't exposed there. We'll restrict to protein-coding via Geneformer vocab in build_gene_map. For now, keep all 'gene' biotype.

    # Standardize cell_id (we use soma_joinid prefixed with the dataset short id for uniqueness, but soma_joinid is already unique within a Census version — use it directly).
    adata.obs["cell_id"] = adata.obs["soma_joinid"].astype(str).values
    adata.obs_names = adata.obs["cell_id"].astype(str).values
    # Use Ensembl ID (feature_id) as the canonical var index.
    adata.var["ensembl_id"] = adata.var["feature_id"].values
    adata.var_names = adata.var["ensembl_id"].astype(str).values
    adata.var["hgnc_symbol"] = adata.var["feature_name"].fillna("").values

    # Friendly cell type name (we already have cell_type from Census; cell_type_ontology_term_id is the CL id).
    adata.obs["cell_type_name"] = adata.obs["cell_type"].astype(str).values
    adata.obs["cell_type"] = adata.obs["cell_type_ontology_term_id"].astype(str).values  # rename for contract: .obs["cell_type"] holds CL id

    # log1p normalize. Census X is raw counts.
    print("normalizing: total_count → 10k, log1p ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print(f"final AnnData: {adata.shape}, donors={adata.obs['donor_id'].nunique()}, cell_types={adata.obs['cell_type'].nunique()}")

    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(CELLS_H5AD, compression="gzip")
    print(f"wrote → {CELLS_H5AD}")

    (raw_dir / "PROVENANCE.md").write_text(
        f"Source: cellxgene_census v{cellxgene_census.__version__}, census_version='{CENSUS_VERSION}'\n"
        f"Filter: {OBS_FILTER}\n"
        f"Var filter: {VAR_FILTER}\n"
        f"Target cells: {TARGET_CELLS}; achieved: {adata.n_obs}\n"
        f"Donors: {adata.obs['donor_id'].nunique()}; datasets: {adata.obs['dataset_id'].nunique()}\n"
        f"Cell types (CL): {adata.obs['cell_type'].nunique()}\n"
        f"Seed: {SEED}\n"
        f"Access date (UTC): {datetime.now(UTC).isoformat(timespec='seconds')}\n"
    )

    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_donors": int(adata.obs["donor_id"].nunique()),
        "n_cell_types": int(adata.obs["cell_type"].nunique()),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
