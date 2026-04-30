"""Build donor-stratified K-fold splits. PROBE-20 audit lives here.

Every cell from a given donor lives in exactly one fold.
"""
from __future__ import annotations

import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd

from morpheus.paths import CELLS_H5AD, SPLITS


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    obs = adata.obs[["cell_id", "donor_id"]].copy().reset_index(drop=True)

    rng = np.random.default_rng(args.seed)
    donors = obs["donor_id"].astype(str).unique()
    rng.shuffle(donors)

    # Greedy bin packing: put each donor (largest first) into the fold with the smallest current cell count.
    sizes = obs["donor_id"].astype(str).value_counts().to_dict()
    donors_by_size = sorted(donors, key=lambda d: -sizes[d])
    fold_counts = np.zeros(args.k, dtype=int)
    donor_to_fold: dict[str, int] = {}
    for d in donors_by_size:
        f = int(np.argmin(fold_counts))
        donor_to_fold[d] = f
        fold_counts[f] += sizes[d]

    obs["fold"] = obs["donor_id"].astype(str).map(donor_to_fold).astype(int)
    splits = obs[["cell_id", "fold"]].copy()
    splits["cell_id"] = splits["cell_id"].astype(str)
    splits.to_parquet(SPLITS, index=False)
    print(f"wrote → {SPLITS}")

    # Audit: every donor in exactly one fold.
    audit = obs.groupby("donor_id")["fold"].nunique()
    bad = audit[audit > 1]
    if len(bad) > 0:
        raise AssertionError(f"DONOR LEAKAGE: {len(bad)} donors span multiple folds: {bad.head().to_dict()}")
    print(f"PROBE-20 audit OK: all {len(audit)} donors live in exactly one fold")
    print("fold sizes:", fold_counts.tolist())
    print(json.dumps({"k": args.k, "seed": args.seed, "fold_sizes": fold_counts.tolist(), "n_donors": len(donors)}, indent=2))


if __name__ == "__main__":
    main()
