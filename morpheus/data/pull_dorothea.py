"""PROBE-6: pull DoRothEA edges via decoupler. Vocab filter applied later in build_gene_map."""
from __future__ import annotations

import json
from datetime import UTC, datetime

import decoupler as dc
import pandas as pd

from morpheus.paths import DATA_RAW, DOROTHEA_EDGES


def _normalize(df: pd.DataFrame, source: str) -> pd.DataFrame:
    df = df.rename(columns={"source": "tf_symbol", "target": "target_symbol"})
    if "mor" not in df.columns:
        if "weight" in df.columns:
            df["mor"] = df["weight"].apply(lambda w: 1 if w > 0 else (-1 if w < 0 else 0)).astype(int)
        else:
            df["mor"] = 0
    if "confidence" not in df.columns:
        df["confidence"] = "A"
    df["source"] = source
    return df[["tf_symbol", "target_symbol", "confidence", "mor", "source"]].drop_duplicates()


def main() -> None:
    print("pulling DoRothEA (human, levels A/B/C) via decoupler.op.dorothea ...")
    dor = _normalize(dc.op.dorothea(organism="human", levels=["A", "B", "C"], verbose=True), "dorothea")
    print(f"DoRothEA edges: {len(dor)}")

    print("pulling CollecTRI via decoupler.op.collectri ...")
    col = _normalize(dc.op.collectri(organism="human", verbose=True), "collectri")
    print(f"CollecTRI edges: {len(col)}")

    df = pd.concat([dor, col], ignore_index=True).drop_duplicates(subset=["tf_symbol", "target_symbol", "source"])
    print("union shape:", df.shape)
    print("source dist:\n", df["source"].value_counts().to_string())
    print("confidence dist:\n", df["confidence"].value_counts().to_string())

    raw_dir = DATA_RAW / "dorothea"
    raw_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(raw_dir / "edges_union.parquet", index=False)
    (raw_dir / "PROVENANCE.md").write_text(
        f"Sources: dc.op.dorothea(human, levels=A/B/C) ∪ dc.op.collectri(human)\n"
        f"decoupler version: {dc.__version__}\n"
        f"Access date (UTC): {datetime.now(UTC).isoformat(timespec='seconds')}\n"
        f"DoRothEA edges: {len(dor)}; CollecTRI edges: {len(col)}; union: {len(df)}\n"
    )

    tmp_path = DOROTHEA_EDGES.with_suffix(".prevocab.parquet")
    df.to_parquet(tmp_path, index=False)
    print(f"wrote pre-vocab edges → {tmp_path}")

    # Also dump column schema for reference.
    print(json.dumps({"columns": list(df.columns), "n_rows": len(df)}, indent=2))


if __name__ == "__main__":
    main()
