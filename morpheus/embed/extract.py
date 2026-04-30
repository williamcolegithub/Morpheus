"""PROBE-9 + PROBE-10 + PROBE-11: Geneformer cell + gene + attention extraction.

- Cell embeddings: mean-pool of last-hidden-state over non-pad tokens.
- Gene embeddings: input token embeddings (the learned lookup table itself).
- Attention sample: stratified-by-cell-type sample of full attention tensors (~1k cells).

Streaming: iterates cells.h5ad in chunks, never densifies the full matrix.
Supports `--weights random_init` for the random-init control (PROBE-17).
"""
from __future__ import annotations

import argparse
import json
import pickle
import time
from datetime import UTC, datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.sparse import csr_matrix
from tqdm import tqdm
from transformers import BertConfig, BertForMaskedLM, BertModel

from morpheus.paths import CELLS_H5AD, EMBEDDINGS, GENEFORMER_DIR

PAD_TOKEN_ID = 0
MAX_LEN = 2048
DEFAULT_VARIANT = "Geneformer-V1-10M"


def _device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _load_dicts():
    base = GENEFORMER_DIR / "geneformer" / "gene_dictionaries_30m"
    with open(base / "token_dictionary_gc30M.pkl", "rb") as f:
        ensembl2token = {k: int(v) for k, v in pickle.load(f).items()}
    with open(base / "gene_median_dictionary_gc30M.pkl", "rb") as f:
        ensembl2median = {k: float(v) for k, v in pickle.load(f).items()}
    return ensembl2token, ensembl2median


def _load_model(weights: str, variant: str = DEFAULT_VARIANT, attn_implementation: str = "sdpa") -> BertModel:
    config_path = GENEFORMER_DIR / variant / "config.json"
    config = BertConfig.from_pretrained(config_path)
    config.attn_implementation = attn_implementation
    if weights == "random_init":
        model = BertForMaskedLM(config)
    else:
        model = BertForMaskedLM.from_pretrained(GENEFORMER_DIR / variant, attn_implementation=attn_implementation)
    return model.bert


def _tokenize_chunk(
    cp10k_chunk: np.ndarray,
    var_tokens: np.ndarray,
    var_medians: np.ndarray,
    var_in_vocab: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Tokenize a (chunk, n_genes) CP10K matrix.

    Returns (input_ids: int64 (chunk, MAX_LEN), attention_mask: int64).
    """
    chunk_size = cp10k_chunk.shape[0]
    vocab_tokens = var_tokens[var_in_vocab]
    vocab_medians = var_medians[var_in_vocab]
    scaled = cp10k_chunk[:, var_in_vocab] / np.clip(vocab_medians, 1e-9, None)
    n_vocab = scaled.shape[1]
    take = min(MAX_LEN, n_vocab)
    # Top-K indices per row via argpartition, then sort within those K.
    if take < n_vocab:
        part = np.argpartition(-scaled, kth=take - 1, axis=1)[:, :take]
    else:
        part = np.tile(np.arange(n_vocab), (chunk_size, 1))
    rows = np.arange(chunk_size)[:, None]
    top_vals = scaled[rows, part]
    order = np.argsort(-top_vals, axis=1)
    sorted_idx = part[rows, order]
    sorted_vals = top_vals[rows, order]
    valid = sorted_vals > 0  # only nonzero genes count
    tok_matrix = vocab_tokens[sorted_idx]  # (chunk, take)
    input_ids = np.zeros((chunk_size, MAX_LEN), dtype=np.int64)
    attention_mask = np.zeros((chunk_size, MAX_LEN), dtype=np.int64)
    input_ids[:, :take] = np.where(valid, tok_matrix, 0)
    attention_mask[:, :take] = valid.astype(np.int64)
    return input_ids, attention_mask


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--weights", default=DEFAULT_VARIANT, help="model variant dir name OR 'random_init'")
    p.add_argument("--out", default=None, help="output dir under data/embeddings/. Defaults to model name.")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=1000, help="cells per tokenization chunk")
    p.add_argument("--attention-sample", type=int, default=200, help="how many cells to keep full attention for")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-cells", type=int, default=None, help="optional cap for smoke testing")
    args = p.parse_args()

    device = _device()
    print(f"device: {device}")

    out_name = args.out or ("geneformer" if args.weights == DEFAULT_VARIANT else ("geneformer_random_init" if args.weights == "random_init" else args.weights))
    out_dir = EMBEDDINGS / out_name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading model: weights={args.weights}, out={out_dir}")
    model = _load_model(args.weights, attn_implementation="sdpa")
    model.eval().to(device)
    embed_dim = model.config.hidden_size
    print(f"  hidden_size={embed_dim}, layers={model.config.num_hidden_layers}, heads={model.config.num_attention_heads}")

    # Save the gene token embeddings (input embedding table) — used by Layer 2/3 probes.
    word_emb = model.embeddings.word_embeddings.weight.detach().cpu().to(torch.float32).numpy()
    gene_index_rows = []
    ensembl2token, _ = _load_dicts()
    for ensembl_id, tok in ensembl2token.items():
        if ensembl_id.startswith("ENSG"):
            gene_index_rows.append({"row": int(tok), "ensembl_id": ensembl_id})
    gene_index = pd.DataFrame(gene_index_rows).sort_values("row").reset_index(drop=True)
    np.save(out_dir / "genes.npy", word_emb)
    gene_index.to_parquet(out_dir / "gene_index.parquet", index=False)
    print(f"  saved genes.npy {word_emb.shape} + gene_index.parquet ({len(gene_index)} rows)")

    print("loading cells.h5ad header ...")
    adata = ad.read_h5ad(CELLS_H5AD, backed="r")
    n_cells = adata.n_obs if args.max_cells is None else min(args.max_cells, adata.n_obs)
    n_genes = adata.n_vars
    print(f"  n_cells={n_cells:,}, n_genes={n_genes:,}")

    # Build per-var arrays once.
    var_ensembl = adata.var_names.values  # ensembl IDs
    var_tokens = np.array([ensembl2token.get(e, -1) for e in var_ensembl], dtype=np.int64)
    var_in_vocab = var_tokens >= 2  # 0/1 are pad/mask
    _, ensembl2median = _load_dicts()
    var_medians = np.array([ensembl2median.get(e, 1.0) for e in var_ensembl], dtype=np.float64)
    print(f"  vars in Geneformer vocab: {var_in_vocab.sum():,} / {n_genes:,}")

    # Pre-allocate cell embedding memmap.
    cells_path = out_dir / "cells.npy"
    cells_arr = np.lib.format.open_memmap(
        cells_path, mode="w+", dtype=np.float32, shape=(n_cells, embed_dim)
    )
    cell_ids_out: list[str] = []

    # Attention sampling: pick stratified subset of cells up front.
    rng = np.random.default_rng(args.seed)
    obs_meta = adata.obs[["cell_id", "cell_type"]].iloc[:n_cells].copy().reset_index(drop=True)
    if args.attention_sample > 0:
        per_type = max(1, args.attention_sample // max(1, obs_meta["cell_type"].nunique()))
        att_idx_list = []
        for _, grp in obs_meta.groupby("cell_type"):
            take = min(per_type, len(grp))
            att_idx_list.append(rng.choice(grp.index.to_numpy(), size=take, replace=False))
        att_idx = np.sort(np.concatenate(att_idx_list))[: args.attention_sample]
        att_idx_set = set(int(i) for i in att_idx)
    else:
        att_idx_set = set()

    att_records: list[dict] = []  # accumulate {cell_id, attention np.ndarray[L,H,T,T]}

    t0 = time.time()
    with torch.no_grad():
        for chunk_start in tqdm(range(0, n_cells, args.chunk_size), desc="chunks"):
            chunk_end = min(chunk_start + args.chunk_size, n_cells)
            X_chunk = adata.X[chunk_start:chunk_end]
            if isinstance(X_chunk, csr_matrix):
                X_chunk_dense = X_chunk.toarray().astype(np.float32)
            else:
                X_chunk_dense = np.asarray(X_chunk, dtype=np.float32)
            # Undo log1p → CP10K.
            cp10k = np.expm1(X_chunk_dense)
            input_ids, attention_mask = _tokenize_chunk(cp10k, var_tokens, var_medians, var_in_vocab)

            # Run model in mini-batches.
            for bs_start in range(0, chunk_end - chunk_start, args.batch_size):
                bs_end = min(bs_start + args.batch_size, chunk_end - chunk_start)
                ids_np = input_ids[bs_start:bs_end]
                mask_np = attention_mask[bs_start:bs_end]
                # Dynamic truncation: shrink to max actual length in batch.
                max_len = int(mask_np.sum(axis=1).max()) if mask_np.size else 1
                max_len = max(max_len, 1)
                ids_np = ids_np[:, :max_len]
                mask_np = mask_np[:, :max_len]
                ids = torch.from_numpy(ids_np).to(device)
                mask = torch.from_numpy(mask_np).to(device)
                out = model(input_ids=ids, attention_mask=mask)
                last = out.last_hidden_state
                m = mask.unsqueeze(-1).to(last.dtype)
                summed = (last * m).sum(dim=1)
                lengths = m.sum(dim=1).clamp(min=1.0)
                pooled = (summed / lengths).cpu().to(torch.float32).numpy()
                cells_arr[chunk_start + bs_start : chunk_start + bs_end] = pooled

            cell_ids_out.extend(obs_meta["cell_id"].iloc[chunk_start:chunk_end].astype(str).tolist())

    cells_arr.flush()
    del cells_arr

    cell_index = pd.DataFrame({"row": np.arange(len(cell_ids_out), dtype=np.int64), "cell_id": cell_ids_out})
    cell_index.to_parquet(out_dir / "cell_index.parquet", index=False)

    # Attention pass — small, eager, dynamic-truncated per cell. Saved as a list of variable-shape arrays
    # to avoid the 40 GB cost of padding everything to (T=2048).
    if att_idx_set:
        ATT_MAX_TOKENS = 256  # cap token axis to keep size manageable; first 256 ranked tokens carry the bulk of signal
        print(f"attention pass: {len(att_idx_set)} cells with eager attention (ATT_MAX_TOKENS={ATT_MAX_TOKENS}) ...")
        del model
        if device.type == "mps":
            torch.mps.empty_cache()
        att_model = _load_model(args.weights, attn_implementation="eager").eval().to(device)
        sorted_att_idx = sorted(att_idx_set)
        att_records = []
        with torch.no_grad():
            for i_start in tqdm(range(0, len(sorted_att_idx), args.batch_size), desc="attention"):
                gi_batch = sorted_att_idx[i_start : i_start + args.batch_size]
                X_rows = adata.X[gi_batch]
                if isinstance(X_rows, csr_matrix):
                    X_dense = X_rows.toarray().astype(np.float32)
                else:
                    X_dense = np.asarray(X_rows, dtype=np.float32)
                cp10k = np.expm1(X_dense)
                ids_b, mask_b = _tokenize_chunk(cp10k, var_tokens, var_medians, var_in_vocab)
                # Dynamic truncate to min(actual_max_len_in_batch, ATT_MAX_TOKENS).
                max_len = min(int(mask_b.sum(axis=1).max()) if mask_b.size else 1, ATT_MAX_TOKENS)
                max_len = max(max_len, 1)
                ids_b = ids_b[:, :max_len]
                mask_b = mask_b[:, :max_len]
                ids = torch.from_numpy(ids_b).to(device)
                mask = torch.from_numpy(mask_b).to(device)
                out = att_model(input_ids=ids, attention_mask=mask, output_attentions=True)
                att_stacked = torch.stack(out.attentions, dim=1).to(torch.float16).cpu().numpy()  # (B, L, H, T, T)
                for k, gi in enumerate(gi_batch):
                    att_records.append({
                        "cell_id": str(obs_meta["cell_id"].iloc[gi]),
                        "attention": att_stacked[k],
                        "input_ids": ids_b[k].copy(),
                    })
        # All entries have the same T (since we cap at ATT_MAX_TOKENS and pad to batch max), but T can differ across batches.
        # Save as object array of (cell_id, attention, input_ids) tuples — readable via np.load(allow_pickle=True).
        att_cell_ids = np.array([r["cell_id"] for r in att_records], dtype=object)
        att_attention = np.array([r["attention"] for r in att_records], dtype=object)
        att_input_ids = np.array([r["input_ids"] for r in att_records], dtype=object)
        np.savez_compressed(
            out_dir / "attention_sample.npz",
            cell_ids=att_cell_ids,
            attention=att_attention,
            input_ids=att_input_ids,
        )
        sizes = [r["attention"].shape for r in att_records[:5]]
        print(f"  saved attention_sample.npz: n={len(att_records)}; first shapes {sizes}")

    elapsed = time.time() - t0
    meta = {
        "model": out_name,
        "weights_source": args.weights,
        "variant": DEFAULT_VARIANT,
        "embedding_dim": int(embed_dim),
        "n_cells": int(n_cells),
        "vocab_genes": int(var_in_vocab.sum()),
        "extracted_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "device": str(device),
        "batch_size": int(args.batch_size),
        "chunk_size": int(args.chunk_size),
        "attention_sampling": {"strategy": "stratified_by_cell_type", "n": len(att_records), "seed": args.seed},
        "elapsed_sec": round(elapsed, 2),
    }
    (out_dir / "META.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote → {out_dir}")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
