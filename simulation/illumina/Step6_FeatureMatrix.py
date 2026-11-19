#!/usr/bin/env python3
import os
import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack, save_npz
from sklearn.feature_selection import chi2
from sklearn.model_selection import GroupKFold
import argparse
import sys

#############################################
# 1. HELPER FUNCTIONS
#############################################
def load_unique_kmers(unique_kmers_file):
    """Load unique k-mers (one per line or tab-delimited) and build {kmer->index}."""
    with open(unique_kmers_file, 'r') as f:
        unique_kmers_list = [line.split('\t')[0].strip().upper() for line in f if line.strip()]
    if not unique_kmers_list:
        raise ValueError(f"No k-mers found in {unique_kmers_file}")
    return {kmer: idx for idx, kmer in enumerate(unique_kmers_list)}

def extract_label_and_sample_id_from_parquet(parquet_filename):
    """
    Extract label and sample_id from a Parquet filename.
    Expected filename: 'batch_{batchnum}_{sample_id}.parquet'
    The 'sample_id' must contain 'Resistant' or 'Susceptible' to derive label.
    """
    base = os.path.basename(parquet_filename)
    base_no_ext = os.path.splitext(base)[0]
    parts = base_no_ext.split('_', 2)
    if len(parts) < 3:
        raise ValueError(f"Filename format unexpected: {parquet_filename}")
    sample_id = parts[2]
    if "Resistant" in sample_id:
        label = 1
    elif "Susceptible" in sample_id:
        label = 0
    else:
        # Allow unknown, but warn
        print(f"[WARN] Could not infer label from filename (no 'Resistant'/'Susceptible'): {parquet_filename}")
        label = -1
    return label, sample_id

#############################################
# 2. BUILD FEATURE MATRIX FROM PARQUET
#############################################
def create_feature_matrix_from_parquet(parquet_folder, unique_kmers_file, output_file):
    """
    Reads all Parquet files (columns: 'kmer','presence'), aggregates per sample,
    builds a sparse presence/absence matrix (columns=global unique k-mers), and
    saves matrix + labels to NPZs basing paths on `output_file`.
    """
    if not os.path.isdir(parquet_folder):
        raise ValueError(f"Parquet folder does not exist: {parquet_folder}")

    unique_kmers_index = load_unique_kmers(unique_kmers_file)
    unique_kmers_keys = set(unique_kmers_index.keys())

    parquet_files = glob.glob(os.path.join(parquet_folder, "*.parquet"))
    if not parquet_files:
        raise ValueError(f"No Parquet files found in {parquet_folder}")

    feature_matrix_dict = {}  # sample_id -> sparse row vector
    labels_dict = {}          # sample_id -> label
    sample_ids = []           # preserve order

    print(f"[INFO] Found {len(parquet_files)} parquet files in {parquet_folder}")

    for pf in parquet_files:
        try:
            label, sample_id = extract_label_and_sample_id_from_parquet(pf)
        except Exception as e:
            print(f"[WARN] Skipping file (name parse failed): {pf} | {e}")
            continue

        try:
            df = pd.read_parquet(pf)
        except Exception as e:
            print(f"[WARN] Skipping file (read failed): {pf} | {e}")
            continue

        if not {"kmer", "presence"}.issubset(df.columns):
            print(f"[WARN] Skipping file (missing columns): {pf}")
            continue

        # Clean and filter
        df["kmer"] = df["kmer"].astype(str).str.strip().str.upper()
        df = df[df["kmer"].isin(unique_kmers_keys)]
        if df.empty:
            continue

        # Map to column indices
        indices = df["kmer"].map(unique_kmers_index).values
        values = df["presence"].astype(np.uint8).values
        row = csr_matrix((values, ([0]*len(indices), indices)),
                         shape=(1, len(unique_kmers_index)), dtype=np.uint8)

        if sample_id in feature_matrix_dict:
            feature_matrix_dict[sample_id] += row
        else:
            feature_matrix_dict[sample_id] = row
            labels_dict[sample_id] = label
            sample_ids.append(sample_id)

    if not feature_matrix_dict:
        raise ValueError("After filtering, no data remained to build the feature matrix.")

    # Convert counts -> binary presence/absence
    for sid, row in feature_matrix_dict.items():
        if row.nnz > 0:
            row.data[:] = 1

    # Stack in sample_ids order
    feature_matrix = vstack([feature_matrix_dict[sid] for sid in sample_ids], format="csr")
    labels = np.array([labels_dict[sid] for sid in sample_ids], dtype=np.int8)
    sample_ids_array = np.array(sample_ids, dtype=object)

    # Save NPZs alongside the requested output_file base
    parent_dir = os.path.dirname(output_file)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)

    npz_file = output_file.replace(".csv", ".npz")
    labels_file = output_file.replace(".csv", "_labels.npz")

    save_npz(npz_file, feature_matrix)
    np.savez(labels_file, labels=labels, sample_ids=sample_ids_array)

    print(f"[INFO] Feature matrix saved to {npz_file}")
    print(f"[INFO] Labels+sample_ids saved to {labels_file}")

    return feature_matrix, labels, sample_ids_array, unique_kmers_index

#############################################
# 3. CV SPLIT & PER-FOLD SELECTION
#############################################
def split_data_and_process(
    feature_matrix,
    labels,
    sample_ids,
    output_folder,
    unique_kmers_index,
    top_k=1500,
    n_folds=5
):
    """
    Per fold:
      1) GroupKFold split (group by sample base to avoid leak across sims).
      2) chi2 on TRAIN only -> pick top_k features.
      3) Save train/test CSVs with those columns.
    """
    os.makedirs(output_folder, exist_ok=True)

    n_features_total = feature_matrix.shape[1]
    if n_features_total == 0:
        raise ValueError("Feature matrix has zero columns.")

    top_k = int(min(top_k, n_features_total))
    if top_k <= 0:
        raise ValueError(f"top_k={top_k} is not valid for {n_features_total} features.")

    reverse_index = {v: k for k, v in unique_kmers_index.items()}

    # Group by sample name without the trailing _sim# (best-effort)
    sample_names = np.array([sid.rsplit('_', 1)[0] for sid in sample_ids], dtype=object)

    gkf = GroupKFold(n_splits=n_folds)
    fold_number = 1

    for train_idx, test_idx in gkf.split(feature_matrix, labels, groups=sample_names):
        X_train_full = feature_matrix[train_idx]
        y_train      = labels[train_idx]
        X_test_full  = feature_matrix[test_idx]
        y_test       = labels[test_idx]

        chi_scores, _ = chi2(X_train_full, y_train)
        # Handle degenerate case (all zero chi2) by still picking indices deterministically
        order = np.argsort(chi_scores)
        top_feature_indices = np.sort(order[-top_k:])

        X_train_selected = X_train_full[:, top_feature_indices].toarray()
        X_test_selected  = X_test_full[:, top_feature_indices].toarray()

        selected_kmer_names = [reverse_index[i] for i in top_feature_indices]

        train_df = pd.DataFrame(X_train_selected, columns=selected_kmer_names)
        train_df["Label"] = y_train
        train_df["Sample_ID"] = sample_ids[train_idx]

        test_df = pd.DataFrame(X_test_selected, columns=selected_kmer_names)
        test_df["Label"] = y_test
        test_df["Sample_ID"] = sample_ids[test_idx]

        train_file = os.path.join(output_folder, f"fold_{fold_number}_train.csv")
        test_file  = os.path.join(output_folder, f"fold_{fold_number}_test.csv")
        train_df.to_csv(train_file, index=False)
        test_df.to_csv(test_file, index=False)
        print(f"[INFO] Fold {fold_number}: saved\n  {train_file}\n  {test_file}")
        fold_number += 1

#############################################
# 4. MAIN
#############################################
def main():
    parser = argparse.ArgumentParser(
        description="Build feature matrix from Parquet files and perform CV split with per-fold χ² selection."
    )
    parser.add_argument("--parquet_folder", required=True,
                        help="Directory containing Parquet files (batch_*_SAMPLE.parquet)")
    parser.add_argument("--unique_kmers_file", required=True,
                        help="File with globally unique k-mers (one per line or tab-delimited)")
    parser.add_argument("--output_file", required=True,
                        help="Base name for saved .npz files (e.g., .../final_feature_matrix.csv)")
    parser.add_argument("--final_output_folder", required=True,
                        help="Folder to store fold CSVs and the full-data CSV")
    parser.add_argument("--top_k", type=int, default=1500,
                        help="Number of top k-mers to select per fold (default: 1500)")
    parser.add_argument("--n_folds", type=int, default=5,
                        help="Number of cross-validation folds (default: 5)")
    # NEW: optional marker file for Snakemake
    parser.add_argument("--marker_file", default=None,
                        help="Optional marker file to write on success (for Snakemake)")

    args = parser.parse_args()

    try:
        # 1) Build the full feature matrix once
        feature_matrix, labels, sample_ids, unique_kmers_index = create_feature_matrix_from_parquet(
            args.parquet_folder, args.unique_kmers_file, args.output_file
        )

        # 2) Ensure final_output_folder exists (for fold CSVs and the full-data CSV)
        os.makedirs(args.final_output_folder, exist_ok=True)

        # 3) Split into folds, do per-fold χ² selection, and write CSVs
        split_data_and_process(
            feature_matrix=feature_matrix,
            labels=labels,
            sample_ids=sample_ids,
            output_folder=args.final_output_folder,
            unique_kmers_index=unique_kmers_index,
            top_k=args.top_k,
            n_folds=args.n_folds
        )

        # 4) Compute χ² on all internal samples to build the final training set
        print("\n[INFO] Computing χ² on all internal samples to pick final top_k features...")
        chi_scores_full, _ = chi2(feature_matrix, labels)
        order_full = np.argsort(chi_scores_full)
        top_k_full = int(min(args.top_k, feature_matrix.shape[1]))
        top_feature_indices_full = np.sort(order_full[-top_k_full:])

        reverse_index = {v: k for v, k in enumerate(unique_kmers_index.keys())}  # not used; fix properly
        # Rebuild reverse_index correctly:
        reverse_index = {v: k for k, v in unique_kmers_index.items()}
        final_kmer_names = [reverse_index[i] for i in top_feature_indices_full]

        X_full_selected = feature_matrix[:, top_feature_indices_full].toarray()
        full_df = pd.DataFrame(X_full_selected, columns=final_kmer_names)
        full_df["Label"] = labels
        full_df["Sample_ID"] = sample_ids

        full_output_file = os.path.join(args.final_output_folder, "all_internal_topk.csv")
        full_df.to_csv(full_output_file, index=False)
        print(f"[INFO] Saved final full-data top_{top_k_full} features to {full_output_file}")

        # 5) Touch optional marker file for Snakemake
        if args.marker_file:
            out_dir = os.path.dirname(args.marker_file)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(args.marker_file, "w") as f:
                f.write("done\n")
            print(f"[INFO] Wrote marker file: {args.marker_file}")

    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(2)

if __name__ == "__main__":
    main()
