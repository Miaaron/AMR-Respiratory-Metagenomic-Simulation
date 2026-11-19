#!/usr/bin/env python3
import os
import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack, save_npz
from sklearn.feature_selection import chi2
from sklearn.model_selection import GroupKFold
import argparse

#############################################
# 1. HELPER FUNCTIONS
#############################################
def load_unique_kmers(unique_kmers_file):
    """Load and clean unique k-mers from a text file (one per line)."""
    with open(unique_kmers_file, 'r') as f:
        unique_kmers_list = [line.split('\t')[0].strip().upper() for line in f]
    return {kmer: idx for idx, kmer in enumerate(unique_kmers_list)}

def extract_label_and_sample_id_from_parquet(parquet_filename):
    """
    Extract label and sample_id from a Parquet filename.
    Expected filename format: "batch_{batch_num}_{sample_id}.parquet"
    The sample_id is expected to contain "Resistant" or "Susceptible".
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
        label = -1
    return label, sample_id

#############################################
# 2. BUILD FEATURE MATRIX FROM PARQUET
#############################################
def create_feature_matrix_from_parquet(parquet_folder, unique_kmers_file, output_file):
    """
    Reads all Parquet files from parquet_folder (each with columns 'kmer' and 'presence'),
    aggregates data per sample (using the sample id extracted from the filename),
    builds a sparse presence/absence feature matrix (columns correspond to unique kmers),
    and saves the matrix along with labels and sample IDs.
    """
    unique_kmers_index = load_unique_kmers(unique_kmers_file)
    feature_matrix_dict = {}  # sample_id -> sparse row vector
    labels_dict = {}          # sample_id -> label
    sample_ids = []           # to preserve order

    parquet_files = glob.glob(os.path.join(parquet_folder, "*.parquet"))
    if not parquet_files:
        raise ValueError(f"No Parquet files found in {parquet_folder}")

    for pf in parquet_files:
        try:
            label, sample_id = extract_label_and_sample_id_from_parquet(pf)
        except Exception as e:
            print(f"Error parsing filename {pf}: {e}")
            continue

        try:
            df = pd.read_parquet(pf)
        except Exception as e:
            print(f"Error reading {pf}: {e}")
            continue

        df["kmer"] = df["kmer"].str.strip().str.upper()
        mask = df["kmer"].isin(unique_kmers_index)
        df = df.loc[mask]
        if df.empty:
            continue

        indices = df["kmer"].map(unique_kmers_index).values
        values = df["presence"].values
        row = csr_matrix((values, ([0]*len(indices), indices)),
                         shape=(1, len(unique_kmers_index)))

        if sample_id in feature_matrix_dict:
            feature_matrix_dict[sample_id] += row
        else:
            feature_matrix_dict[sample_id] = row
            labels_dict[sample_id] = label
            sample_ids.append(sample_id)

    # Convert counts to binary (presence/absence)
    for sid, row in feature_matrix_dict.items():
        if row.nnz > 0:
            row.data[:] = 1

    # Stack rows in the order of sample_ids
    feature_matrix = vstack([feature_matrix_dict[sid] for sid in sample_ids])
    labels = np.array([labels_dict[sid] for sid in sample_ids])
    sample_ids_array = np.array(sample_ids)

    # Save the full feature matrix and labels
    npz_file = output_file.replace(".csv", ".npz")
    labels_file = output_file.replace(".csv", "_labels.npz")
    save_npz(npz_file, feature_matrix)
    np.savez(labels_file, labels=labels, sample_ids=sample_ids_array)
    print(f"Feature matrix saved to {npz_file}")
    print(f"Labels and sample IDs saved to {labels_file}")

    return feature_matrix, labels, sample_ids_array

#############################################
# 3. CROSS-VALIDATION SPLIT & PER-FOLD SELECTION
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
    For each fold:
      1. Split into train/test.
      2. Compute χ²–scores on the TRAIN set only, pick top_k k-mers.
      3. Build train/test DataFrames with those top_k columns.
      4. Save CSVs for each fold.
    """
    reverse_index = {v: k for k, v in unique_kmers_index.items()}
    sample_names = np.array([sid.rsplit('_', 1)[0] for sid in sample_ids])
    gkf = GroupKFold(n_splits=n_folds)

    fold_number = 1
    for train_idx, test_idx in gkf.split(feature_matrix, labels, groups=sample_names):
        X_train_full = feature_matrix[train_idx]
        y_train = labels[train_idx]
        X_test_full = feature_matrix[test_idx]
        y_test = labels[test_idx]

        chi_scores, _ = chi2(X_train_full, y_train)
        top_feature_indices = np.argsort(chi_scores)[-top_k:]
        top_feature_indices = np.sort(top_feature_indices)

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
        print(f"Fold {fold_number}: saved {train_file} (train) and {test_file} (test).")
        fold_number += 1

#############################################
# 4. MAIN PIPELINE
#############################################
def main():
    parser = argparse.ArgumentParser(
        description="Build feature matrix from Parquet files and perform CV split with per-fold χ² selection."
    )
    parser.add_argument(
        "--parquet_folder", required=True,
        help="Directory containing Parquet files"
    )
    parser.add_argument(
        "--unique_kmers_file", required=True,
        help="File with globally unique k-mers (one per line or tab-delimited)"
    )
    parser.add_argument(
        "--output_file", required=True,
        help="Base name for saved .npz files (e.g., feature_matrix.csv)"
    )
    parser.add_argument(
        "--final_output_folder", required=True,
        help="Folder to store final CSV outputs (folds) – also used for final full-data output"
    )
    parser.add_argument(
        "--top_k", type=int, default=1500,
        help="Number of top k-mers to select per fold (default: 1500)."
    )
    parser.add_argument(
        "--n_folds", type=int, default=5,
        help="Number of cross-validation folds (default: 5)."
    )
    args = parser.parse_args()

    # Ensure the parent directory of output_file exists
    parent_dir = os.path.dirname(args.output_file)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)

    # 1) Build the full feature matrix once
    feature_matrix, labels, sample_ids = create_feature_matrix_from_parquet(
        args.parquet_folder,
        args.unique_kmers_file,
        args.output_file
    )

    # 2) Load unique_kmers_index (to recover k-mer names later)
    unique_kmers_index = load_unique_kmers(args.unique_kmers_file)

    # 3) Ensure final_output_folder exists (for all fold CSVs and the full-data CSV)
    os.makedirs(args.final_output_folder, exist_ok=True)

    # 4) Split into folds, do per-fold χ² selection, and write CSVs
    split_data_and_process(
        feature_matrix=feature_matrix,
        labels=labels,
        sample_ids=sample_ids,
        output_folder=args.final_output_folder,
        unique_kmers_index=unique_kmers_index,
        top_k=args.top_k,
        n_folds=args.n_folds
    )

    # 5) Compute χ² on 100% of internal data to build the final training set
    print("\nComputing χ² on all internal samples to pick final top_k features...")
    chi_scores_full, _ = chi2(feature_matrix, labels)
    top_feature_indices_full = np.argsort(chi_scores_full)[-args.top_k:]
    top_feature_indices_full = np.sort(top_feature_indices_full)

    # Reverse index to map column indices back to k-mer strings
    reverse_index = {v: k for k, v in unique_kmers_index.items()}
    final_kmer_names = [reverse_index[i] for i in top_feature_indices_full]

    # Build a dense matrix of all internal samples on those top_k features
    X_full_selected = feature_matrix[:, top_feature_indices_full].toarray()
    full_df = pd.DataFrame(X_full_selected, columns=final_kmer_names)
    full_df["Label"] = labels
    full_df["Sample_ID"] = sample_ids

    # Save the full-data CSV
    full_output_file = os.path.join(args.final_output_folder, "all_internal_topk.csv")
    full_df.to_csv(full_output_file, index=False)
    print(f"Saved final, full-dataset (100%) top_{args.top_k} features to {full_output_file}")

if __name__ == "__main__":
    main()
