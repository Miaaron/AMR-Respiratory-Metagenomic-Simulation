#!/usr/bin/env python3
import os
import shutil
import argparse

def require_file(path, desc):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required {desc}: {path}")

def wipe_dir_contents(path, keep_files=None):
    keep_files = set(keep_files or [])
    if not os.path.isdir(path):
        print(f"Skip: {path} is not a directory")
        return

    for name in os.listdir(path):
        if name in keep_files:
            print(f"Keeping: {os.path.join(path, name)}")
            continue

        full = os.path.join(path, name)
        try:
            if os.path.isdir(full) and not os.path.islink(full):
                shutil.rmtree(full)
                print(f"Deleted directory: {full}")
            else:
                os.remove(full)
                print(f"Deleted file: {full}")
        except Exception as e:
            print(f"Could not delete {full}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Final cleanup after pipeline completion")
    parser.add_argument("--sim_dir", required=True)
    parser.add_argument("--batch_dir", required=True)
    parser.add_argument("--temp_dir", required=True)
    parser.add_argument("--merge_dir", required=True)
    parser.add_argument("--final_dir", required=True)
    parser.add_argument("--kmc_pathogen_dir", required=False)
    parser.add_argument("--keep_unique_kmers", action="store_true")
    parser.add_argument("--cleanup_marker", required=True)
    args = parser.parse_args()

    # Strong safety checks
    require_file(os.path.join(args.final_dir, "fold_1_train.csv"), "final fold_1_train.csv")
    require_file(os.path.join(args.final_dir, "fold_1_test.csv"), "final fold_1_test.csv")
    require_file(os.path.join(args.final_dir, "all_internal_topk.csv"), "final all_internal_topk.csv")

    print("Final outputs found. Starting cleanup.")
    print(f"Keeping merge directory untouched: {args.merge_dir}")
    print(f"Keeping final directory untouched: {args.final_dir}")

    # Delete simulation intermediates
    wipe_dir_contents(args.sim_dir)

    # Delete batching intermediates (the big parquet files)
    wipe_dir_contents(args.batch_dir)

    # Delete KMC temp/intermediate files used during batching
    wipe_dir_contents(args.temp_dir)

    # Optionally keep only unique_kmers.txt in KMCPathogens
    if args.kmc_pathogen_dir:
        if args.keep_unique_kmers:
            wipe_dir_contents(args.kmc_pathogen_dir, keep_files={"unique_kmers.txt"})
        else:
            wipe_dir_contents(args.kmc_pathogen_dir)

    # Write marker
    with open(args.cleanup_marker, "w") as f:
        f.write("cleanup_done\n")

    print(f"Cleanup complete. Marker written to {args.cleanup_marker}")

if __name__ == "__main__":
    main()