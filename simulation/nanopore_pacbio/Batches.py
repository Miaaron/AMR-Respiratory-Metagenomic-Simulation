#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse
from scipy.sparse import csr_matrix, save_npz
import numpy as np
import pandas as pd

def run_kmc(kmc_path, input_file, output_db, temp_dir, kmer_size=13):
    """
    Run KMC on a single FASTQ file.
    """
    command = [
        kmc_path,
        f"-k{kmer_size}",
        "-ci1",  # Include k-mers that occur at least once
        "-fq",   # Input file is in FASTQ format
        input_file,
        output_db,
        temp_dir,
    ]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running KMC for {input_file}: {e.stderr.decode('utf-8')}")
        return False
    return True

def extract_unique_kmers(kmc_dump_path, kmc_db, temp_kmer_txt):
    """
    Extract unique k-mers from a KMC database.
    """
    command = [kmc_dump_path, kmc_db, temp_kmer_txt]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        unique_kmers = set()
        with open(temp_kmer_txt, 'r') as temp_file:
            for line in temp_file:
                # Each line is expected to have a k-mer and its count; we take only the k-mer.
                kmer, _ = line.strip().split()
                unique_kmers.add(kmer)
        return unique_kmers
    except subprocess.CalledProcessError as e:
        print(f"Error running kmc_dump for {kmc_db}: {e.stderr.decode('utf-8')}")
        return set()

def process_sample(sample_name, fastq_file, kmc_path, kmc_dump_path, temp_dir, pathogen_kmers, batch_size, output_dir, kmer_size=13):
    """
    Process a single FASTQ file to generate presence/absence data in batches and save as Parquet files.
    """
    kmc_output_db = os.path.join(temp_dir, f"{sample_name}_kmc_db")
    temp_kmer_txt = os.path.join(temp_dir, f"{sample_name}_kmers.txt")
    
    if run_kmc(kmc_path, fastq_file, kmc_output_db, temp_dir, kmer_size):
        sample_kmers = extract_unique_kmers(kmc_dump_path, kmc_output_db, temp_kmer_txt)
        if sample_kmers:
            for i in range(0, len(pathogen_kmers), batch_size):
                batch = pathogen_kmers[i:i + batch_size]
                print(f"Processing batch {i // batch_size + 1} for sample {sample_name}")
                try:
                    presence_absence = [1 if kmer in sample_kmers else 0 for kmer in batch]
                    batch_df = pd.DataFrame({
                        'kmer': batch,
                        'presence': presence_absence
                    })
                    output_file = os.path.join(output_dir, f"batch_{i // batch_size + 1}_{sample_name}.parquet")
                    batch_df.to_parquet(output_file, index=False)
                    print(f"Saved Parquet for {sample_name}, batch {i // batch_size + 1} to {output_file}")
                except Exception as e:
                    print(f"Error processing batch {i // batch_size + 1} for {sample_name}: {e}")
        else:
            print(f"No k-mers extracted from {kmc_output_db}")
    else:
        print(f"Failed to run KMC on {fastq_file}")
    # Write done marker
    marker_file = os.path.join(output_dir, f"{sample_name}_batch_done.txt")
    with open(marker_file, "w") as f:
        f.write("done\n")
    print(f"Batching complete for {sample_name}, marker written to {marker_file}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create batches (Parquet files) from a merged FASTQ file using KMC")
    parser.add_argument("--input", required=True, help="Input merged FASTQ file for a sample")
    parser.add_argument("--sample_name", required=True, help="Sample name (used to name output files)")
    parser.add_argument("--output_dir", required=True, help="Output directory where Parquet files will be saved")
    parser.add_argument("--kmc", required=True, help="Path to the kmc executable")
    parser.add_argument("--kmc_dump", required=True, help="Path to the kmc_dump executable")
    parser.add_argument("--temp_dir", required=True, help="Temporary directory for kmc intermediate files")
    parser.add_argument("--pathogen_kmers_file", required=True, help="File containing pathogen k-mers (unique)")
    parser.add_argument("--batch_size", type=int, default=2000000, help="Batch size for splitting k-mers")
    parser.add_argument("--kmer_size", type=int, default=13, help="k-mer size (default: 13)")
    args = parser.parse_args()

    # Load pathogen k-mers from file
    with open(args.pathogen_kmers_file, 'r') as f:
        pathogen_kmers = [line.strip().split()[0] for line in f]

    os.makedirs(args.output_dir, exist_ok=True)
    
    process_sample(args.sample_name, args.input, args.kmc, args.kmc_dump, args.temp_dir, pathogen_kmers, args.batch_size, args.output_dir, args.kmer_size)
