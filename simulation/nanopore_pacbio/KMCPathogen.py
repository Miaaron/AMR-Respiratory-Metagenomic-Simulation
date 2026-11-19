#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse

def run_kmc(kmc_path, input_file, output_db, temp_dir, kmer_size=13):
    """
    Run KMC to count k-mers in a given FASTA file.
    """
    command = [
        kmc_path,
        f"-k{kmer_size}",
        "-ci1",  # Include all k-mers with count >= 1
        "-fm",   # Specify FASTA format
        input_file,
        output_db,
        temp_dir
    ]
    try:
        print(f"Running KMC on {input_file}...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"KMC completed successfully for {input_file}: {result.stdout.decode('utf-8')}")
    except subprocess.CalledProcessError as e:
        print(f"Error running KMC for {input_file}: {e.stderr.decode('utf-8')}")
        return False
    return True

def extract_unique_kmers(kmc_dump_path, kmc_db, temp_kmer_txt):
    """
    Extract k-mers from a KMC database and save them to a temporary text file.
    """
    command = [kmc_dump_path, kmc_db, temp_kmer_txt]
    try:
        print(f"Extracting k-mers to {temp_kmer_txt}...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"kmc_dump completed successfully for {temp_kmer_txt}: {result.stdout.decode('utf-8')}")
    except subprocess.CalledProcessError as e:
        print(f"Error running kmc_dump for {temp_kmer_txt}: {e.stderr.decode('utf-8')}")
        return False
    return True

def process_fasta_file(kmc_path, kmc_dump_path, fasta_file, temp_dir, global_unique_kmers, kmer_size=13):
    """
    Process a single FASTA file to extract unique k-mers and add them to the global set.
    """
    sample_name = os.path.basename(fasta_file).replace(".fasta", "")
    kmc_output_db = os.path.join(temp_dir, f"{sample_name}_kmc_db")
    temp_kmer_txt = os.path.join(temp_dir, f"{sample_name}_kmers.txt")

    if run_kmc(kmc_path, fasta_file, kmc_output_db, temp_dir, kmer_size):
        if extract_unique_kmers(kmc_dump_path, kmc_output_db, temp_kmer_txt):
            with open(temp_kmer_txt, "r") as temp_file:
                for line in temp_file:
                    kmer = line.split()[0]  # Take the first column (the k-mer)
                    global_unique_kmers.add(kmer)
            print(f"Added k-mers from {temp_kmer_txt} to global unique set")
        else:
            print(f"Failed to extract k-mers from {kmc_output_db}")
    else:
        print(f"Failed to run KMC on {fasta_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract unique k-mers from FASTA files using KMC")
    parser.add_argument("--kmc", required=True, help="Path to the kmc executable")
    parser.add_argument("--kmc_dump", required=True, help="Path to the kmc_dump executable")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTA files")
    parser.add_argument("--output_file", required=True, help="Output file for the unique k-mers")
    parser.add_argument("--temp_dir", required=True, help="Temporary directory for intermediate kmc files")
    parser.add_argument("--kmer_size", type=int, default=13, help="K-mer size (default: 13)")
    args = parser.parse_args()

    # Ensure temp and output directories exist
    os.makedirs(args.temp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    # Get list of FASTA files
    fasta_files = glob.glob(os.path.join(args.input_dir, "*.fasta"))
    global_unique_kmers = set()

    for fasta_file in fasta_files:
        process_fasta_file(args.kmc, args.kmc_dump, fasta_file, args.temp_dir, global_unique_kmers, args.kmer_size)

    # Write the globally unique k-mers to the final output file
    with open(args.output_file, "w") as output_file:
        for kmer in global_unique_kmers:
            output_file.write(f"{kmer}\t1\n")
    print(f"KMC processing completed. Total unique k-mers: {len(global_unique_kmers)}")

if __name__ == "__main__":
    main()
