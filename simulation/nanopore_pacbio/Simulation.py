#!/usr/bin/env python3
import os
import glob
import subprocess
import re
import random
import argparse
from concurrent.futures import ProcessPoolExecutor

def run_pbsim2(fasta_file, coverage, output_dir, hmm_model_path, sim_number, seed):
    # Build output prefix; PBSIM will add its own suffixes.
    prefix = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(fasta_file))[0]}_sim{sim_number}")
    pbsim_command = [
        '/PATH/TO/ENV/bin/pbsim2/bin/pbsim',
        '--depth', str(coverage),
        '--hmm_model', hmm_model_path,
        '--prefix', prefix,
        '--seed', str(seed),
        '--difference-ratio', '23:31:46',
        fasta_file  # PBSIM requires a reference FASTA as the last argument.
    ]
    result = subprocess.run(pbsim_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout, result.stderr

def calculate_expected_bases(coverage, genome_length):
    """Calculate expected number of bases given coverage and genome length."""
    return coverage * genome_length

def calculate_total_bases(fastq_filename):
    """Calculate the total number of bases in a FASTQ file."""
    total_bases = 0
    with open(fastq_filename, 'r') as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 1:
                total_bases += len(line.strip())
    return total_bases

def simulate_sample(sample_dir, output_dir, hmm_model_path, total_bases, pathogen_percentage, sim_count_resistant, sim_count_susceptible):
    # Regular expression to identify numeric filenames (used to mark the pathogen file)
    numeric_filename_pattern = re.compile(r'^\d+')
    pathogen_file = None
    genome_lengths = {}

    # Iterate over all FASTA files in the sample directory
    for fasta_file in glob.glob(os.path.join(sample_dir, '*.fasta')):
        # Assume the pathogen file has a filename starting with a number.
        if re.match(numeric_filename_pattern, os.path.basename(fasta_file).split('.')[0]):
            pathogen_file = fasta_file
        with open(fasta_file, 'r') as f:
            seq = ''.join([line.strip() for line in f if not line.startswith('>')])
            genome_lengths[fasta_file] = len(seq)

    if pathogen_file is None:
        raise ValueError(f"Pathogen file not found in {sample_dir}.")

    print("\nGenome lengths:")
    for fasta_file, length in genome_lengths.items():
        print(f"{fasta_file}: {length} bases")

    # Calculate pathogen and commensal bases, then coverage.
    pathogen_bases = total_bases * pathogen_percentage
    commensal_bases = total_bases - pathogen_bases
    pathogen_coverage = pathogen_bases / genome_lengths[pathogen_file]
    print(f"\nPathogen file: {pathogen_file}")
    print(f"Pathogen coverage: {pathogen_coverage:.4f}x")
    print(f"Expected bases for pathogen: {calculate_expected_bases(pathogen_coverage, genome_lengths[pathogen_file]):.0f}")

    print(f"\nTotal remaining bases for other species: {commensal_bases} bases")
    other_species_files = [f for f in genome_lengths if f != pathogen_file]
    total_other_species_length = sum(genome_lengths[f] for f in other_species_files)
    even_coverages = {}
    total_allocated_bases = 0
    for fasta_file in other_species_files:
        species_length = genome_lengths[fasta_file]
        allocated_bases = random.uniform(0.9, 1.1) * (commensal_bases * (species_length / total_other_species_length))
        total_allocated_bases += allocated_bases
        even_coverages[fasta_file] = allocated_bases / species_length
        print(f"{fasta_file}: Allocated bases: {allocated_bases:.0f}, Coverage: {even_coverages[fasta_file]:.4f}x")
    scaling_factor = commensal_bases / total_allocated_bases
    for fasta_file in even_coverages:
        even_coverages[fasta_file] *= scaling_factor
        print(f"{fasta_file}: Adjusted coverage: {even_coverages[fasta_file]:.4f}x")
    # Combine coverages for pathogen and commensals.
    coverages = {pathogen_file: pathogen_coverage}
    coverages.update(even_coverages)

    # Determine simulation counts per FASTA file based on filename.
    sample_label = os.path.basename(sample_dir)  # gets "Resistant" or "Susceptible"
    if "Resistant" in sample_label:
        default_sim_count = sim_count_resistant
    elif "Susceptible" in sample_label:
        default_sim_count = sim_count_susceptible
    else:
        default_sim_count = sim_count_susceptible  # fallback

    simulation_counts = {f: default_sim_count for f in genome_lengths}


    print("\nUsing simulation counts per file:")
    for f, count in simulation_counts.items():
        print(f"{os.path.basename(f)}: {count}")

    # Ensure output directory exists.
    os.makedirs(output_dir, exist_ok=True)
    # Run simulations for each FASTA file.
    for fasta_file, cov in coverages.items():
        sim_count = simulation_counts[fasta_file]
        for sim_number in range(1, sim_count + 1):
            seed = sim_number  # Modify this if you want a more complex seed.
            print(f"Submitting simulation for {os.path.basename(fasta_file)}: sim {sim_number} with coverage {cov:.4f}x, seed: {seed}")
            stdout, stderr = run_pbsim2(fasta_file, cov, output_dir, hmm_model_path, sim_number, seed)
            print("Stdout:", stdout)
            print("Stderr:", stderr)

    # Write a marker file to indicate simulation completion.
    marker_file = os.path.join(output_dir, "simulation_done.txt")
    with open(marker_file, "w") as marker:
        marker.write("Simulation completed.\n")
    print("Simulation completed for sample:", sample_dir)

def main():
    parser = argparse.ArgumentParser(description="Simulate a sample using PBSIM2")
    parser.add_argument("--sample_dir", required=True, help="Path to the sample directory containing FASTA files")
    parser.add_argument("--output_dir", required=True, help="Output directory for simulation results (per sample)")
    parser.add_argument("--hmm_model", required=True, help="Path to the PBSIM2 HMM model file")
    parser.add_argument("--total_bases", type=int, required=True, help="Total number of bases (e.g., 1000000000)")
    parser.add_argument("--pathogen_percentage", type=float, required=True, help="Fraction of bases for the pathogen (e.g., 0.05)")
    parser.add_argument("--simulation_count_resistant", type=int, required=True, help="Number of simulations for Resistant FASTA files")
    parser.add_argument("--simulation_count_susceptible", type=int, required=True, help="Number of simulations for Susceptible FASTA files")
    args = parser.parse_args()

    simulate_sample(
        args.sample_dir,
        args.output_dir,
        args.hmm_model,
        args.total_bases,
        args.pathogen_percentage,
        args.simulation_count_resistant,
        args.simulation_count_susceptible
    )

if __name__ == "__main__":
    main()
