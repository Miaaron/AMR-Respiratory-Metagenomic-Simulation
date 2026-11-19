#!/usr/bin/env python3
import subprocess
import os
import glob
import argparse
import random
import json
from concurrent.futures import ProcessPoolExecutor

def get_simulation_number(folder_name):
    return int(folder_name.rsplit('_', 1)[-1][3:])  # Assumes suffix like _sim3

def get_reference_name(genome_file_path):
    return os.path.basename(genome_file_path).split('.')[0]

def remove_header(coverage_file_path):
    with open(coverage_file_path, 'r') as file:
        lines = file.readlines()
    with open(coverage_file_path, 'w') as file:
        file.writelines(lines[1:])

def simulate_one(folder_name, input_base, output_base):
    folder_path = os.path.join(input_base, folder_name)
    simulation_output_path = os.path.join(output_base, folder_name)
    os.makedirs(simulation_output_path, exist_ok=True)

    coverage_file_path = glob.glob(os.path.join(folder_path, "*.txt"))[0]
    remove_header(coverage_file_path)
    genome_file_paths = glob.glob(os.path.join(folder_path, "*.fasta"))

    seed_record = {}

    for genome_file_path in genome_file_paths:
        reference_name = get_reference_name(genome_file_path)
        expected_output = os.path.join(simulation_output_path, f"*_{reference_name}.fastq.gz")
        existing_files = glob.glob(expected_output)

        if not existing_files:
            seed_value = random.randint(0, 99999)
            seed_record[f"{folder_name}_{reference_name}"] = seed_value

            output_file_base = os.path.join(simulation_output_path, f"{folder_name}_{reference_name}")
            command = (
                f"iss generate --coverage_file {coverage_file_path} --genomes {genome_file_path} "
                f"--model NovaSeq --seed {seed_value} --output {output_file_base}"
            )
            print("Executing command:", command)
            process = subprocess.run(command, shell=True, capture_output=True, text=True)

            if process.stderr:
                print(f"Errors for {folder_name} {reference_name}: {process.stderr}")
            else:
                print(f"Output for {folder_name} {reference_name}: {process.stdout}")
        else:
            print(f"Output for {folder_name} {reference_name} already exists. Skipping.")

    return seed_record

def main(input_base, output_base, marker_file):
    folder_names = [f for f in os.listdir(input_base) if os.path.isdir(os.path.join(input_base, f))]
    comprehensive_seed_record = {}

    with ProcessPoolExecutor(max_workers=50) as executor:
        futures = [executor.submit(simulate_one, folder, input_base, output_base) for folder in folder_names]
        for future in futures:
            record = future.result()
            comprehensive_seed_record.update(record)

    with open('simulation_seed_record.json', 'w') as file:
        json.dump(comprehensive_seed_record, file, indent=4)

    with open(marker_file, "w") as f:
        f.write("done\n")

    print("Seed record saved to simulation_seed_record.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallel ISS read simulation")
    parser.add_argument("--input_base", required=True, help="Folder containing populated ISS directories")
    parser.add_argument("--output_base", required=True, help="Output directory for ISS simulated reads")
    parser.add_argument("--marker_file", required=True, help="File to mark completion")
    args = parser.parse_args()

    main(args.input_base, args.output_base, args.marker_file)
