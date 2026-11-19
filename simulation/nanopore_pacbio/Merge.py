# === Updated Merge.py ===
import os
import glob
import re
import argparse

def modify_fastq_header(line, sim_number, species_name):
    if line.startswith('@'):
        parts = line.split(maxsplit=1)
        modified_identifier = f"{parts[0]}_{species_name}_sim{sim_number}"
        modified_line = f"{modified_identifier} {' '.join(parts[1:])}\n" if len(parts) > 1 else f"{modified_identifier}\n"
        return modified_line
    return line

def count_bases_in_fastq(file_path):
    total_bases = 0
    total_reads = 0
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                total_bases += len(line.strip())
                total_reads += 1
    return total_bases, total_reads

def merge_fastq_files(sample_dir, target_directory):
    sample_name = os.path.basename(os.path.normpath(sample_dir))
    os.makedirs(target_directory, exist_ok=True)
    fastq_files = glob.glob(os.path.join(sample_dir, '*.fastq'))

    sim_groups = {}
    pattern = re.compile(r'^(.*?)_sim(\d+)_')
    for file_path in fastq_files:
        file_name = os.path.basename(file_path)
        match = pattern.search(file_name)
        if match:
            species_name = match.group(1)
            sim_number = match.group(2)
            sim_groups.setdefault(sim_number, []).append((file_path, species_name))
        else:
            print(f"Warning: No species or sim number found in file name: {file_name}")

    for sim_number, file_info in sim_groups.items():
        merged_file_name = f"{sample_name}_sim{sim_number}.fastq"
        merged_file_path = os.path.join(target_directory, merged_file_name)

        print(f"\nMerging files for sim{sim_number} into {merged_file_name}.")
        total_bases_before_merge = 0
        total_reads_before_merge = 0

        with open(merged_file_path, 'w') as merged_file:
            for file_path, species_name in file_info:
                bases, reads = count_bases_in_fastq(file_path)
                total_bases_before_merge += bases
                total_reads_before_merge += reads
                with open(file_path, 'r') as f:
                    read_lines = []
                    for line_num, line in enumerate(f):
                        read_lines.append(line)
                        if line_num % 4 == 3:
                            if len(read_lines) == 4:
                                seq_line = read_lines[1].strip()
                                qual_line = read_lines[3].strip()
                                if len(seq_line) != len(qual_line):
                                    print(f"Warning: Length mismatch in {file_path} for read {read_lines[0].strip()}")
                                read_lines[0] = modify_fastq_header(read_lines[0], sim_number, species_name)
                                read_lines[2] = "+" + read_lines[0][1:]
                                merged_file.write("".join(read_lines))
                                read_lines = []

        total_bases_after_merge, total_reads_after_merge = count_bases_in_fastq(merged_file_path)
        print(f"Merged file saved at: {merged_file_path}")
        print(f"Total reads before merge: {total_reads_before_merge}, Bases before merge: {total_bases_before_merge}")
        print(f"Total reads after merge: {total_reads_after_merge}, Bases after merge: {total_bases_after_merge}")

    # Write marker file
    with open(os.path.join(target_directory, f"{sample_name}_merge_done.txt"), 'w') as f:
        f.write("done\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge FASTQ files for a sample based on simulation number")
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    merge_fastq_files(args.input_dir, args.output_dir)