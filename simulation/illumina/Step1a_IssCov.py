import os
import re
from glob import glob
from collections import defaultdict

def extract_chunks_by_simulation(log_path):
    with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    chunks = re.split(r'(?:^|\n)Stderr: :::: Simulation parameters :::', content)
    return [chunk.strip() for chunk in chunks if 'prefix :' in chunk]

def extract_sim_metadata(chunk):
    match = re.search(r'prefix\s*:\s*.*/(sample\d+_[^/]+)/[^/]+_sim(\d+)', chunk)
    if match:
        sample_prefix = match.group(1)
        sim_num = match.group(2)
        return f"{sample_prefix}_sim{sim_num}"
    else:
        print("[⚠] Skipped chunk — no prefix match.")
        return None

def extract_coverage_data(chunk):
    ref_map = dict(re.findall(r'ref\.(\d+)\s+\(len:\d+\)\s+:\s+([^\s]+)', chunk))
    depths = re.findall(r'Simulation stats \(ref\.(\d+)\).*?depth\s*:\s*([\d.]+)', chunk, re.DOTALL)
    return [(ref_map[rid], depth) for rid, depth in depths if rid in ref_map]

def write_coverage(output_dir, sim_name, coverage_data):
    sim_folder = os.path.join(output_dir, sim_name)
    os.makedirs(sim_folder, exist_ok=True)
    outpath = os.path.join(sim_folder, "coverage.txt")
    with open(outpath, 'w') as f:
        f.write("Reference\tCoverage\n")
        for ref, cov in coverage_data:
            f.write(f"{ref}\t{cov}\n")
    print(f"[✓] Written: {outpath}")

def process_all_logs(logs_dir, output_dir):
    all_logs = glob(os.path.join(logs_dir, "*.log"))
    sim_dict = defaultdict(list)

    for log_path in all_logs:
        chunks = extract_chunks_by_simulation(log_path)
        for chunk in chunks:
            sim_name = extract_sim_metadata(chunk)
            if sim_name:
                coverage = extract_coverage_data(chunk)
                sim_dict[sim_name].extend(coverage)

    for sim_name, entries in sim_dict.items():
        write_coverage(output_dir, sim_name, entries)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--logs_dir", required=True, help="Directory with simulation logs")
    parser.add_argument("--output_dir", required=True, help="Where to write coverage.txt files")
    parser.add_argument("--marker_file", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    process_all_logs(args.logs_dir, args.output_dir)

    with open(args.marker_file, "w") as f:
        f.write("done\n")
