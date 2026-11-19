import os
import glob
import argparse
import subprocess
import pandas as pd


def run_kmc(kmc_path, input_file, output_db, temp_dir, kmer_size=13):
    command = [
        kmc_path, f"-k{kmer_size}", "-ci1", "-fq",
        input_file, output_db, temp_dir
    ]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running KMC for {input_file}: {e.stderr.decode()}")
        return False
    return True


def extract_unique_kmers(kmc_dump_path, kmc_db, temp_kmer_txt):
    command = [kmc_dump_path, kmc_db, temp_kmer_txt]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        unique_kmers = set()
        with open(temp_kmer_txt, "r") as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    unique_kmers.add(parts[0])
        return unique_kmers
    except subprocess.CalledProcessError as e:
        print(f"Error running kmc_dump for {kmc_db}: {e.stderr.decode()}")
        return None


def process_sample(sample_name, fastq_file, kmc_path, kmc_dump_path, temp_dir,
                   pathogen_kmers, batch_size, output_dir, kmer_size=13):
    db = os.path.join(temp_dir, f"{sample_name}_kmc_db")
    txt = os.path.join(temp_dir, f"{sample_name}_kmers.txt")

    ok = run_kmc(kmc_path, fastq_file, db, temp_dir, kmer_size)
    if not ok:
        return False

    sample_kmers = extract_unique_kmers(kmc_dump_path, db, txt)
    if sample_kmers is None:
        return False

    for i in range(0, len(pathogen_kmers), batch_size):
        batch = pathogen_kmers[i:i + batch_size]
        presence = [1 if kmer in sample_kmers else 0 for kmer in batch]
        df = pd.DataFrame({"kmer": batch, "presence": presence})
        out_file = os.path.join(output_dir, f"batch_{i // batch_size + 1}_{sample_name}.parquet")
        df.to_parquet(out_file, index=False)

    with open(os.path.join(output_dir, f"{sample_name}_batch_done.txt"), "w") as f:
        f.write("done\n")

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--pathogen_kmers_file", required=True)
    parser.add_argument("--kmc", required=True)
    parser.add_argument("--kmc_dump", required=True)
    parser.add_argument("--temp_dir", required=True)
    parser.add_argument("--batch_size", type=int, default=2000000)
    parser.add_argument("--kmer_size", type=int, default=13)
    parser.add_argument("--marker_file", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.temp_dir, exist_ok=True)

    with open(args.pathogen_kmers_file) as f:
        pathogen_kmers = [line.split()[0] for line in f if line.strip()]

    fastq_files = glob.glob(os.path.join(args.input_dir, "*.fastq"))
    if not fastq_files:
        raise RuntimeError(f"No FASTQ files found in input_dir: {args.input_dir}")

    failed_samples = []

    for fq in fastq_files:
        sample_name = os.path.basename(fq).replace(".fastq", "")
        success = process_sample(
            sample_name=sample_name,
            fastq_file=fq,
            kmc_path=args.kmc,
            kmc_dump_path=args.kmc_dump,
            temp_dir=args.temp_dir,
            pathogen_kmers=pathogen_kmers,
            batch_size=args.batch_size,
            output_dir=args.output_dir,
            kmer_size=args.kmer_size
        )
        if not success:
            failed_samples.append(sample_name)

    if failed_samples:
        raise RuntimeError(
            f"Step5_Batches failed for {len(failed_samples)} sample(s): "
            + ", ".join(failed_samples[:10])
            + (" ..." if len(failed_samples) > 10 else "")
        )

    with open(args.marker_file, "w") as f:
        f.write("done\n")