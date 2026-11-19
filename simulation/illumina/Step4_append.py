import os
import argparse

def concatenate_pairs(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Traverse each subfolder (one per sample)
    for sample_folder in os.listdir(input_dir):
        sample_path = os.path.join(input_dir, sample_folder)
        if not os.path.isdir(sample_path):
            continue

        # Look for R1 and R2 inside the subfolder
        r1 = None
        r2 = None
        for file in os.listdir(sample_path):
            if file.endswith("_R1.fastq"):
                r1 = os.path.join(sample_path, file)
            elif file.endswith("_R2.fastq"):
                r2 = os.path.join(sample_path, file)

        # Skip if either file is missing
        if not r1 or not r2:
            print(f"[WARNING] Skipping {sample_folder} — R1 or R2 missing.")
            continue

        # Output file path
        output_file = os.path.join(output_dir, f"{sample_folder}_concatenated.fastq")

        # Concatenate R1 and R2
        with open(output_file, 'wb') as wfd:
            for f in [r1, r2]:
                with open(f, 'rb') as fd:
                    wfd.writelines(fd)

        print(f"[INFO] Appended: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--marker_file", required=True)
    args = parser.parse_args()

    concatenate_pairs(args.input_dir, args.output_dir)

    with open(args.marker_file, "w") as f:
        f.write("done\n")
