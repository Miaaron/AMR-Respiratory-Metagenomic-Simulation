import os
import glob
import argparse

def merge_fastq_files(base_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    subdirs = [f.path for f in os.scandir(base_dir) if f.is_dir()]

    for subdir in subdirs:
        folder_name = os.path.basename(subdir)
        sample_output_dir = os.path.join(output_dir, folder_name)
        os.makedirs(sample_output_dir, exist_ok=True)

        r1_files = sorted(glob.glob(os.path.join(subdir, '*_R1.fastq')))
        r2_files = sorted(glob.glob(os.path.join(subdir, '*_R2.fastq')))

        merged_r1_file = os.path.join(sample_output_dir, f"{folder_name}_R1.fastq")
        merged_r2_file = os.path.join(sample_output_dir, f"{folder_name}_R2.fastq")

        with open(merged_r1_file, 'wb') as merged_r1:
            for r1_file in r1_files:
                with open(r1_file, 'rb') as f:
                    merged_r1.write(f.read())

        with open(merged_r2_file, 'wb') as merged_r2:
            for r2_file in r2_files:
                with open(r2_file, 'rb') as f:
                    merged_r2.write(f.read())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--marker_file", required=True)
    args = parser.parse_args()
    merge_fastq_files(args.input_dir, args.output_dir)
    with open(args.marker_file, "w") as f:
        f.write("done\n")
