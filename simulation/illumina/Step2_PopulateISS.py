import os
import shutil
import argparse

def copy_files(src_base_dir, dest_base_dir):
    for src_dir_name in os.listdir(src_base_dir):
        src_parts = src_dir_name.split('_')
        src_sample_id = src_parts[0]
        src_resistance_tag = src_parts[1] if len(src_parts) > 1 else None

        for dest_dir_name in os.listdir(dest_base_dir):
            dest_parts = dest_dir_name.split('_')
            dest_sample_id = dest_parts[0]
            dest_resistance_tag = dest_parts[1] if len(dest_parts) > 1 else None

            if src_sample_id == dest_sample_id and src_resistance_tag == dest_resistance_tag:
                src_dir_path = os.path.join(src_base_dir, src_dir_name)
                dest_dir_path = os.path.join(dest_base_dir, dest_dir_name)
                for file_name in os.listdir(src_dir_path):
                    src_file_path = os.path.join(src_dir_path, file_name)
                    dest_file_path = os.path.join(dest_dir_path, file_name)
                    if not os.path.exists(dest_file_path):
                        shutil.copy2(src_file_path, dest_file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--src_dir", required=True)
    parser.add_argument("--dest_dir", required=True)
    parser.add_argument("--marker_file", required=True)
    args = parser.parse_args()
    copy_files(args.src_dir, args.dest_dir)
    with open(args.marker_file, "w") as f:
        f.write("done\n")
