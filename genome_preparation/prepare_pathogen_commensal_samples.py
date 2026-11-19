#!/usr/bin/env python3
import os
import shutil
from pathlib import Path

# ---------- CONFIG ----------
# Common commensal FASTA files (all copied into every sample)
COMMENSAL_PATH = Path("/PATH/TO/DATASET/CommensalS/UnpackedSeqT")

# Root that holds species folders with Resistant*/Susceptible* subfolders
# e.g., .../Respiratory/Nodupes/Klebsiella_Pneumoniae/Resistantamikacin/cleaned_FASTA/*.fa
NODUPES_ROOT = Path("/PATH/TO/DATASET/Processnodupes/Respiratory/Nodupes")

# Where you want the final structure (Nanopore as requested)
# e.g., /PATH/TO/WORKDIR/Samples/ABamikacin/AB/sample1_Resistantamikacin_AcinetobacterBaumannii/
OUTPUT_ROOT = Path("/PATH/TO/WORKDIR/Samples")

# Map species folder name -> (abbr, pretty_name_without_spaces)
SPECIES_MAP = {
    "Acinetobacter_Baumannii": ("AB", "AcinetobacterBaumannii"),
    "Klebsiella_Pneumoniae":   ("KP", "KlebsiellaPneumoniae"),
    "Pseudomonas_Aeruginosa":  ("PA", "PseudomonasAeruginosa"),
}

# File extensions to treat as FASTA
FA_EXTS = {".fa", ".fna", ".fasta", ".fsa", ".fas"}

# Set to True to symlink commensals instead of copying (saves space)
SYMLINK_COMMENSALS = False
# ---------------------------

def list_fasta_files(folder: Path):
    if not folder.exists():
        return []
    return sorted([p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in FA_EXTS])

def safe_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    shutil.copy2(src, dst)

def safe_link(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    os.symlink(src, dst)

def main():
    commensals = list_fasta_files(COMMENSAL_PATH)
    if not commensals:
        raise SystemExit(f"No commensal FASTA files found in {COMMENSAL_PATH}")

    for species_dir_name, (abbr, pretty) in SPECIES_MAP.items():
        species_root = NODUPES_ROOT / species_dir_name
        if not species_root.exists():
            print(f"[skip] Missing species folder: {species_root}")
            continue

        # Iterate Resistant*/Susceptible* antibiotic categories
        for cat_dir in sorted([d for d in species_root.iterdir() if d.is_dir()]):
            cat_name = cat_dir.name  # e.g., Resistantamikacin
            if not (cat_name.startswith("Resistant") or cat_name.startswith("Susceptible")):
                continue

            # Derive antibiotic name and top-level output (e.g., ABamikacin/AB)
            antibiotic = cat_name.replace("Resistant", "").replace("Susceptible", "")
            antibiotic_lc = antibiotic.lower()
            out_top = OUTPUT_ROOT / f"{abbr}{antibiotic_lc}" / abbr
            out_top.mkdir(parents=True, exist_ok=True)

            # Pathogen files live here
            cleaned = cat_dir / "cleaned_FASTA"
            pathogen_files = list_fasta_files(cleaned)
            if not pathogen_files:
                print(f"[warn] No pathogen FASTA in {cleaned}")
                continue

            # Create one sample folder per pathogen file
            for i, pf in enumerate(pathogen_files, start=1):
                sample_dir_name = f"sample{i}_{cat_name}_{pretty}"
                sample_dir = out_top / sample_dir_name
                sample_dir.mkdir(parents=True, exist_ok=True)

                # Copy/symlink all commensals
                for cf in commensals:
                    dst = sample_dir / cf.name
                    if SYMLINK_COMMENSALS:
                        safe_link(cf, dst)
                    else:
                        safe_copy(cf, dst)

                # Copy the pathogen file
                safe_copy(pf, sample_dir / pf.name)

            print(f"[ok] {species_dir_name} / {cat_name}: "
                  f"made {len(pathogen_files)} samples -> {out_top}")

    print("\nAll pathogen samples combined with commensal files successfully.")

if __name__ == "__main__":
    main()
