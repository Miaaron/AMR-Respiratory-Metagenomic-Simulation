# Simulated respiratory metagenomic sequencing for AMR prediction

This repository documents the **genome preparation and simulated sequencing workflow** used in the accompanying antimicrobial-resistance (AMR) prediction manuscript. It is intentionally focused on **dataset construction**, not downstream model fitting, KARGA/ResFinder evaluation, or figure generation.

## Pathogens and antibiotics

| Pathogen | Antibiotics |
| --- | --- |
| *Acinetobacter baumannii* | amikacin |
| *Klebsiella pneumoniae* | amikacin, cefepime, ceftazidime, ciprofloxacin, levofloxacin, piperacillin |
| *Pseudomonas aeruginosa* | amikacin, ceftazidime, ciprofloxacin, levofloxacin |

These are 11 pathogen–antibiotic combinations. Resistant/susceptible phenotype labels and pathogen genome sequence data were obtained from **BV-BRC**. The study used a background panel of **87 commensal species** informed by NCBI BioProject **PRJNA945482** and additional respiratory microbiome literature; commensal reference assemblies were retrieved from **NCBI/GenBank**. The `source_genome_ids/` directory provides BV-BRC pathogen genome IDs paired with antibiotic resistance phenotype labels and the available NCBI commensal sequence accessions.

## Data availability

The original genome FASTA assemblies, simulated FASTQ datasets, and intermediate k-mer matrices are **not distributed in this repository because of their size**. Users must obtain the reference genomes from BV-BRC and NCBI, prepare the reference community, and run the simulations using their own local paths and installed software. To reconstruct the original dataset, use the source identifiers in `source_genome_ids/` rather than choosing new assemblies only by species name.

### Source-genome identifiers

- `source_genome_ids/pathogen_genome_ids.csv` - **2,878 pathogen–antibiotic–phenotype records**, covering all **11** analyzed combinations. Columns: `pathogen`, `antibiotic`, `phenotype`, `bvbrc_genome_id`. These are **BV-BRC genome identifiers**, not necessarily NCBI assembly accession numbers. The same genome can occur in multiple antibiotic subsets; the file contains **791 unique genome IDs**. Counts describe the original source FASTA directories, not the final filtered isolates, simulated replicates, or model training observations.
- `source_genome_ids/commensal_sequence_accessions.csv` - **790 NCBI sequence/contig accessions associated with 85 named species**, recovered from saved `.fna.gz` genome headers. These are **sequence** accessions, not necessarily genome **assembly** accessions. The manuscript describes an 87-species commensal panel; **the current accession inventory does not establish the identity of the remaining two species**. Do not treat this file as a complete 87-species assembly manifest without reconciliation.

The FASTA sequences and simulation outputs are not included because of their size. External BV-BRC and NCBI resources may change, and the original tool/reference versions and unrecorded aspects of community randomization can prevent bit-for-bit regeneration.

The `genome_preparation/` directory contains **one study-specific script**:

- `prepare_pathogen_commensal_samples.py` - copies one pathogen genome and the commensal genome FASTAs into per-sample directories used by the simulation workflows.

The original genome-downloading helper scripts are intentionally **not included**. Retrieve the pathogen genomes from **BV-BRC** using `source_genome_ids/pathogen_genome_ids.csv`, and obtain the commensal references from **NCBI/GenBank**, guided by `source_genome_ids/commensal_sequence_accessions.csv` and the study methods. The latter is a sequence-accession inventory, not a complete assembly-level retrieval manifest (85 species represented versus 87 in the manuscript). Downloading by species name alone may return a different assembly. Decompress `.fna.gz` files to FASTA before preparing samples. The preparation script still requires users to set their local directories and follow its original source-file naming assumptions.

## Simulation workflows

- `simulation/illumina/`: InSilicoSeq (ISS), NovaSeq error model, 150-bp paired-end data. The included `Snakefile` and `config.yaml` also implement downstream read processing and k-mer feature generation, reflecting the original workflow.
- `simulation/nanopore_pacbio/`: PBSIM2 long-read simulation, followed by the original workflow's merging and feature preparation. **Nanopore and PacBio used the same basic PBSIM2 pipeline**, but require the correct platform-specific HMM model and `--difference-ratio` option:
  - Nanopore: `23:31:46`
  - PacBio: `6:50:54`

**Important:** `simulation/nanopore_pacbio/Simulation.py` retains the source Nanopore setting `23:31:46` as a literal in its PBSIM2 command. For a PacBio run, change that argument to `6:50:54` and select the appropriate PacBio HMM model. The shared example `config.yaml` currently references `R94.model` and must be adapted accordingly. These are instructions, not an assertion that the included source file will automatically switch platform.

The study examined two pathogen abundances (5% and 10%) and used a nominal microbial sequencing budget of 1,000,000,000 bases per metagenome. **Configuration values differ between organism–antibiotic subsets**; the example `config.yaml` files come from one original subset and must not be interpreted as the values for every run. Replicate counts for resistant and susceptible samples were chosen to address class imbalance.

### Random seeds

**Illumina:** `Step2b_simulate.py` selects an ISS random seed for each reference/simulation and writes it to `simulation_seed_record.json`. The `illumina_seed_records/` directory contains 11 collected subset seed-record files compressed individually as `.json.gz` to make GitHub browser uploads practical. They preserve the original **Illumina** seed assignments; no seed values were changed. To inspect a file without decompressing it, use `gzip -cd illumina_seed_records/PAamikacin_simulation_seed_record.json.gz | head`. To restore the original JSON files, run `gunzip -k illumina_seed_records/*.json.gz` from the repository root (the `.gz` originals remain in place). These files are provenance records, not additional simulation parameters automatically loaded by the pipeline.

**Nanopore and PacBio (PBSIM2):** In `Simulation.py`, `seed = sim_number` and the simulator is invoked with `--seed <seed>`. Thus PBSIM2 seeds are determined by replicate numbers (1, 2, 3, ...). The actual maximum replicate number depends on the subset configuration. The commensal coverage allocation also uses Python's unseeded `random.uniform(0.9, 1.1)`; **the PBSIM2 seed alone does not make community composition fully deterministic**. The supplied log audit checked Nanopore settings, not every PacBio run.

### Running a subset

Install the software invoked by the scripts, including Python dependencies, Snakemake, PBSIM2 (long reads), InSilicoSeq (Illumina), and KMC for the downstream k-mer steps. The `config.yaml` files use `/PATH/TO/...` placeholders; replace these with real input/output directories, HMM-model paths, and executable paths. For example:

```bash
cd simulation/nanopore_pacbio
# Edit config.yaml and verify Simulation.py's PBSIM2 executable,
# --difference-ratio, and HMM model before running.
snakemake -n
snakemake --cores 4
```

Likewise, for Illumina:

```bash
cd simulation/illumina
# Edit config.yaml for the chosen subset and ensure its prerequisites exist.
snakemake -n
snakemake --cores 4
```

These are example commands for the original Snakemake workflows, **not validated end-to-end installation instructions**. The Illumina workflow depends on the genome-specific coverage estimates generated from the corresponding long-read simulation logs. Some original rules perform deletion/cleanup; read the `Snakefile` and verify configured directories before running.

## Scope and limitations

This code is provided to document how sequencing datasets were generated. It has not been retested end-to-end after replacing machine-specific paths with placeholders. Included simulation and sample-preparation scripts were retained apart from earlier machine-specific path substitutions and descriptive renaming of the sample-preparation script. The two genome-download helper scripts are intentionally excluded. The supplied pathogen genome IDs and commensal sequence accessions improve source tracing, but full reproducibility additionally requires reconciliation of the 87-species commensal panel, exact assembly mappings where needed, phenotype definitions, exact tool versions, source data, and the original configuration for each organism–antibiotic setting.
