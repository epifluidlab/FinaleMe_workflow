# FinaleMe Snakemake Workflow

This Snakemake workflow automates the prediction of DNA methylation from cell-free DNA (cfDNA) whole-genome sequencing (WGS) data using FinaleMe. It supports per-sample processing through methylation prediction and an optional downstream tissue-of-origin (TOO) analysis for all samples.

## Key Features

*   **FinaleMe Integration:** Implements the core steps of the FinaleMe analysis pipeline:
    1.  Feature extraction from BAM files.
    2.  HMM model training.
    3.  CpG methylation prediction (decoding).
    4.  Conversion of predictions to BigWig format.
    5.  Optional: Tissue-of-Origin (TOO) analysis using predicted methylation levels.
*   **Input Compatibility:** Designed for coordinate-sorted and indexed BAM files.
*   **Configurability:** Workflow behavior and parameters are controlled via a YAML configuration file.
*   **Parallelization:** Snakemake enables multi-core processing for suitable steps.
*   **SLURM Integration:** Can be adapted for job submission to SLURM clusters (requires a SLURM profile for Snakemake).

## Installation and Setup

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/epifluidlab/FinaleMe_workflow
    cd FinaleMe_workflow
    ```

2.  **Create Conda Environment:**
    Set up a Conda environment with the necessary software dependencies.
    ```bash
    # Install dependencies onto a Conda environment
    conda env create -f environment.yml
    conda activate finaleme_workflow
    ```

3.  **Download FinaleMe JARs:**
    Download the main FinaleMe JAR (`FinaleMe-VERSION-jar-with-dependencies.jar`). Auxiliary JARs can be found in the `lib/` directory of this repository, which should be present upon cloning.
    *   Main JAR: [FinaleMe v0.58.1 Release](https://github.com/epifluidlab/FinaleMe/releases/tag/v.0.58.1)

## Dependencies

This workflow relies on the following tools. It's recommended to install them via Conda (see `environment.yml`).

*   **Java:** (Tested with OpenJDK 1.8.0).
*   **Snakemake:** Workflow management system.
*   **Perl:** Required for the `bedpredict2bw.b37.pl` script (tested with v5.26.3).
*   **UCSC Tools:** Specifically `bedGraphToBigWig` (tested with v4).
    *   Install via Bioconda: `conda install -c bioconda ucsc-bedgraphtobigwig`
*   **Bedtools:** (Tested with v2.29.2).
    *   Install via Bioconda: `conda install -c bioconda bedtools`
*   **Samtools:** For preparing reference genome files (e.g., `faidx`).

## Required Data

Before running the workflow, ensure you have the following data, and their paths are correctly specified in the configuration file:

*   **Input BAM files:** Coordinate-sorted and indexed (`.bai`) BAM files for each sample. An example has been provided in the `input` directory
*   **Reference Genome (2bit format):** E.g., `hg19.2bit`. Can be downloaded from UCSC or converted using `faToTwoBit`.
*   **CpG Motif BedGraph:** E.g., `CG_motif.hg19.common_chr.pos_only.bedgraph`. Available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7779198.svg)](https://doi.org/10.5281/zenodo.7779198).
*   **Exclude Regions BED:** Regions to mask (dark regions), e.g., `wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed`.
*   **Methylation Prior BigWig:** E.g., `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`. Available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7779198.svg)](https://doi.org/10.5281/zenodo.7779198).
*   **Chromosome Sizes File:** E.g., `hg19.chrom.sizes`. Can be generated using `samtools faidx` and `awk`.
*   **FinaleMe Scripts:**
    *   `bedpredict2bw.b37.pl` (for tissue of origin analysis only) - see the `scripts` directory
    <!-- *   `TissueOfOriginExampleScript.R` (for Step 5, if enabled). -->
*   **(Optional - For Tissue of Origin Analysis):**
    *   `autosome_1kb_intervals.UCSC.cpgIsland_plus_shore.b37.bed`: Bed file with 1kb intervals: Download the UCSC.cpgIsland annotation file from UCSC genome browser, keep the autosomes, and generate 1kb non-overlapped windows
    *   Reference methylomes for TOO analysis [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14013719.svg)](https://doi.org/10.5281/zenodo.14013719).

## Configuration (`params.yaml`)

The workflow is controlled by a `params.yaml` file. Check `params.yaml` in this repository for an example of how to configure this workflow. 

## Quick Start

1.  **Prepare Data and Config:**
    *   Organize your input BAMs, supplementary files, and FinaleMe JARs/scripts as per the paths in your `params.yaml`.
    *   Ensure your `params.yaml` is correctly filled out.

2.  **Run Snakemake:**
    Navigate to the directory containing the `Snakefile` and `params.yaml`.
    ```bash
    # Activate the conda environment
    conda activate finaleme_workflow

    # Dry-run to check the workflow plan
    snakemake -n --configfile params.yaml

    # Execute the workflow (adjust --cores and --jobs as needed)
    # --cores: Total number of CPU cores Snakemake can use.
    # --jobs: Maximum number of concurrent jobs (rules) to run.

    # Note: FinaleMe hasn't been updated to utilize multithreading yet, so jobs should ideally match the number of cores
    snakemake --configfile params.yaml --cores <number_of_cores> --jobs <number_of_jobs>
    ```

3.  **SLURM Execution (Optional):**
    This workflow is SLURM compatible if you're using this workflow in an HPC environment as a job.
    ```bash
    # Using the sample SLURM profile provided in this repository
    snakemake --configfile params.yaml --profile slurm_profile > snakemake.log 2>&1 &
    ```
## Workflow Structure and Output

*   **Input BAMs:** Located in the directory specified by `input_dir`.
*   **Supplementary Files:** Reference genomes, annotations, etc., are in `supplement_dir`.
*   **Main Output:** Processed files for each sample are written to subdirectories within `output_dir/{sample_name}/`.
    *   `{sample_name}.CpgMultiMetricsStats.details.bed.gz`: Extracted features (Step 1).
    *   `{sample_name}.finaleme.model`: Trained HMM model (Step 2).
    *   `{sample_name}.finaleme.prediction.bed.gz`: Raw methylation predictions (Step 3).
    *   `{sample_name}.finaleme.cov.b37.bw`: Coverage BigWig (Step 4).
    *   `{sample_name}.finaleme.methy_count.b37.bw`: Methylation count BigWig (Step 4).
*   **Tissue of Origin Output (if `too_enabled: True`):**
    Files related to the TOO analysis are placed in `output_dir/tissue_of_origin/`.
    *   `cfdna.methy_summary.cmd.txt`: Intermediate command file for AlignMultiWig.
    *   `cfdna.names_order.txt`: Intermediate sample order file.
    *   `output.add_value.methy.bed.gz`: Merged methylation data across samples.
    *   `tissue_of_origin_results.tsv` (or similar): Final results from the R script.

## Notes

*   Ensure sufficient memory (`-Xmx` Java options in `params.yaml`) is allocated for each Java step, especially for HMM training and decoding, based on your dataset size and system resources.
*   Log files for each step are generated in `output/{sample_name}/logs/` or `output_dir/tissue_of_origin/logs/`.
