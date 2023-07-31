# HiCt Sequencing Pipeline

This repository contains a comprehensive pipeline for processing raw sequencing data and performing variant calling on a specific genomic region (RBD or 3CL-Pro) of the SARS-CoV-2 virus. The pipeline is designed to handle multiple samples efficiently, generating consensus sequences and reports for each sample. It utilizes popular bioinformatics tools such as BWA, Samtools, Bcftools, and BBmap to ensure accurate and reliable results.

## Introduction

### Software Dependencies

- BWA (Burrows-Wheeler Aligner): A software package for mapping low-divergent sequences against a large reference genome.
- Samtools: A suite of utilities for interacting with and post-processing SAM/BAM files.
- Bcftools: A utility for variant calling and manipulating VCF (Variant Call Format) files.
- BBmap: A toolset for efficiently aligning DNA reads to a reference genome.

### R Dependencies

- gsheet: A package for reading Google Sheets and creating data frames from them.
- kableExtra: A package for creating complex tables and customizing them.
- staplr: A package for flexible PDF table generation from R data frames.

## Quick Start

To get started with the pipeline, follow these steps:

1. Clone this repository to your local machine.
2. Install the required software dependencies (BWA, Samtools, Bcftools, BBmap) and ensure they are accessible in your system's PATH.
3. Install the R dependencies (gsheet, kableExtra, staplr) in your R environment.
4. Define the necessary parameters in the provided configuration file (config.yaml) to suit your analysis requirements.

## Configuration

The pipeline's behavior can be customized using the configuration file (config.yaml). The configuration file allows you to specify various parameters, including:

- Minimum amplicon length: The minimum length of the amplicons to be considered during read alignment.
- Minimum coverage: The minimum coverage required to consider a position for variant calling.
- Minimum allele frequency: The minimum frequency required to consider a variant.
- Quality thresholds: Various quality thresholds for read filtering, mapping quality, and base quality.
- Google Sheet URL: The URL of the Google Sheet containing metadata for the samples.

## Define Functions

The pipeline utilizes several custom functions to execute specific tasks efficiently. These functions include:

- `run_conda_command`: A function to run a command with a specific Conda environment.
- `make_consensus`: A function to generate consensus sequences for each sample based on variant calling results.
- `major_variants`: A function to filter and save variant calling results based on specified criteria.
- `rename_consensus`: A function to rename the consensus sequence headers to match the sample ID.
- `extract_stat`: A function to extract statistics from the mapping stat file and update the reference.
- `make_report`: A function to generate a comprehensive report for each sample, including general information, coverage plots, and mutation data.

## Execution

The pipeline follows a step-by-step process for each sample:

1. Read metadata from the provided Google Sheet and extract relevant sample information.
2. Record versions of essential tools (BWA, Samtools, Bcftools, BBmap) used for future reference.
3. Execute the necessary commands for each sample, including read alignment, variant calling, consensus sequence generation, and report generation.

## Custom Commands

In addition to functions, the pipeline includes custom commands that offer advanced analysis capabilities:

- `unique_mutations`: A command to generate a reference file of unique nucleotide mutations for coding mutations. This is useful for tracking genetic variations over time.
- `lineage_prevalency`: A command to determine the prevalence of each lineage for specific mutations. This analysis helps in understanding the distribution of mutations across different lineages.

## Usage and Output

After successful execution of the pipeline, you will obtain the following outputs:

- Consensus sequences: FASTA files containing the consensus sequences for each sample.
- Variant calling results: VCF files containing the variant calling results for each sample.
- Coverage plots: PDF files displaying coverage plots for each sample.
- Sample reports: PDF files containing comprehensive reports for each sample, including key metrics, coverage plots, and mutation information.

## Define Parameters

In this section, you can modify various parameters used throughout the code, such as minimum amplicon length, minimum coverage, minimum allele frequency, etc. Adjust these values based on your analysis requirements.

## Define Functions

Several custom functions are defined in this section to carry out specific tasks within the pipeline. These functions include `run_conda_command`, `make_consensus`, `major_variants`, `rename_consensus`, `extract_stat`, `make_report`, `unique_mutations`, and `lineage_prevalency`. Each function serves a crucial role in processing the data and generating reports.

## Initialization

This part initializes the necessary file architecture and reads the metadata from a Google Sheet. Ensure you have provided the correct URL for the Google Sheet containing metadata. The code creates relevant directories and extracts sample names, dates, and VSPs from the metadata.

## Recording Version History

The pipeline records the versions of essential tools used, including BWA, Samtools, Bcftools, and BBmap. The version information is saved in the reference file for future reference.

## Execute commands for each sample

In this section, the code executes the necessary commands for each sample. It aligns the reads to the reference, converts SAM files to BAM files, filters reads based on size, quality, and mapping quality, performs variant calling, generates consensus sequences, and generates reports. This part iterates through each sample to process the data efficiently.

## Custom Functions

### run_conda_command

This function is used to run a command with a specific Conda environment. It writes a bash script, activates the Conda environment, executes the provided command, and deactivates the environment.

### make_consensus

This function generates consensus sequences for each sample based on the variant calling results. It reads the SAMtools pileup, calculates base occurrences, and assigns the consensus nucleotide.

### major_variants

This function filters and saves variant calling results based on specified criteria such as minimum coverage, allele frequency, mapping quality, and base quality.

### rename_consensus

This function renames the consensus sequence headers to match the sample ID.

### extract_stat

This function extracts statistics from the mapping stat file and updates the reference with information on read alignment and coverage.

### make_report

This function generates a report for each sample. It creates three pages: the first page contains general sample information, the second page displays coverage plots, and the third page includes mutation information.

## Custom Commands

### unique_mutations

This command generates a reference file of unique nucleotide mutations for coding mutations. It combines data from the "mutations.csv" file and determines unique mutations based on gene and type.

### lineage_prevalency

This command determines the prevalence of a lineage for each specific mutation. It analyzes the "mutations.csv" file and generates a matrix with mutations as columns and lineages as rows, representing the prevalence of each mutation in each lineage.

Please note that this README is only a brief overview of the code's functionalities. For a more in-depth understanding, refer to the code comments and the functions' implementations.

## Additional Notes

- The pipeline has been optimized for the RBD region of the SARS-CoV-2 virus, but it can be adapted for other genomic regions with minor modifications.
- The pipeline is designed to handle multiple samples simultaneously, ensuring efficient processing and analysis of large datasets.

For more detailed instructions and explanations, please refer to the code comments, function implementations, and example files provided in the repository.

Happy variant calling!

