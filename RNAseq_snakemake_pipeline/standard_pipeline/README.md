# RNA-seq snakemake pipeline - standard pipeline

This is a snakemake pipeline for RNA-seq analysis of repetitive element expression. This is a basic pipeline that process RNA-seq data according to our standard workflow and generates combined count tables containing raw counts from TEtranscripts and Telescope. TEtranscripts is a software that quantifies repetitive element expression at the subfamily level, while Telescope quantifies repetitive element expression at the locus-specific level. This pipeline generates two count tables, one for TEtranscripts and one for Telescope, containing gene + subfamily TE raw counts and gene + locus TE raw counts, respectively for all samples. This pipeline is designed to be run on a SLURM cluster.

The standard pipeline steps are as follows:

1. Initial QC - FastQC
2. Adapter/read trimming & 2nd pass QC - TrimGalore
3. Read alignment - STAR
4. QC report generation - MultiQC
5. Gene and TE subfamily read calling - TEtranscripts
6. Locus-specific read calling - Telescope
7. Combine count files from TEtranscripts and Telescope - combine_counts

This pipeline assumes paired end data!

This pipeline is designed to be run on a SLURM cluster. The cluster configuration file is `RNAseq.standard.Snakemake.cluster.config.yaml`. The pipeline configuration file is `RNAseq.standard.Snakemake.config.yaml`.

## Setup:

### Create a snakemake pipeline environment

#### 1. Download or copy the snakemake.yaml file

#### 2. Load miniconda

`ml miniconda`

#### 3. Create the new environment

`conda env create -n snakemake -f snakemake.yaml`

#### 4. Activate the envirionment

`conda activate snakemake`

### Raw data requirements:

- Raw files must be in the format: `<sample name>_<read>.<fastq file extension>`
  - `<sample name>` Can be anything you want, but do not have any hyphens in the name! Underscores or dots are fine.
  - `<read>` R1 for read 1 and R2 for read 2.
  - `<fastq file extension>` This is set in the config file. Must be gzipped fastq files (e.g., fastq.gz, fq.gz).

### Folder structure:

- Raw data should be located in folder called `raw_files`.
- `RNAseq.standard.Snakefile`, `RNAseq.standard.Snakemake.cluster.config.yaml`, and `RNAseq.standard.Snakemake.config.yaml` must be in the working directory.
- All environment yaml files should be located in a folder called `envs`.
- Environment yaml files must be in the envs/ folder. The environment yaml files are:
  - fastqc.yaml
  - trimgalore.yaml
  - tetranscripts.yaml
  - telescope.yaml
  - combine_counts.yaml
- All other subfolders will be created automatically.

## To run the pipeline:

1. Start a new tmux session with the command `tmux new-session -s <session name>`
2. Activate the snakemake enviroinment with the command `conda activate snakemake`
3. Navigate to the working directory where the Snakefile, cluster config, and config files are located as well as the raw data and envs subfolders.
4. Run the snakemake command:

```bash
snakemake -s RNAseq.standard.Snakefile -j 100 --use-conda --configfile RNAseq.standard.Snakemake.config.yaml --cluster-config RNAseq.standard.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
```

Snakemake will create two tsv files located in `results/` containing the raw counts from TEtranscripts and Telescope. These are called `tetranscripts_counts.tsv` and `telescope_counts.tsv`. These count tables are ready for downstream analysis, for example, differential expression analysis.
