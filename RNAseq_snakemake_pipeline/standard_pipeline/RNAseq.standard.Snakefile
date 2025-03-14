#type: ignore
"""
This is a snakemake pipeline for RNA-seq analysis of repetitive element expression. This is the standard pipeline.

The pipeline steps are as follows:

1. Initial QC - FastQC
2. Adapter/read trimming and 2nd pass QC on trimmed data - TrimGalore
3. Read alignment - STAR
4. Alignment QC - CollectRnaSeqMetrics from Picard
5. QC report generation - MultiQC
6. Gene and TE subfamily read calling - TEtranscripts
7. Locus-specific TE read calling - Telescope
8. Locus-specific TE read calling - TElocal
9. Combine count files from TEtranscripts, Telescope, and TElocal - combine_counts

REQUIREMENTS:
Raw data:
- This pipeline assumes paired end reads
- Raw files must be in the format <sample name>_<read>.<fastq file extension>
    <sample name> - can be anything you want
    <read> - R1 for read 1 and R2 for read 2
    <fastq file extension> - set in the config file. Must be gzipped fastq files (e.g., fastq.gz)

Folder structure:
- Raw data should be located in folder called raw_files
- This file, RNAseq.standard.Snakemake.cluster.config.yaml, and RNAseq.standard.Snakemake.config.yaml must be in the working directory
- Environment yaml files must be in the envs/ folder. The environment yaml files are:
    - fastqc.yaml
    - trimGalore.yaml
    - picard.yaml
    - tetranscripts.yaml
    - telescope.yaml
    - telocal.yaml
    - combine_counts.yaml
- All subfolders will be created automatically

TO RUN THE PIPELINE:
1. start a new tmux session with the command tmux new-session -s <session name>
2. Run the following command:
snakemake -s RNAseq.standard.Snakefile -j 100 --configfile RNAseq.standard.Snakemake.config.yaml --cluster-config RNAseq.standard.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
3. Snakemake will create four count tables located in the results/ folder:
    - telescope_counts.tsv # this file is ready for input into DESeq2
    - tetranscripts_counts.tsv # this file is ready for input into DESeq2
    - telocal_counts.tsv # this file is ready for input into DESeq2
    - telocal_counts_annotated.tsv # this file contains the RepeatMasker annotations and chromosome positions
"""

import os
from os import path

# Load config files
configfile: "RNAseq.standard.Snakemake.config.yaml"

# Specify the local rules
localrules: MultiQC, create_count_file_list, combine_counts, all

#### Validate the config file ###############################################
# Check that the config file exists
assert 'working_dir' in config, "Missing field in configfile -- 'working_dir'"
assert 'raw_file_extension' in config, "Missing field in configfile -- 'raw_file_extension'"
assert 'read_length' in config, "Missing field in configfile -- 'read_length'"
assert 'TEtrx_strandedness' in config, "Missing field in configfile -- 'TEtrx_strandedness'"
assert 'Tel_strandedness' in config, "Missing field in configfile -- 'Tel_strandedness'"
assert 'Picard_strandedness' in config, "Missing field in configfile -- 'Picard_strandedness'"
assert 'STAR_genomeDir' in config, "Missing field in configfile -- 'STAR_genomeDir'"
assert 'STAR_GTF' in config, "Missing field in configfile -- 'STAR_GTF'"
assert 'REF_FLAT' in config, "Missing field in configfile -- 'REF_FLAT'"
assert 'rRNA_interval_list' in config, "Missing field in configfile -- 'rRNA_interval_list'"
assert 'TEtrx_TE' in config, "Missing field in configfile -- 'TEtrx_TE'"
assert 'Telescope_GTF' in config, "Missing field in configfile -- 'Telescope_GTF'"
assert 'TEloc_TE' in config, "Missing field in configfile -- 'TEloc_TE'"

# Check that each field is filled in and properly formatted
# Working dir
assert len(config['working_dir']) > 0, "config file: Please provide a working directory (working_dir)."
assert path.exists(config['working_dir']), "config file: working_dir " + config['working_dir'] + " does not exist."
assert config['working_dir'].endswith('/'), "config file: working_dir must end in '/'!"
# Raw file extension
assert len(config['raw_file_extension']) > 0, "config file: Please provide a raw file extension (raw_file_extension)."
assert config['raw_file_extension'].startswith('.') == False, "config file: raw_file_extension should not start with '.'"
if config['raw_file_extension'].endswith('gz') == False: #just warn if it doesnt end in .gz -- files might still be gzipped
    print("WARNING: config file: raw_file_extension does not end in 'gz'. Raw files must be gzipped!")
# Read length
assert config['read_length'] > 0, "config file: Please provide a read length (read_length)."
# TEtranscripts strandedness
assert config['TEtrx_strandedness'] in ['no', 'forward', 'reverse'], "config file: TEtranscripts strandedness must be 'no', 'forward', or 'reverse'."
# Telescope strandedness
assert config['Tel_strandedness'] in ['None', 'FR', 'RF'], "config file: Telescope strandedness must be 'None', 'FR', or 'RF'."
# Picard strandedness
assert config['Picard_strandedness'] in ['NONE', 'FIRST_READ_TRANSCRIPTION_STRAND', 'SECOND_READ_TRANSCRIPTION_STRAND'], "config file: Picard strandedness must be 'NONE', 'FIRST_READ_TRANSCRIPTION_STRAND', or 'SECOND_READ_TRANSCRIPTION_STRAND'."
# STAR genomeDir
assert len(config['STAR_genomeDir']) > 0, "config file: Please provide a STAR genomeDir (STAR_genomeDir)."
assert path.exists(config['STAR_genomeDir']), "config file: STAR_genomeDir " + config['STAR_genomeDir'] + " does not exist."
# STAR GTF
assert len(config['STAR_GTF']) > 0, "config file: Please provide a STAR GTF file (STAR_GTF)."
assert path.exists(config['STAR_GTF']), "config file: STAR_GTF " + config['STAR_GTF'] + " does not exist."
# REF_FLAT
assert len(config['REF_FLAT']) > 0, "config file: Please provide a REF_FLAT file (REF_FLAT)."
assert path.exists(config['REF_FLAT']), "config file: REF_FLAT " + config['REF_FLAT'] + " does not exist."
# rRNA_intervals
assert len(config['rRNA_interval_list']) > 0, "config file: Please provide a rRNA interval list file (rRNA_interval_list)."
assert path.exists(config['rRNA_interval_list']), "config file: rRNA_interval_list " + config['rRNA_interval_list'] + " does not exist."
# TEtranscripts TE
assert len(config['TEtrx_TE']) > 0, "config file: Please provide a TEtranscripts TE file (TEtrx_TE)."
assert path.exists(config['TEtrx_TE']), "config file: TEtrx_TE " + config['TEtrx_TE'] + " does not exist."
# Telescope GTF
assert len(config['Telescope_GTF']) > 0, "config file: Please provide a Telescope GTF file (Telescope_GTF)."
assert path.exists(config['Telescope_GTF']), "config file: Telescope_GTF " + config['Telescope_GTF'] + " does not exist."
# TElocal TE
assert len(config['TEloc_TE']) > 0, "config file: Please provide a TElocal TE file (TEloc_TE)."
assert path.exists(config['TEloc_TE']), "config file: TEloc_TE " + config['TEloc_TE'] + " does not exist."

# Store some of the config values as variables
working_dir = config['working_dir']
raw_file_ext = config['raw_file_extension']

# Create output directories for each rule
for rule in ["FastQC", "TrimGalore", "STAR", "CollectRnaSeqMetrics", "TEtranscripts", "Telescope", "TElocal"]:
    os.makedirs(f"outputs/{rule}", exist_ok=True)

############################################################################################################

print("Successfully validated the configuration files! Starting the pipeline...")

#### Get sample IDs ##########################################################
# Store the sample IDs in a list
sample_ids = glob_wildcards(working_dir + "raw_files/{sample}_{read}" + raw_file_ext).sample
# Remove the duplicate sample IDs
sample_ids = list(set(sample_ids))
print("Sample IDs: " + str(sample_ids))

# Validate the sample IDs
assert len(sample_ids) > 0, "No sample IDs found in the raw files directory."

#### Rules ###################################################################
rule all:
    input:
        expand(working_dir + "FastQC/{sample}_{read}_fastqc.html", sample=sample_ids, read=["R1", "R2"]),
        expand(working_dir + "FastQC/{sample}_{read}_fastqc.zip", sample=sample_ids, read=["R1", "R2"]),
        expand(working_dir + "trimgalore/{sample}_R1_val_1.fq.gz", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_R2_val_2.fq.gz", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_R1_val_1_fastqc.zip", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_R2_val_2_fastqc.zip", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_R1_val_1_fastqc.html", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_R2_val_2_fastqc.html", sample=sample_ids),
        expand(working_dir + "trimgalore/{sample}_{read}.fastq.gz_trimming_report.txt", sample=sample_ids, read=["R1", "R2"]),
        expand(working_dir + "star/{sample}/{sample}_Aligned.out.bam", sample=sample_ids),
        expand(working_dir + "star/{sample}/sorted_{sample}_Aligned.out.bam", sample=sample_ids),
        expand(working_dir + "star/{sample}/{sample}_Log.final.out", sample=sample_ids),
        expand(working_dir + "star/{sample}/{sample}_Log.out", sample=sample_ids),
        expand(working_dir + "star/{sample}/{sample}_Log.progress.out", sample=sample_ids),
        expand(working_dir + "star/{sample}/{sample}_SJ.out.tab", sample=sample_ids),
        expand(working_dir + "CollectRnaSeqMetrics/{sample}_rnaseq_metrics.txt", sample=sample_ids),
        expand(working_dir + "TEtranscripts/{sample}-tetranscripts.cntTable", sample=sample_ids),
        expand(working_dir + "telescope/{sample}-TE_counts.tsv", sample=sample_ids),
        expand(working_dir + "TElocal/{sample}-telocal.cntTable", sample=sample_ids),
        working_dir + "multiqc/multiqc_report.html",
        working_dir + "lists/combined_count_files.txt",
        working_dir + "results/telescope_counts.tsv",
        working_dir + "results/tetranscripts_counts.tsv",
        working_dir + "results/telocal_counts.tsv",
        working_dir + "results/telocal_counts_annotated.tsv"

    shell:
        "echo '--- Snakemake Pipeline Completed Successfully! ---'"

rule FastQC:
    message: "Running FastQC for {wildcards.sample}_{wildcards.read}"

    input:
        working_dir + 'raw_files/{sample}_{read}.' + raw_file_ext
    
    output:
        working_dir + 'FastQC/{sample}_{read}_fastqc.html',
        working_dir + 'FastQC/{sample}_{read}_fastqc.zip'
    
    log: working_dir + 'logs/FastQC/FastQC.{sample}_{read}.log'
    
    params:
        output_dir = working_dir + 'FastQC'
    
    conda:
        config['envs']['fastqc']
    
    shell:
        """
        # Create the output directory
        mkdir -p {params.output_dir}
        
        # Run FastQC
        fastqc \
            --outdir {params.output_dir} \
            {input} \
            &> {log}
        """

rule TrimGalore:
    message: "Running TrimGalore for {wildcards.sample}"
    
    input:
        r1 = working_dir + 'raw_files/{sample}_R1.' + raw_file_ext,
        r2 = working_dir + 'raw_files/{sample}_R2.' + raw_file_ext
    
    output:
        working_dir + 'trimgalore/{sample}_R1_val_1.fq.gz',
        working_dir + 'trimgalore/{sample}_R2_val_2.fq.gz',
        working_dir + 'trimgalore/{sample}_R1_val_1_fastqc.zip',
        working_dir + 'trimgalore/{sample}_R2_val_2_fastqc.zip',
        working_dir + 'trimgalore/{sample}_R1_val_1_fastqc.html',
        working_dir + 'trimgalore/{sample}_R2_val_2_fastqc.html',
        working_dir + 'trimgalore/{sample}_R1.fastq.gz_trimming_report.txt',
        working_dir + 'trimgalore/{sample}_R2.fastq.gz_trimming_report.txt'
    
    log:
        working_dir + 'logs/trimgalore/{sample}.log'
    
    params:
        output_dir = working_dir + 'trimgalore',
    
    conda:
        config['envs']['trimgalore']
    
    shell:
        """
        # Create the output directory
        mkdir -p {params.output_dir}

        # Run TrimGalore
        trim_galore \
            --paired \
            --gzip \
            --fastqc \
            --output_dir {params.output_dir} \
            {input.r1} {input.r2} \
            &> {log}
        """

rule STAR:
    message: "Running STAR for {wildcards.sample}"
    
    input: 
        read1 = working_dir + 'trimgalore/{sample}_R1_val_1.fq.gz',
        read2 = working_dir + 'trimgalore/{sample}_R2_val_2.fq.gz'
    
    output:
        bam = working_dir + 'star/{sample}/{sample}_Aligned.out.bam',
        sorted_bam = working_dir + 'star/{sample}/sorted_{sample}_Aligned.out.bam',
        final_log = working_dir + 'star/{sample}/{sample}_Log.final.out',
        log_out = working_dir + 'star/{sample}/{sample}_Log.out',
        progress = working_dir + 'star/{sample}/{sample}_Log.progress.out',
        sj = working_dir + 'star/{sample}/{sample}_SJ.out.tab',

    log: 
        working_dir + 'logs/star/{sample}.log'
    
    params:
        outdir = working_dir + 'star/{sample}/',
        prefix = working_dir + 'star/{sample}/{sample}_',
        genome_dir = config['STAR_genomeDir'],
        gtf = config['STAR_GTF'],
        read_length = config['read_length'],
        threads = 16,
        multimap_nmax = 100,
        winanchor_multimap = 200
        
    shell:
        """
        # Create the output directory
        mkdir -p {params.outdir}
        
        # Load the STAR module
        ml star/2.7.0e

        # Run STAR
        STAR \
            --runThreadN {params.threads} \
            --genomeDir {params.genome_dir} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.read_length} \
            --readFilesIn {input.read1} {input.read2} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --winAnchorMultimapNmax {params.winanchor_multimap} \
            --outFilterMultimapNmax {params.multimap_nmax} \
            --outFileNamePrefix {params.prefix} \
            &> {log}
        
        # Sort the BAM file using samtools. This is necessary for Picard's CollectRnaSeqMetrics.
        ml samtools
        samtools sort -o {output.sorted_bam} {output.bam}
        """

rule CollectRnaSeqMetrics:
    message: "Running CollectRnaSeqMetrics for {wildcards.sample}"
    
    input:
        bam = working_dir + 'star/{sample}/sorted_{sample}_Aligned.out.bam',
    
    output:
        working_dir + 'CollectRnaSeqMetrics/{sample}_rnaseq_metrics.txt',
    
    log:
        working_dir + 'logs/CollectRnaSeqMetrics/{sample}_rnaseq_metrics.log'
    
    params:
        outdir = working_dir + 'CollectRnaSeqMetrics',
        ref_flat = config['REF_FLAT'],
        rRNA_intervals = config['rRNA_interval_list'],
        strandedness = config['Picard_strandedness']
    
    conda:
        config['envs']['picard']
    
    shell:
        """
        # Create the output directory
        mkdir -p {params.outdir}

        # Run CollectRnaseqMetrics
        picard CollectRnaSeqMetrics \
            --INPUT {input.bam} \
            --OUTPUT {output} \
            --REF_FLAT {params.ref_flat} \
            --RIBOSOMAL_INTERVALS {params.rRNA_intervals} \
            --STRAND_SPECIFICITY {params.strandedness} \
            &> {log}
        """

rule MultiQC:
    message: "Running MultiQC"
    
    input:
        # Gather the FastQC reports from raw reads
        fastqc_html = expand(working_dir + "FastQC/{sample}_{read}_fastqc.html", sample=sample_ids, read=["R1", "R2"]),
        fastqc_zip = expand(working_dir + "FastQC/{sample}_{read}_fastqc.zip", sample=sample_ids, read=["R1", "R2"]),
        trimmed_read1_html = expand(working_dir + 'trimgalore/{sample}_R1_val_1_fastqc.html', sample=sample_ids),
        trimmed_read1_zip = expand(working_dir + 'trimgalore/{sample}_R1_val_1_fastqc.zip', sample=sample_ids),
        trimmed_read2_html = expand(working_dir + 'trimgalore/{sample}_R2_val_2_fastqc.html', sample=sample_ids),
        trimmed_read2_zip = expand(working_dir + 'trimgalore/{sample}_R2_val_2_fastqc.zip', sample=sample_ids),
        star = expand(working_dir + "star/{sample}/{sample}_Log.final.out", sample=sample_ids),
        picard = expand(working_dir + "CollectRnaSeqMetrics/{sample}_rnaseq_metrics.txt", sample=sample_ids)
    
    output:
        working_dir + "multiqc/multiqc_report.html",
    
    params:
        outdir = working_dir + "multiqc"
    
    log:
        working_dir + "logs/multiqc/multiqc.log"
    
    shell:
        """
        # Load module
        ml multiqc

        # Create output directory
        mkdir -p {params.outdir}
        
        multiqc \
            --force \
            --outdir {params.outdir} \
            {input} \
            &> {log}
        """

rule alignment_summary:
    message: "Creating alignment summary"

    input:
        star_multiqc = working_dir + "multiqc/multiqc_data/multiqc_star.txt",
        picard_multiqc = working_dir + "multiqc/multiqc_data/multiqc_picard_RnaSeqMetrics.txt"

    output:
        summary = working_dir + "results/alignment_summary.csv"

    run:
        import pandas as pd
        
        # Read STAR alignment stats
        star_df = pd.read_csv(input.star_multiqc, sep="\t")
        star_cols = ["Sample", "total_reads", "uniquely_mapped", "uniquely_mapped_percent", "multimapped", "multimapped_percent"]
        star_df = star_df[star_cols]

        # Our data is paired-end, so we need to multiply the read counts by 2. This is because STAR counts a paired-end read as one read.
        star_df["total_reads"] = star_df["total_reads"] * 2
        star_df["uniquely_mapped"] = star_df["uniquely_mapped"] * 2
        star_df["multimapped"] = star_df["multimapped"] * 2

        # Read Picard RNA metrics
        picard_df = pd.read_csv(input.picard_multiqc, sep="\t")
        picard_cols = ["Sample", "PCT_RIBOSOMAL_BASES", "PCT_MRNA_BASES"]
        picard_df = picard_df[picard_cols]

        # Remove "sorted_" prefix from sample names
        picard_df["Sample"] = picard_df["Sample"].str.replace("^sorted_", "", regex=True)

        # Convert rRNA percentage from fraction to percent
        picard_df["% rRNA"] = (picard_df["PCT_RIBOSOMAL_BASES"] * 100).round(2)

        # Convert mRNA percentage from fraction to percent
        picard_df["% mRNA"] = (picard_df["PCT_MRNA_BASES"] * 100).round(2)

        # Drop the original columns
        picard_df.drop(columns=["PCT_RIBOSOMAL_BASES", "PCT_MRNA_BASES"], inplace=True)

        # Merge both tables
        summary_df = star_df.merge(picard_df, on="Sample", how="outer")

        # Save the output
        summary_df.to_csv(output.summary, index=False)

rule TEtranscripts:
    message: "Running TEtranscripts for {wildcards.sample}"
    
    input:
        working_dir + 'star/{sample}/{sample}_Aligned.out.bam',
    
    output:
        working_dir + 'TEtranscripts/{sample}-tetranscripts.cntTable',

    log: 'logs/TEtranscripts/TEtranscripts.{sample}.log'

    params:
        output_dir = working_dir + 'TEtranscripts',
        strandedness = config['TEtrx_strandedness'],
        gtf = config['STAR_GTF'],
        te = config['TEtrx_TE']
    
    conda:
        config['envs']['tetranscripts']

    shell:
        """
        # Create the output directory
        mkdir -p {params.output_dir}

        # Run TEtranscripts
        TEcount --format BAM --mode multi --stranded {params.strandedness} \
        -b {input} \
        --GTF {params.gtf} \
        --TE {params.te} \
        --project TEtranscripts/{wildcards.sample}-tetranscripts \
        &> {log}
        """

rule Telescope:
    message: "Running Telescope for {wildcards.sample}"

    input:
        working_dir + 'star/{sample}/{sample}_Aligned.out.bam',
    
    output:
        working_dir + 'telescope/{sample}-TE_counts.tsv',
        working_dir + 'telescope/{sample}-checkpoint.npz',
        working_dir + 'telescope/{sample}-run_stats.tsv'
    
    log: 'logs/telescope/telescope.{sample}.log'

    params:
        output_dir = working_dir + 'telescope',
        strandedness = config['Tel_strandedness'],
        gtf = config['Telescope_GTF']
    
    conda:
        config['envs']['telescope']
    
    shell:
        """
        # Create the output directory
        mkdir -p {params.output_dir}

        # Run Telescope
        telescope assign {input} {params.gtf} \
        --stranded_mode {params.strandedness} \
        --outdir {params.output_dir} \
        --exp_tag {params.output_dir}/{wildcards.sample} \
        &> {log}
        """

rule TElocal:
    message: "Running TElocal for {wildcards.sample}"
    
    input:
        working_dir + 'star/{sample}/{sample}_Aligned.out.bam',
    
    output:
        working_dir + 'TElocal/{sample}-telocal.cntTable',

    log: 'logs/TElocal/TElocal.{sample}.log'

    params:
        output_dir = working_dir + 'TElocal',
        strandedness = config['TEtrx_strandedness'],
        gtf = config['STAR_GTF'],
        te = config['TEloc_TE']
    
    conda:
        config['envs']['telocal']

    shell:
        """
        # Create the output directory
        mkdir -p {params.output_dir}

        # Run TElocal
        TElocal --mode multi --stranded {params.strandedness} \
        -b {input} \
        --GTF {params.gtf} \
        --TE {params.te} \
        --project TElocal/{wildcards.sample}-telocal \
        &> {log}
        """

rule create_count_file_list:
    message: "Creating a list of count files for combining"
    
    input:
        telescope_files = expand(working_dir + "telescope/{sample}-TE_counts.tsv", sample=sample_ids),
        tetranscripts_files = expand(working_dir + "TEtranscripts/{sample}-tetranscripts.cntTable", sample=sample_ids),
        telocal_files = expand(working_dir + "TElocal/{sample}-telocal.cntTable", sample=sample_ids)
    
    output:
        combined_list = working_dir + "lists/combined_count_files.txt"
    
    params:
        output_dir = working_dir + "lists"
    
    run:
        # Create the output directory
        os.makedirs(params.output_dir, exist_ok=True)

        # Write both telescope and tetranscripts file paths to a single file
        with open(output.combined_list, 'w') as f:
            # Write telescope file paths
            for file in input.telescope_files:
                f.write(f"{path.abspath(file)}\n")
            
            # Write tetranscripts file paths
            for file in input.tetranscripts_files:
                f.write(f"{path.abspath(file)}\n")
            
            # Write telocal file paths
            for file in input.telocal_files:
                f.write(f"{path.abspath(file)}\n")

rule combine_counts:
    message: "Combining count files"
    
    input:
        combined_list = working_dir + "lists/combined_count_files.txt"
    
    output:
        telescope_counts = working_dir + "results/telescope_counts.tsv",
        tetranscripts_counts = working_dir + "results/tetranscripts_counts.tsv",
        telocal_counts = working_dir + "results/telocal_counts.tsv",
        telocal_counts_annt = working_dir + "results/telocal_counts_annotated.tsv"
    
    params:
        output_dir = working_dir + "results"
    
    conda:
        config['envs']['combine_counts']

    shell:
        """
        # Make the output directory
        mkdir -p {params.output_dir}
        
        # Run for telescope counts
        python scripts/combine_counts.py {input.combined_list} {params.output_dir} -mode telescope

        # Run for tetranscripts counts
        python scripts/combine_counts.py {input.combined_list} {params.output_dir} -mode tetranscripts

        # Run for telocal counts
        python scripts/combine_counts.py {input.combined_list} {params.output_dir} -mode telocal
        """
