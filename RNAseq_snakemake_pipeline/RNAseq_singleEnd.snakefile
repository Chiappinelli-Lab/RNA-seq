"""
This is a snake pipeline for single-end RNA-seq analysis of repetitive element expression. The Pipeline steps are as follows:

1. Initial QC - FastQC
2. Adapter/read trimming -  cutadapt
3. 2nd pass QC on trimmed data - FastQC
4. Read alignment - STAR
5. Read calling - telescope
6. Differential expression - DESeq2

REQUIREMENTS:
Raw data:
-Raw files must be in the format <sample name>_<replicate>.<fastq file extension>
	<sample name> - can be anything you want
	<replicate> - a single number (all replicaetes for a given sample must be sequential, i.e. 1,2,3, NOT 1,3,4,)
	<fastq file extension> - set in the config file. Must be gzipped fastq files (e.g., fastq.gz, fq.gz)!
-All samples must have the same number of replicates

Folder structure:
-Raw data should be located in folder called raw_files
-This file, RNAseq.Snakemake.cluster.config.yaml, and RNAseq.Snakemake.config.yaml must be in the working directory
-You need to create a new directory called outputs (sbatch outputs will go here)
-All other subfolders will be created automatically

TO RUN THE PIPELINE:
1. start a new tmux session with the command tmux new-session -s <session name>
2. Run the following command:
snakemake -s RNAseq.Snakefile -j 100 --configfile RNAseq.Snakemake.config.yaml --cluster-config RNAseq.Snakemake.cluster.config.yaml --cluster "sbatch -o {cluster.output} -e {cluster.err} -p {cluster.p} -N {cluster.N} -J {cluster.jobName} -t {cluster.time} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type}"
3. Snakemake will create two sample_table.txt files located in telescope/ and TEtranscripts/. Add columns to these files with
	any extra information about each sample. These will appear as extra columns in the final all.samples.DESeq2.tibble.tsv
	files.
"""

import numpy as np
import os.path
from os import path

localrules: Telescope_DESeq, make_sample_tables, Combine_tables_telescope, Combine_tables_TEtranscripts

#### validate the config file ######################################################################
# check that each field exists
assert 'working_dir' in config, "Missing field in configfile -- 'working_dir'"
assert 'raw_file_extension' in config, "Missing field in configfile -- 'raw_file_extension'"
assert 'num_replicates' in config, "Missing field in configfile -- 'num_replicates'"
assert 'control_sample' in config, "Missing field in configfile -- 'control_sample'"
assert 'fw_adapter' in config, "Missing field in configfile -- 'fw_adapter'"
assert 'rev_adapter' in config, "Missing field in configfile -- 'rev_adapter'"

# Check that each field is filled in and properly formatted
#working_dir
assert len(config['working_dir']) > 0, "config file: Please provide a working directory (working_dir)."
assert path.exists(config['working_dir']), "config file: working_dir " + config['working_dir'] + " does not exist."
assert config['working_dir'].endswith('/'), "config file: working_dir must end in '/'!"
#raw_file_extension
assert len(config['raw_file_extension']) > 0, "config file: Please provide a raw file extension (raw_file_extension)."
assert config['raw_file_extension'].startswith('.') == False, "config file: raw_file_extension should not start with '.'"
if config['raw_file_extension'].endswith('gz') == False: #just warn if it doesnt end in .gz -- files might still be gzipped
	print("WARNING: config file: raw_file_extension does not end in 'gz'. Raw files must be gzipped!")
#num_replicates
assert isinstance(config['num_replicates'], int), "config file: num_replicates must be an integer number."
assert config['num_replicates'] > 0, "config file: num_repliates must be at least 1 (single samples have 1 replicate)."
#control_sample
assert len(config['control_sample']) > 0, "config file: Please designate a control sample (control_sample)."
#fw_adapter
#to-do: check that adapter sequence contains only ATGC characters (overkill?)
assert len(config['fw_adapter']) > 0, "config file: Please provide a forward adapter sequence (fw_adapter)."
#STAR_genomeDir
assert len(config['STAR_genomeDir']) > 0, "config file: Please provide a STAR genome Directory (STAR_genomeDir)."
assert path.exists(config['STAR_genomeDir']), "config file: STAR_genomeDir " + config['STAR_genomeDir'] + " does not exist."
#working_dir
assert len(config['STAR_GTF']) > 0, "config file: Please provide a STAR annotation file path (STAR_GTF)."
assert path.exists(config['STAR_GTF']), "config file: STAR_GTF " + config['STAR_GTF'] + " does not exist."
#working_dir
assert len(config['TEtrx_GTF']) > 0, "config file: Please provide a TEtranscripts gene annotation file (TEtrx_GTF)."
assert path.exists(config['TEtrx_GTF']), "config file: TEtrx_GTF " + config['TEtrx_GTF'] + " does not exist."
#working_dir
assert len(config['TEtrx_TE']) > 0, "config file: Please provide a TEtranscripts RE annotation file (TEtrx_TE)."
assert path.exists(config['TEtrx_TE']), "config file: TEtrx_TE " + config['TEtrx_TE'] + " does not exist."
#working_dir
assert len(config['telescope_GTF']) > 0, "config file: Please provide a telescope annotation file (telescope_GTF)."
assert path.exists(config['telescope_GTF']), "config file: telescope_GTF " + config['telescope_GTF'] + " does not exist."

#put some of the config values into variables for convenience
working_dir = config['working_dir']
raw_file_ext = config['raw_file_extension']
replicates = list(range(1, int(config['num_replicates']) + 1))
controlSample = config['control_sample']


# Check if outputs directory exists. If not, create it.
outputs_path = working_dir + "outputs"
if not path.exists(outputs_path):
	print(outputs_path + " directory does not exist!")
	os.mkdir(outputs_path) # hoping that python will handle the error itslef if the directory can not be created
	print ("Successfully created the directory %s " % outputs_path)

####################################################################################################

print("Starting Snakemake...")

#### Get sample IDs ################################################################################
SAMPLE_IDS = glob_wildcards(working_dir + 'raw_files/{sample}_{replicate}.' + raw_file_ext).sample
# this will grab multiple of the same SAMPLE_ID (one for each replicate). We want just the unique ones
# convert to a set (only takes unique values)
SAMPLE_ID_set = set(SAMPLE_IDS)
# convert back to a list
SAMPLE_IDS = list(SAMPLE_ID_set)
# create list minus control
SAMPLE_IDS_no_control = list()
for sample in SAMPLE_IDS:
    if controlSample in sample:
         continue
    SAMPLE_IDS_no_control.append(sample)

print("SAMPLE_IDS = " + str(SAMPLE_IDS))
print("SAMPLE_IDS_no_control = " + str(SAMPLE_IDS_no_control))

#validate SAMPLE_IDS exist
assert len(SAMPLE_IDS) > 1, "No samples found!"

#### Rules ########################################################################################
rule all:
	input:
		expand(working_dir + 'FastQC/{sample}_{replicate}_fastqc.html', sample = SAMPLE_IDS, replicate = replicates),
		expand(working_dir + 'FastQC_2/{sample}_{replicate}.cutadapt.q20.minlen1_fastqc.html', sample = SAMPLE_IDS, replicate = replicates),
		# expand(working_dir + 'cutadapt/{sample}_{replicate}.cutadapt.q20.minlen1.' + raw_file_ext, sample = SAMPLE_IDS, replicate = replicates),
		# expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.{file}.Aligned.out.bam', file = SAMPLE_IDS),
		# expand(working_dir + 'telescope/{sample}_{replicate}-telescope_report.tsv', sample = SAMPLE_IDS, replicate = replicates),
		# expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS),
		# expand(working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt', sample = SAMPLE_IDS)
		working_dir + 'all.samples.telescope.DESeq2.tibble.tsv',
		working_dir + 'all.samples.TEtranscripts.DESeq2.tibble.tsv'


rule FastQC:
	input:
		working_dir + 'raw_files/{sample}_{replicate}.' + raw_file_ext,

	output:
		working_dir + 'FastQC/{sample}_{replicate}_fastqc.html',
		working_dir + 'FastQC/{sample}_{replicate}_fastqc.zip'
	log: 'logs/FastQC.{sample}_{replicate}.log'
	shell:
		'''
		ml fastQC/0.11.8 &&
				fastqc -o FastQC {input}
		'''
rule CutAdapt:
	input:
		read1 = working_dir + 'raw_files/{sample}_{replicate}.' + raw_file_ext,
	output:
		read1 = working_dir + 'cutadapt/{sample}_{replicate}.cutadapt.q20.minlen1.' + raw_file_ext,
		report = working_dir + 'cutadapt/{sample}_{replicate}.cutadapt.report.txt'
	log: 'logs/CutAdapt.{sample}_{replicate}.log'
	shell:
		'''
		ml python/2.7.16

		cutadapt -a {config[fw_adapter]} -q 20 --minimum-length 1 -o {output.read1} {input.read1} > cutadapt/{wildcards.sample}_{wildcards.replicate}.cutadapt.report.txt
		'''
rule FastQC_pass2:
	input:
	  working_dir + 'cutadapt/{sample}_{replicate}.cutadapt.q20.minlen1.' + raw_file_ext

	output:
	  working_dir + 'FastQC_2/{sample}_{replicate}.cutadapt.q20.minlen1_fastqc.html',
	  working_dir + 'FastQC_2/{sample}_{replicate}.cutadapt.q20.minlen1_fastqc.zip'
	log: 'logs/FastQC_2.{sample}_{replicate}.log'
	shell:
		'''
		ml fastQC/0.11.8 &&
			fastqc -o FastQC_2 {input}
		'''

rule STAR:
	input:
		read1 = working_dir + 'cutadapt/{sample}_{replicate}.cutadapt.q20.minlen1.' + raw_file_ext,
	output:
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Aligned.out.bam',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Log.final.out',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Log.out',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Log.progress.out',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.SJ.out.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/exonGeTrInfo.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/exonInfo.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/geneInfo.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/sjdbInfo.txt',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/sjdbList.fromGTF.out.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/sjdbList.out.tab',
		working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}._STARgenome/transcriptInfo.tab'
	log: 'logs/STAR.{sample}_{replicate}.log'
	shell:
		'''
		ml star/2.7.0e
		STAR --runThreadN 16 --genomeDir {config[STAR_genomeDir]} \
		--sjdbGTFfile {config[STAR_GTF]} \
		--sjdbOverhang 100 \
		--readFilesIn {input.read1} \
		--readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix RNAseq.STAR/RNAseq.STAR.{wildcards.sample}_{wildcards.replicate}.
		'''


rule TEtranscripts:
	input:
		treatment_files = expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.{{sample}}_{replicate}.Aligned.out.bam', replicate = replicates),
		control_files = expand(working_dir + 'RNAseq.STAR/RNAseq.STAR.' + controlSample + '_{replicate}.Aligned.out.bam', replicate = replicates)
	output:
		working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt',
		working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq.cntTable',
		working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_DESeq2.R',
		working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_sigdiff_gene_TE.txt'
	log: 'logs/TEtranscripts.{sample}.log'
	conda: 'environment.yaml'
	shell:
		'''
		ml python/2.7.16
		ml R/3.4.4
		#ml gcc/8.1.0 no longer need to re-load gcc after loading R in pegasus
		ml xz/5.2.5
		TEtranscripts --format BAM --mode multi --stranded reverse -t {input.treatment_files} -c {input.control_files} \
		--GTF {config[TEtrx_GTF]} \
		--TE {config[TEtrx_TE]} \
		--project TEtranscripts/{wildcards.sample}.TEtranscripts.DESeq

		'''

rule Telescope:
	input: working_dir + 'RNAseq.STAR/RNAseq.STAR.{sample}_{replicate}.Aligned.out.bam'
	output: working_dir + 'telescope/{sample}_{replicate}-telescope_report.tsv'
	log: 'logs/Telescope.{sample}_{replicate}.log'
	conda: 'environment.yaml'
	shell:
		'''
		telescope assign {input} {config[telescope_GTF]} \
		--outdir telescope/ --exp_tag {wildcards.sample}_{wildcards.replicate}

		'''

rule Telescope_DESeq:
	input:
		treatment_reports = expand(working_dir + 'telescope/{{sample}}_{replicate}-telescope_report.tsv', replicate = replicates),
		control_reports = expand(working_dir + 'telescope/' + controlSample + '_{replicate}-telescope_report.tsv', replicate = replicates),
		script = 'scripts/make.telescope.DESeq2.input.filter.baseMean.10.py',
	output:
		treat_files_list = temp(working_dir + 'telescope/{sample}.treat_files.txt'),
		cntrl_files_list = temp(working_dir + 'telescope/{sample}.cntrl_files.txt'),
		cntTable = working_dir + 'telescope/{sample}.telescope.count.table.tsv',
		DESeq2_script = working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.Rscript.R',
		DESeq2_output = working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv'
	params:
		cntrl_sample = controlSample,
		workingDir = working_dir
	log: 'logs/Telescope_DESeq.{sample}.log'
	conda: 'environment.yaml'
	shell:
		'''
		ml python
		ml R/3.4.4
		#ml gcc/8.1.0 no longer need to re-load gcc after loading R in pegasus

		ls {params.workingDir}telescope/*{wildcards.sample}*telescope_report.tsv > {output.treat_files_list}
		ls {params.workingDir}telescope/*{params.cntrl_sample}*telescope_report.tsv > {output.cntrl_files_list}

		#  make the DESeq count table and generate a DESeq script
		python {input.script} {output.cntrl_files_list} {output.treat_files_list} {config[telescope_GTF]} {params.workingDir}telescope/ -o {wildcards.sample}.telescope.count.table
		# run the generated script
		Rscript {output.DESeq2_script}

		'''
rule make_sample_tables:
	input:
	output:
		telescope_table = working_dir + 'telescope/sample_table.txt',
		TEtranscripts_table = working_dir + 'TEtranscripts/sample_table.txt'
	log: 'logs/make_sample_tables.log'
	run:
		## telescope sample table #####
		# table header
		table_str = "DESeq_output_file" + "\t" + "Sample_name" + "\n"

		for sampleName in SAMPLE_IDS:
			table_str = table_str + sampleName + ".telescope.count.table.DESeq2.tsv" "\t" + sampleName + "\n"

		print("Generated sample_table for telescope. Add columns to this file with optional sample data to be added to the final " +
			"all.samples.telescope.DESeq2.tibble.tsv output")
		print(table_str)
		sample_table_file = open(output.telescope_table, "w")
		sample_table_file.write(table_str)
		################################

		## TEtranscripts sample table ##
		# table header
		table_str = "DESeq_output_file" + "\t" + "Sample_name" + "\n"

		for sampleName in SAMPLE_IDS:
			table_str = table_str + sampleName + ".TEtranscripts.DESeq_gene_TE_analysis.txt" "\t" + sampleName + "\n"

		print("Generated sample_table for TEtranscripts. Add columns to this file with optional sample data to be added to the final " +
			"all.samples.telescope.DESeq2.tibble.tsv output")
		print(table_str)
		sample_table_file = open(output.TEtranscripts_table, "w")
		sample_table_file.write(table_str)
		################################


rule Combine_tables_telescope:
	input:
		telescope_files = expand(working_dir + 'telescope/{sample}.telescope.count.table.DESeq2.tsv', sample = SAMPLE_IDS_no_control),
		sample_table = working_dir + 'telescope/sample_table.txt',
		script = 'scripts/combine.tables.plot.data.py'
	output: working_dir + 'all.samples.telescope.DESeq2.tibble.tsv'
	params:
		workingDir = working_dir
	log: 'logs/Combine_tables_telescope.log'
	shell:
		'''
		ml python
		python {input.script} \
		{params.workingDir}telescope/ \
		-n all.samples.telescope -x .count.table.DESeq2.tsv -a {config[telescope_GTF]}  \
		-p TRUE -s {input.sample_table}
		'''
rule Combine_tables_TEtranscripts:
	input:
		TEtranscripts_files = expand(working_dir + 'TEtranscripts/{sample}.TEtranscripts.DESeq_gene_TE_analysis.txt', sample = SAMPLE_IDS_no_control),
		sample_table = working_dir + 'TEtranscripts/sample_table.txt',
		script = 'scripts/combine.tables.plot.data.py'
	output: working_dir + 'all.samples.TEtranscripts.DESeq2.tibble.tsv'
	params:
		workingDir = working_dir
	log: 'logs/Combine_tables_TEtranscripts.log'
	shell:
		'''
		ml python
		python {input.script} \
		{params.workingDir}TEtranscripts/ \
		-n all.samples.TEtranscripts -x .DESeq_gene_TE_analysis.txt   \
		-p TRUE -s {input.sample_table}
		'''
