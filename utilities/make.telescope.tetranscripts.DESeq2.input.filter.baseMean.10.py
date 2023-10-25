#!/usr/bin/python
# Last updated: 2023.09.15 by Kevin Nestler
import sys  # Allows you to use system variabls such as ARGV
import os   # Allows us to use paths in our code for file locations
import argparse # This can let you take user input to fine tune variables
import re   # This lets you use regular expressions
import pandas # This lets you join tables together without dictionaries
import textwrap # This lets me clear indentation from the script output I make
# This script will take in 2 files/lists that specify sample file unique IDs:
# 1. Control Samples
# 2. Treatment/Experimental Samples
# 3. Annotation file Path -- Assumes GTF!!! This is the annotation file for the Telescope reports only.
# 4. Output directory
# 5. Optional input will be a tag for output file names and a flag for handling NA values.
# 6. Required input will be a flag for mode: tetranscripts or telescope.

# The input files should provide the absolute file path for each Telescope report file and TEtranscripts count table.
# The path for each Telescope report file should end in "-telescope_report.tsv".
# The path for each TEtranscripts count table should end in "-tetranscripts.cntTable".
# Please make sure the sample names match before the hyphen in the file name! 
# This is what the script will use to match the samples.

# The mode flag gives you the option to generate a single count table and DESeq2 script for:
# 1. TEtranscripts data only -- This will only use the TEtranscripts data to generate a count table and DESeq2 script.  
# 2. Telescope data only -- This will use the canonical gene IDs from the TEtranscripts data and the Telescope data to generate a count table and DESeq2 script.

# Example Telescope report format for an individual file:
# [jimcdonald@mgpc telescope.troubleshoot]$ head TCGA-04-1331-01A.bt2.tele.troubleshoot-telescope_report.tsv
# ## RunInfo      version:1.0.3   annotated_features:28513        total_fragments:93648604        pair_mapped:73382346    pair_mixed:17869347     single_mapped:0 unmapped:2396911        unique:55889111 ambig:35362582  overlap_unique:471738overlap_ambig:1165872
# transcript      transcript_length       final_count     final_conf      final_prop      init_aligned    unique_count    init_best       init_best_random        init_best_avg   init_prop
# __no_feature    0       1077218 1064019.00      0.586   1113337 0       907823  936650  936621.33       0.494
# L1FLnI_20q13.12e        10494   86711   86306.00        0.0172  90696   84418   87714   88011   88011.25        0.0172
# HARLEQUIN_1q32.1        6394    23910   23302.00        0.0233  25386   10822   24131   24285   24270.49        0.0235
# L1FLnI_8q21.11w 10016   9488    9288.00 0.000318        23823   9123    11021   11103   11098.60        0.000349
# L1FLnI_6p22.3e  9889    8056    8048.00 0.00895 9188    7999    8085    8088    8087.32 0.00895

# Example TEtranscripts count table format for an individual file:
# [knestler@log004 raw_files]$ head CNMC_D_760_2-tetranscripts.cntTable
# gene/TE CNMC_D_760_2Aligned.out.bam          ## it ultimately doesn't matter what the sample name is here because it will be replaced with the sample name of the respective file before the hyphen
# "ENSG00000000003.12"    429
# "ENSG00000000005.5"     0
# "ENSG00000000419.10"    745
# "ENSG00000000457.11"    691
# "ENSG00000000460.14"    474
# "ENSG00000000938.10"    0
# "ENSG00000000971.13"    9
# "ENSG00000001036.11"    495
# "ENSG00000001084.8"     809

# This script will go through each Telescope report TSV file and/or TEtranscripts count table
# and compile data for the counts into a single count file for DESeq2. 
# Missing entries will be given a 0 value. It will output the count file 
# as well as a DESeq2 script that is prepared to run on that count file.

# Get user arguments/usage
parser = argparse.ArgumentParser(description="Create a count table to input into DESeq2 from Telescope report files and TEtranscripts count tables. Make a DESeq2 script to process that data.")
parser.add_argument("cntrl_files", type=str, help="File/list containing absolute file paths for control samples.")
parser.add_argument("treat_files", type=str, help="File/list containing absolute file paths for treated or experimental samples.")
parser.add_argument("out_dir", type=str, help="Directory for output count tables and DESeq secripts.")
parser.add_argument("-a", "--annotation", type=str, help="Absolute file path to GTF file used by Telescope during the telescope align step.")
parser.add_argument("-o", "--out_name", type=str, default="combined.count.table", help="Name for output files.")
parser.add_argument("-na", "--NA_value", type=str, default="zero", choices=['zero', 'exclude'], help="How to handle NA values -- set to zero or exclude.")
parser.add_argument("-mode", required=True, type=str, choices=['tetranscripts', 'telescope'], help="Generate DESeq2 R script and count table from tetranscripts and telescope data.")
args = parser.parse_args()

def make_annotation_table(annotation_file):
    # Initialize the data frame
    annotation_frame = pandas.DataFrame(columns = ['transcript'])
    locus_dict = dict()

    # Make a dictionary to find all unique locus IDs from the gene_id field of the annotation file
    # Print the unique ones to the data frame
    with open(annotation_file) as annotation_file:
        for line in annotation_file:
            locus_regex = re.search('gene_id\s?"([a-zA-Z0-9_.\\(\\)]+)";\s?transcript_id', line)
            if locus_regex:
                locus = locus_regex.group(1)
                if locus not in locus_dict:
                    locus_dict[locus] = "found"
                    new_row = pandas.DataFrame({'transcript': [locus]}) # Previous line used append to add a new row to the data frame, but this is deprecated. Now, we use concat instead.
                    annotation_frame = pandas.concat([annotation_frame, new_row], ignore_index=True)
            else:
                # Some lines in Matthew's annotation file give the gene_id for ERVs: "### ERV316A3_1p36.33 ###" on a separate line
                # I will skip these with the regex. However, you can put a warning here if you want instead of just continuing
                #print("Warning! This line lacks proper gene_id format:\n\t", line)
                continue
    return annotation_frame

# Take a new report file and add it to the count table using pandas merge
def join_telescope_tables(telescope_report_filename, data_frame):
    file_base_name = os.path.basename(telescope_report_filename)
    report = pandas.read_csv(telescope_report_filename, sep='\t', skiprows=1)
    slim_report = pandas.DataFrame(report, columns=['transcript', 'final_conf']) # If you want to have user specified column, replace this variable with passed input from args
    slim_report = pandas.DataFrame.rename(slim_report, columns={'final_conf' : file_base_name})
    data_frame = pandas.merge(data_frame, slim_report, on='transcript', how='left')
    return data_frame

# Take a new tetranscripts report file and add it to the count table using pandas merge
def join_tetranscripts_tables(tetranscripts_filename, data_frame):
    file_base_name = os.path.basename(tetranscripts_filename)
    report = pandas.read_csv(tetranscripts_filename, sep='\t', index_col=0)
    report.columns.values[0] = file_base_name
    data_frame = pandas.concat([data_frame, report], axis=1)
    return data_frame

# Merge the telescope and tetranscripts data frames
def join_telescope_tetranscripts_dataframes(telescope_output_data_frame, tetranscripts_output_data_frame, data_frame):
    # reset the index of the tetranscripts data frame
    tetranscripts_output_data_frame = tetranscripts_output_data_frame.reset_index()
    # make sure the genes are strings in the tetranscripts output data frame
    tetranscripts_output_data_frame.iloc[:,0] = tetranscripts_output_data_frame.iloc[:,0].astype(str)
    # keep only the ENSG IDs from the tetranscripts output data frame
    tetranscripts_output_data_frame = tetranscripts_output_data_frame[tetranscripts_output_data_frame.iloc[:,0].str.contains("ENSG")]
    for dataframe in [telescope_output_data_frame, tetranscripts_output_data_frame]:
        dataframe.columns.values[0] = 'transcript'
        dataframe.set_index('transcript', inplace=True)
        dataframe.columns = dataframe.columns.str.rsplit("-", n=1).str[0]
        dataframe.reset_index(inplace=True)
    data_frame = pandas.concat([tetranscripts_output_data_frame, telescope_output_data_frame], axis=0)
    return data_frame

# Generate the R script to run DESeq2 on the combined count tables
def make_Rscript(count_out_path, cntrl_count, treat_count, DESeq2_out_path):
    script_text = \
'''library("DESeq2")
# Load in the data
countdata <- read.table("{count_file}", header=TRUE)
#Then create a data frame matching the sample names and treatments so that the model can be made
names <- colnames(countdata)
treatment <- c(rep("control", {cntrl_count}), rep("experimental", {treat_count}))
coldata <- data.frame(treatment, row.names=names)
# Pre-filtering
countdata <- countdata[rowSums(countdata) >=10,]
# Make the DESeq object
cntTable <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
deseq <- DESeq(cntTable)
# Unless you have replicates, expect a warning that the samples were counted as replicates for purposes
# of determining dispersion, making the differential expression suspect and good only for exploratory purposes.
results <- results(deseq)
# Just to have a record in the ERR/OUT FILES, check that you have the samples labeled as treatment_<treatment name>_vs_<control treatment name>
resultsNames(deseq)
# Make sure this is correct with the relevel command anyway.
cntTable$treatment <- relevel(cntTable$treatment, "control")
# Double check the right output
resultsNames(deseq)
deseq <- DESeq(cntTable)
# Then make the final data for a specific comparison.
results <- results(deseq, contrast = c("treatment", "experimental", "control") )
# make the results a data frame
results <- as.data.frame(results)
# Add the gene names back in
results <- results %>% mutate(name = rownames(results))
# reorder the results columns so that the name column is first
results <- results[,c(7,1,2,3,4,5,6)]
# Take a look at the results
head(results)
# Write the results to a tsv file
write.table(as.data.frame(results), file="{DESeq2_out_path}", sep="\\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        '''.format(count_file = count_out_path.name, cntrl_count = cntrl_count, treat_count = treat_count, DESeq2_out_path = DESeq2_out_path)
    script_dedent = textwrap.dedent(script_text)
    return script_dedent

def main():
    # Prep the files and other variables
    count_out_path = open( "".join([ args.out_dir, args.out_name, ".tsv" ]), "w", newline='')
    script_out_path = open( "".join([ args.out_dir, args.out_name, ".DESeq2.Rscript.R" ]), "w", newline='')
    DESeq2_out_path = "".join([ args.out_dir, args.out_name, ".DESeq2.tsv" ])
    cntrl_file_list = open(args.cntrl_files)
    treat_file_list = open(args.treat_files)
    telescope_report_suffix = 'telescope_report.tsv'
    tetranscripts_report_suffix = 'tetranscripts.cntTable'
    telescope_cntrl_count = 0
    telescope_treat_count = 0
    tetranscripts_cntrl_count = 0
    tetranscripts_treat_count = 0
    files_processed = False

    # If mode is the combined mode
    if args.mode == "telescope":

        # Capture the telescope annotations in a data frame. That I will pass as the base for adding additional data.
        telescope_output_data_frame = make_annotation_table(args.annotation)
        telescope_annotation_length = len(telescope_output_data_frame)
        tetranscripts_output_data_frame = pandas.DataFrame()

        # Process the data files into one table
        for filename in cntrl_file_list:
            filename = filename.rstrip()
            if filename.endswith(telescope_report_suffix):
                telescope_output_data_frame = join_telescope_tables(filename, telescope_output_data_frame)
                telescope_cntrl_count += 1
                files_processed = True
            elif filename.endswith(tetranscripts_report_suffix):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)
                tetranscripts_cntrl_count += 1
                files_processed = True
            continue
        for filename in treat_file_list:
            filename = filename.rstrip()
            if filename.endswith(telescope_report_suffix):
                telescope_output_data_frame = join_telescope_tables(filename, telescope_output_data_frame)
                telescope_treat_count += 1
                files_processed = True
            elif filename.endswith(tetranscripts_report_suffix):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)
                tetranscripts_treat_count += 1
                files_processed = True
            continue

        # Make sure some files were processed and no entries were lost from the data frame
        assert telescope_cntrl_count > 0
        assert telescope_treat_count > 0
        assert tetranscripts_cntrl_count > 0
        assert tetranscripts_treat_count > 0
        assert len(telescope_output_data_frame) == telescope_annotation_length
        assert files_processed == True
        cntrl_count = int((telescope_cntrl_count + tetranscripts_cntrl_count)/2)
        treat_count = int((telescope_treat_count + tetranscripts_treat_count)/2)
        cntrl_file_list.close
        treat_file_list.close

        # Initialize the combined data frame
        combined_data_frame = pandas.DataFrame()
        # Combine the dataframes
        combined_data_frame = join_telescope_tetranscripts_dataframes(telescope_output_data_frame, tetranscripts_output_data_frame, combined_data_frame)
        # Then write the output
        combined_data_frame = combined_data_frame.rename(columns={'transcript' : ''}) # DESeq doesn't want a name on the first column, so remove it
        if args.NA_value == "exclude":
            combined_data_frame = combined_data_frame.dropna()
        combined_data_frame.to_csv(count_out_path, sep='\t', index=False, na_rep=0)

        # Now make the R script and output it
        script_output = make_Rscript(count_out_path, cntrl_count, treat_count, DESeq2_out_path)
        script_out_path.write(script_output)

    # If mode is the tetranscripts mode
    elif args.mode == "tetranscripts":

        # Initialize the tetranscripts output dataframe.
        tetranscripts_output_data_frame = pandas.DataFrame()

        # Process the data files into one table
        for filename in cntrl_file_list:
            filename = filename.rstrip()
            if filename.endswith(tetranscripts_report_suffix):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)
                tetranscripts_cntrl_count += 1
                files_processed = True
            continue
        for filename in treat_file_list:
            filename = filename.rstrip()
            if filename.endswith(tetranscripts_report_suffix):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)
                tetranscripts_treat_count += 1
                files_processed = True
            continue

        # Make sure some files were processed and no entries were lost from the data frame
        assert tetranscripts_cntrl_count > 0
        assert tetranscripts_treat_count > 0
        assert files_processed == True
        cntrl_count = tetranscripts_cntrl_count
        treat_count = tetranscripts_treat_count
        cntrl_file_list.close
        treat_file_list.close

        # Then write the output
        tetranscripts_output_data_frame.reset_index(inplace=True) # This is necessary to get the gene/TE column back
        tetranscripts_output_data_frame.columns.values[0] = '' # DESeq doesn't want a name on the first column, so remove it
        tetranscripts_output_data_frame.columns = tetranscripts_output_data_frame.columns.str.rsplit("-", n=1).str[0]
        if args.NA_value == "exclude":
            tetranscripts_output_data_frame = tetranscripts_output_data_frame.dropna()
        tetranscripts_output_data_frame.to_csv(count_out_path, sep='\t', index=False, na_rep=0)

        # Now make the R script and output it
        script_output = make_Rscript(count_out_path, cntrl_count, treat_count, DESeq2_out_path)
        script_out_path.write(script_output)
    
    # If mode is not recognized
    else:
        print("Error! Mode not recognized. Please use tetranscripts or telescope.")
        sys.exit()

    # Close the files
    count_out_path.close
    script_out_path.close

if __name__ == "__main__":
    main()
