#!/usr/bin/python

import sys  # Allows you to use system variabls such as ARGV
import os   # Allows us to use paths in our code for file locations
import argparse # This can let you take user input to fine tune variables
import re   # This lets you use regular expressions
import pandas # This lets you join tables together without dictionaries
import csv # try this way to read tibble for graphs
import textwrap # This lets me clear indentation from the script output I make
import numpy as np # For creating arrays
import matplotlib as mpl # For making some simple graphs
import matplotlib.pyplot as plt # For shortcut use of the scripting interface of matplotlib
from matplotlib.backends.backend_pdf import PdfPages # Allows you to save mutliple figures to one PDF
import datetime # This lets me get the time stamp for files without a name
import inspect # Lets me get line numbers

# This script will take in a single directory of files. It will:
# 1. Identify the DESeq output files.
# 2. Create a tibble for manipulation in R containing data from all files.
#    OPTIONAL INPUT: a sample file. Will allow additional data columns in the tibble.
# 3. NOT IMPLEMENTED YET -- It will output some preliminary violin plots and bar graphs to give a summary of the data.
# OPTIONAL INPUT: annotation file in GTF format. Source column must read: "rmsk", "l1base", or "RMSK_hg38_4.0.5"
# CURRENTLY: Telescope ERVS have "rmsk" and LINEs have "l1base". TEtranscripts annotations have "RMSK_hg38_4.0.5".
# Update this as annotation files evolve.
# OPTIONAL INPUT: A sample file with a column labeled 'DESeq_ouput_file' that has additional columns of information to add as identifiers for that file.
# For the sample file, the file name column should contain the full basename of the file as it appears in los of the data directory. Only exact matches work as currently set up.

# Example DESeq2 Data (DESeq data has a slightly different format). NOTE: the column headers are added by me (not present in original file)
# [jimcdonald@login3 jimcdonald]$ head ~/DESeq2.analysis/TCGA.OVCA.RNAseq.telescope.WTv3.DESeq2.tsv
# id    baseMean    log2FoldChange  lfcSE   stat    pvalue  padj
# "ERV316A3_1p36.33"      0       NA      NA      NA      NA      NA
# "HML3_1p36.33"  0.938356720421734       0.713700354593239       0.907375857934334       0.786554268942087       0.431542838721783       NA
# "MER4B_1p36.33" 123.887075283514        0.443213374837082       0.311939061645902       1.42083319895409        0.155365256758875       0.519873822597489

# Example DESeq Data:
# [jimcdonald@login3 RNAseq.A2780.library001.treat3.AZA]$ head RNAseq.A2780.library001.treat3.AZA.DESeq.try2_gene_TE_analysis.txt
# id      baseMean        baseMeanA       baseMeanB       foldChange      log2FoldChange  pval    padj
# ENSG00000000003.12      16.6964689435888        13.7347000153593        19.6582378718184        1.43128265268517        0.517308606627648       0.618865939400659       1
# ENSG00000000005.5       3.21814641581001        3.81519444871091        2.62109838290912        0.68701567328888        -0.541585082425921      0.987921071670796       1
# ENSG00000000419.10      8.67587176965805        6.86735000767963        10.4843935316365        1.52670149619751        0.610418011019129       0.679332671603292       1

# Usage:
    # Minimum input required for the script to run:
        # python make.TCGA.telescope.DESeq2.tibble.py /home/user/scripts/data
    # Use all available optional arguments (short versions). This will output a file with the name 'example_output_filename' in the directory '/home/user/scripts/data/results'. Data will be analyzed using DESeq version 1. Files ending in DESeq.tsv will be analyzed as the input. Probability values from DESeq input will be included in tibble. Threshold p value will be 0.01. a, g, s ________:
        # python make.TCGA.telescope.DESeq2.tibble.py /home/user/scripts/data -d /home/user/scripts/data/results -n example_output_filename -v 1 -x DESeq.tsv -s ______ -a ______ -p TRUE -t 0.01 -g ____

# Get user arguments/usage
parser = argparse.ArgumentParser(description="Create a tibble from DESeq output and create some graphs for preliminary data visualization.")
parser.add_argument("data_dir", type=str, help="Absolute file path for the directory containing the DESeq data to be analzyed.")
parser.add_argument("-d", "--out_dir", type=str, default="unspecified", help="Directory for tibble and graph output. If not specified, output files are saved in current working directory.")
parser.add_argument("-n", "--out_name", type=str, default="unspecified", help="Name for output files. If not specified, the name will be a timestamp.")
parser.add_argument("-v", "--version", type=int, default="2", choices=[1,2], help="DESeq version: 1 = DESeq; 2 = DESeq2. Default = 2.")
parser.add_argument("-x", "--suffix", type=str, default="DESeq2.tsv", help="Tail end of DESeq2 output files. Used to identify files in a directory as correct input. 'Default = DESeq2.tsv'")
parser.add_argument("-s", "--samples", type=str, default="unspecified", help="Optional file containing sample names and other descriptors to be added to the tibble. MUST INCLUDE: FULL INPUT FILE NAME WITH COLUMN LABEL OF 'DESeq_ouput_file'!! Should include a sample name also. If not used, the input file name will be left as the sample name. If included, sample name is assumed to be a column in the sample file.")
parser.add_argument("-a", "--annotation", type=str, default="unspecified", help="Optional file containing RepeatMasker annotation in BED format. If included, RE class added to tibble.")
parser.add_argument("-p", "--probability", type=str, default="FALSE", choices=['TRUE','FALSE'], help="Include probability from DESeq/DESeq2 output? Options: TRUE or FALSE. Default = FALSE.")
parser.add_argument("-t", "--p_threshold", type=float, default=0.05, help="P-value threshold to use for determining significance. Default = 0.05")
parser.add_argument("-g", "--gtf_feature", type=str, default="gene_id", help="Feature from field in last column of GTF used as identifier in the count table. Should be unique. Default is gene_id.")
args = parser.parse_args()

# Set the log2(fold-chage) column number based on DESeq version number
# See examples above for the output format
def set_fold_change_col():
    print("DESeq version for input: {version}".format(version = args.version))
    log2fc_col = 2 # DESeq2 default has the fold change in the 3rd column
    prob_col = 6 # DESeq2 padj value is in the 7th column
    if args.version == 1:
        log2fc_col = 4
        prob_col = 7
    print("Fold change data are in column", log2fc_col+1)
    return log2fc_col, prob_col

# Set the names for output files
def make_output_names(filename, timestamp):
    out_dir = os.getcwd()
    if args.out_dir != "unspecified":
        out_dir = args.out_dir
    print("Output Directory: {dir_name}".format(dir_name = out_dir))

    split_filename = os.path.splitext(filename)
    out_prefix = split_filename[0] # add input filename to tibble name to prevent overwriting
    if args.out_name != "unspecified":
        out_prefix = args.out_name
    # tibble_name = "".join([ out_dir, "/", out_prefix, ".DESeq", str(args.version), ".tibble.tsv" ])
    tibble_name = "".join([ out_dir, "/", out_prefix, ".DESeq", str(args.version), '.', timestamp, ".tibble.tsv" ])
    summary_data_name = "".join([ out_dir, "/", out_prefix, ".DESeq", str(args.version), '.', timestamp, ".summary_data.pdf" ])

    # out_dir_files = os.listdir(out_dir)
    # if os.path.basename(tibble_name) in out_dir_files:
    #     print("WARNING!! Tibble output file already exists in output directory. Adding time stamp to prevent overwriting data.")
    #     tibble_name = "".join([ out_dir, "/", timestamp, ".", out_prefix, ".DESeq", str(args.version), ".tibble.tsv" ])

    print("Tibble output to: {tibblepath}".format(tibblepath = tibble_name))
    print("Summary data output to: {datapath}".format(datapath = summary_data_name))
    return [out_dir, tibble_name, summary_data_name]

# If user specified an annotation file, get the data from it.
def capture_annotation_data(annotation_file):
    annotation_class_dict = dict()
    line_count = 0
    target = args.gtf_feature
    with open(annotation_file) as annotation_file:
        for line in annotation_file:
            line_count += 1
            line_array = line.strip().split("\t")
            locus_id = ""
            re_class = ""

            # Current Telescope annotation files have each new locus given a section label such as '### ERV316A3_1p36.33 ###'.
            # This will fail when checking line_array[2]. Put here to also catch any other files with lines not split properly.
            try:
                line_array[1]
            except IndexError:
                print("WARNING!! Annotation file line {line_number} did not split properly and may not be in GTF format. Consider checking GTF file formating.".format(line_number = line_count))
                continue

            # Now get the locus ID and RE class
            if line_array[0].startswith("chr"):
                # Get the locus ID
                locus_regex = re.search('%s\s?"([a-zA-Z0-9_.\\(\\)-\/]+)";\s?\w+' % target, line_array[8])
                if locus_regex:
                    locus_id = locus_regex.group(1)
                else:
                    print("WARNING!! Annotation file line {line_number} does not contain a properly formatted locus ID. Skipping!!".format(line_number = line_count))
                    continue

                # Get the class
                if line_array[1] == "rmsk": # Telescope annotation for LTRs
                    re_class = "LTR"
                elif line_array[1] == "l1base": # Telescope annotation for LINEs
                    re_class = "LINE"
                elif line_array[1] == "RMSK_hg38_4.0.5": # Full hg38 RepeatMasker annotation GTF file James made has this. Use for TEtranscripts
                    class_regex = re.search('class_id\s?"([a-zA-Z0-9_.?]+)";$', line_array[8])
                    if class_regex:
                        re_class = class_regex.group(1)
                else:
                    print("WARNING!! Annotation file line {line_number} does not match known options for source column. Skipping!!".format(line_number = line_count))
                    continue
            else:
                print("WARNING!! Annotation file line {line_number} does not begin with a chromosome name. Skipping!!".format(line_number = line_count))
                continue

            # Check that the values were defined
            if locus_id == "" or re_class == "": # They should both be defined or skipped. Just double checking for cases I didn't think of.
                print("WARNING!! Annotation file line {line_number} did not define either locus ID or RE class. Skipping!!".format(line_number = line_count))
                continue

            # Now check the dictionary and update dictionary if it isn't in there
            if locus_id not in annotation_class_dict:
                annotation_class_dict[locus_id] = re_class
            else:
                continue
    return annotation_class_dict

def process_sample_file(sample_file):
    sample_df = pandas.read_csv(sample_file, sep='\t')
    if 'DESeq_output_file' not in sample_df.columns:
        print('Sample file provided by user lacks a column "DESeq_output_file" that designates output files. Cannot use this sample file.')
        use_samples = 0
    else:
        print("Head of sample Data Frame:")
        print( sample_df.head() )
        print("\n")
        use_samples = 1
    return sample_df, use_samples

# Take a new report file and add it to the count table using pandas merge
def make_tibble(input_file, sample_df, use_samples, log2fc_col, tibble_handle, re_class_dict, probability_col):
    # Get the index of the data for this file in the sample data frame and check, also checking if it exists in the data frame
    if use_samples == 1: # This passed the check for the file column so check if the file is listed. Requires exact file name match
        input_df_index = list(sample_df[sample_df['DESeq_output_file'].isin([input_file]) == True].index)
        if len(input_df_index) == 0:
            print("WARNING!! The file {filename} does not have an exact match in the sample file. No additional information added for this file.".format(filename = input_file))
            use_samples == 2 # This will trigger addition of NA values to the tibble for these fields so the table still matches with samples that have the information

    # Get each line in the DESeq output, add relevant data, and print to the tibble
    absolute_input_file = "/".join([ args.data_dir, input_file ])
    with open(absolute_input_file) as input_file:
        for line in input_file:
            if 'baseMean' in line or '__no_feature' in line: # Filter out header and one feature I don't want counted from Telescope
                continue
            line_array = line.strip().split("\t")
            log2fc = line_array[log2fc_col]
            locus_id = line_array[0]
            # Telescope output for the gene/transcript ID is in quotes. I need to remove to match regex above.
            # Can't just widen regex as TEtranscripts doesn't have quotes (I checked).
            if locus_id.startswith('"'):
                locus_id = locus_id[1:-1] # Telescope output for the gene/transcript ID is in quotes. I need to remove to match regex above. Can't just widen regex as TEtranscripts doesn't have quotes (I checked).
            # TEtranscripts DESeq/DESeq2 output has REs with names such as (A)n:Simple_repeat:Simple_repeat. Get just the first part.
            TEtrx_annotations = re.match('^(.*):(.*):(.*)', line_array[0])
            if TEtrx_annotations:
                locus_id = TEtrx_annotations.group(1)
            output_line = "\t".join([ input_file.name, locus_id, log2fc ])
            # Add user optional values
            if args.probability == "TRUE":
                output_line = "\t".join([ output_line, line_array[probability_col] ])
            if args.annotation != "unspecified":
                output_line = "\t".join([ output_line, re_class_dict[locus_id] ])
            # Then add in labels that can be used as factors in R to organize graphs by differential expression or significance
            if log2fc == 'NA':
                output_line = "\t".join([ output_line, 'NA' ])
            elif float(log2fc) > 1:
                output_line = "\t".join([ output_line, 'UP' ])
            elif float(log2fc) < -1:
                output_line = "\t".join([ output_line, 'DOWN' ])
            else:
                output_line = "\t".join([ output_line, 'no_change' ])
            if line_array[probability_col] == 'NA':
                output_line = "\t".join([ output_line, 'NA' ])
            elif float(line_array[probability_col]) <= args.p_threshold:
                output_line = "\t".join([ output_line, 'sig' ])
            else:
                output_line = "\t".join([ output_line, 'not_sig' ])
            # Check for the sample file and add the appropriate data as needed.
            if use_samples > 0:
                first_tab_index = output_line.find('\t') # Assuming sample name in file (see help above). So get index of first tab so I can remove filename
                output_line = output_line[first_tab_index+1:] # Then actually remove the file name from the output
                output_line = add_sample_info_to_ouput(output_line, input_file, sample_df, input_df_index, use_samples) # Then add the data from the sample file
            output_line = "".join([ output_line, '\n' ])
            tibble_handle.write(output_line)
    print('Tibble construction for current sample complete. Continuing')

def add_sample_info_to_ouput(output_line, input_file, sample_df, input_df_index, use_samples):
    for column in sample_df.columns:
        if column == 'DESeq_output_file': # Skip the file name I just erased
            continue
        elif use_samples == 2: # If no file name match, add NA for each column (doing it this way accounts for variable amounts of columns)
            output_line = "\t".join([ output_line, 'NA' ])
        else:
            col_string = str(list(sample_df[column][input_df_index]))
            col_string = col_string[2:-2]
            output_line = "\t".join([ output_line, col_string ])
    return output_line

def make_graphs(tibble_file, summary_data_name):
    # Open and read tibble file
    print("In graph subroutine. Tibble file is: ", tibble_file)
    tibble = pandas.read_csv(tibble_file, sep='\t', header=None)

    # Set tibble column index based on what arguments were specified since tibble doesn't have column headers
    if args.probability == "FALSE" and args.annotation == "unspecified":
        DE_trx_col = 3
    elif args.probability == "TRUE" and args.annotation == "unspecified":
        DE_trx_col = 4
    elif args.probability == "FALSE" and args.annotation != "unspecified":
        DE_trx_col = 4
    else:
        DE_trx_col = 5
    DE_transcripts = tibble.iloc[:, [DE_trx_col, DE_trx_col-1]] # dataframe with column 1 = UP, DOWN, or no_change and column 2 = log2fc values

    # Get number of differentially expressed transcripts, number of upregulated transcripts, and number of downregulated transcripts
    n_total_DE_transcripts = len(DE_transcripts[DE_transcripts.iloc[:,0].isin(['UP', 'DOWN'])])
    n_upreg = len(DE_transcripts[DE_transcripts.iloc[:,0].isin(['UP'])])
    n_downreg = len(DE_transcripts[DE_transcripts.iloc[:,0].isin(['DOWN'])])

    # create PDF to save plots in
    with PdfPages(summary_data_name) as pdf:
        ### Create bar chart
        fig, ax = plt.subplots()

        bar_descriptor = ['Upregulated', 'Downregulated', 'Total']
        bar_heights = [n_upreg, n_downreg, n_total_DE_transcripts]
        width = 0.5 # width of bars
        # Make plot
        plt.bar(range(len(bar_descriptor)), bar_heights, color=['red', 'dodgerblue', 'darkviolet'], align='center')
        plt.title('Differentially Expressed Transcripts', fontweight='bold', fontsize='14')
        tibble_basename = os.path.basename(tibble_file) # get basename of tibble file to add to graph
        plt.xlabel('Data from file: {f}'.format(f=tibble_basename), fontsize=8)
        plt.ylabel('Number of Transcripts')
        ax.xaxis.set_tick_params(labelsize=12)
        plt.xticks(range(len(bar_descriptor)), bar_descriptor)
        # Annotate bars with heights
        for rect,i in zip(ax.patches, bar_heights):
            ax.annotate('{:,}'.format(i),
                            xy=(rect.get_x() + rect.get_width() / 2, rect.get_height()),
                            xytext=(0, 2),  # 2 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom', size='8', fontweight='bold')
        plt.tight_layout()
        pdf.savefig()

        # Get log2fc values for each locus ID in DE_transcripts (contains transcripts only if they were up or downregulated)
        log2fc_values = DE_transcripts.iloc[:,1].astype(float)
        log2fc_values.dropna(inplace=True) # Get rid of NaNs or else plot won't work

        ### Create violin plot
        fig2, ax1 = plt.subplots()

        violin = ax1.violinplot(log2fc_values, showmeans=False, showmedians=False, showextrema=False)
        ax1.set_title('Expression Fold Change', fontweight='bold', fontsize='12')
        ax1.set_xlabel('Data from file: {x}'.format(x=tibble_basename), fontsize=8)
        ax1.set_ylabel('Log2(Fold Change)')
        ax1.set_xticks([])

        # format body of violin
        for v in violin['bodies']:
            v.set_facecolor('deepskyblue')
            v.set_edgecolor('black')
            v.set_alpha(0.9)

        # get precentile values and plot box b/w first and third quartiles with white dot on median
        quartile1, median, quartile3 = np.percentile(log2fc_values, [25, 50, 75])
        ax1.scatter(1, median, marker='o', color='white', s=30, zorder=3)
        ax1.vlines(1, quartile1, quartile3, color='k', linestyle='-', lw=7)

        # add thinner line using IQR to mark outliers
        upper_bound = quartile3 + (quartile3 - quartile1) * 1.5
        lower_bound = quartile1 - (quartile3 - quartile1) * 1.5
        ax1.vlines(1, upper_bound, lower_bound, color='k', linestyle='-', lw=1)

        # add text box to label percentile values
        text = '\n'.join(('Minumim = {min:.2f}'.format(min=log2fc_values.min()),
            '$1^{{st}}$ quartile = {q1:.2f}'.format(q1=quartile1),
            'Median = {m:.2f}'.format(m=median),
            '$3^{{rd}}$ quartile = {q3:.2f}'.format(q3=quartile3),
            'Maximum = {max:.2f}'.format(max=log2fc_values.max())))
        box_props = dict(boxstyle='square', facecolor='lightgray')
        ax1.text(0.75, 5.25, text, bbox=box_props)

        plt.tight_layout()
        pdf.savefig()

        plt.close()

    return(summary_data_name)

def main():
    # Get the date and time to use as a unique time stamp for files.
    # I realize the alias part is an extra line of code, but I thought it a neat trick and wanted to use it to try and remember it
    date=datetime.datetime.now()
    timestamp = date.strftime("%Y-%m-%d.%H:%M")

    # Prepare the variables and file names ahead of running the analysis
    log2fc_col, prob_col = set_fold_change_col()
    # out_dir, tibble_name, summary_data_name = make_output_names(timestamp)
    # tibble_handle = open(tibble_name, 'w')

    # If the user specified RepeatMasker annotation file (BEF format), get the class annotation
    re_classes = dict()
    if args.annotation != "unspecified":
        re_classes = capture_annotation_data(args.annotation)

    # Process the sample file
    sample_df = pandas.DataFrame()
    use_samples = 0
    if args.samples != "unspecified":
        sample_df, use_samples = process_sample_file(args.samples) # If data frame is OK, use_samples = 1. Otherwise use_samples = 0.
    else:
        print("No sample file specified by user. Sample name in tibble will be input file name.")
        use_samples = 0

    # Process the DESeq data in the folder and make a tibble data frame.
    # Each sample will create a data frame and then append it to the file.
    file_list = os.listdir(args.data_dir)
    for filename in file_list:
        if filename.endswith(args.suffix):
            out_dir, tibble_name, summary_data_name = make_output_names(filename, timestamp)
            tibble_handle = open(tibble_name, 'w')
            print("Processing: ", filename)
            make_tibble(filename, sample_df, use_samples, log2fc_col, tibble_handle, re_classes, prob_col)
            tibble_handle.close() # close tibble handle before using it to make graphs
            print("Graphing data from: ", tibble_name)
            make_graphs(tibble_name, summary_data_name)
        else:
            print(filename, "is not a DESeq output file. Skipping")
            continue

    # Goal is to use matplotlib to produce a violin plot of expression fold change for each sample.
    # Also to figure out the total number of DE transcripts and output the number up, number down, and total.
    # Would be nice if it could plot the up, down, and total counts for significant genes when applicable.
    # But that will have to wait. Depends on how long testing and other projects take.

if __name__ == "__main__":
    main()
