import sys
import os
import argparse
import pandas

def join_telescope_tables(telescope_report_filename, data_frame):
    file_base_name = os.path.basename(telescope_report_filename)
    report = pandas.read_csv(telescope_report_filename, sep='\t', index_col=0)
    report.columns.values[0] = file_base_name
    data_frame = pandas.concat([data_frame, report], axis=1)
    return data_frame

def join_tetranscripts_tables(tetranscripts_filename, data_frame):
    file_base_name = os.path.basename(tetranscripts_filename)
    report = pandas.read_csv(tetranscripts_filename, sep='\t', index_col=0)
    report.columns.values[0] = file_base_name
    data_frame = pandas.concat([data_frame, report], axis=1)
    return data_frame

def join_telocal_tables(telocal_filename, data_frame):
    file_base_name = os.path.basename(telocal_filename)
    report = pandas.read_csv(telocal_filename, sep='\t', index_col=0)
    report.columns.values[0] = file_base_name
    data_frame = pandas.concat([data_frame, report], axis=1)
    return data_frame

def annotate_telocal_outputs(annotation_file_path, data_frame):
    telocal_annt = pandas.read_csv(annotation_file_path, sep='\t')
    telocal_annt[['Chromosome', 'Coordinates', 'Strand']] = telocal_annt['chromosome:start-stop:strand'].str.split(':', expand=True)
    telocal_annt[['Start', 'Stop']] = telocal_annt['Coordinates'].str.split('-', expand=True)
    telocal_annt = telocal_annt.drop(columns=['chromosome:start-stop:strand', 'Coordinates'])
    data_frame = telocal_annt
    return data_frame

def join_telescope_tetranscripts_dataframes(telescope_output_data_frame, tetranscripts_output_data_frame, data_frame):
    tetranscripts_output_data_frame = tetranscripts_output_data_frame.reset_index()
    telescope_output_data_frame = telescope_output_data_frame.reset_index()
    tetranscripts_output_data_frame.iloc[:,0] = tetranscripts_output_data_frame.iloc[:,0].astype(str)
    tetranscripts_output_data_frame = tetranscripts_output_data_frame[tetranscripts_output_data_frame.iloc[:,0].str.contains("ENS")]
    for dataframe in [telescope_output_data_frame, tetranscripts_output_data_frame]:
        dataframe.columns.values[0] = 'transcript'
        dataframe.set_index('transcript', inplace=True)
        dataframe.columns = dataframe.columns.str.rsplit("-", n=1).str[0]
        dataframe.reset_index(inplace=True)
    data_frame = pandas.concat([tetranscripts_output_data_frame, telescope_output_data_frame], axis=0)
    return data_frame

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Create combined count tables from Telescope report files and TEtranscripts count tables.")
    parser.add_argument("input_files", type=str, help="File containing list of absolute file paths for all samples.")
    parser.add_argument("out_dir", type=str, help="Directory for output count tables.")
    parser.add_argument("-na", "--NA_value", type=str, default="zero", choices=['zero', 'exclude'], help="How to handle NA values -- set to zero or exclude.")
    parser.add_argument("-mode", required=True, type=str, choices=['tetranscripts', 'telescope', 'telocal'], help="Generate count table from tetranscripts or telescope data.")
    args = parser.parse_args()

    # Initialize dataframes
    telescope_output_data_frame = pandas.DataFrame()
    tetranscripts_output_data_frame = pandas.DataFrame()
    telocal_output_data_frame = pandas.DataFrame()
    telocal_annt_output_data_frame = pandas.DataFrame()
    
    # Process input files
    with open(args.input_files) as file_list:
        for filename in file_list:
            filename = filename.rstrip()
            if filename.endswith('TE_counts.tsv'):
                telescope_output_data_frame = join_telescope_tables(filename, telescope_output_data_frame)
            elif filename.endswith('tetranscripts.cntTable'):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)
            elif filename.endswith('telocal.cntTable'):
                telocal_output_data_frame = join_telocal_tables(filename, telocal_output_data_frame)

    # Process and write output based on mode
    if args.mode == "telescope":
        output_path = os.path.join(args.out_dir, "telescope_counts.tsv")
        # Initialize the combined data frame
        combined_data_frame = pandas.DataFrame()
        # Combine the dataframes
        combined_data_frame = join_telescope_tetranscripts_dataframes(telescope_output_data_frame, tetranscripts_output_data_frame, combined_data_frame)
        # Process and write output
        combined_data_frame = combined_data_frame.rename(columns={'transcript' : ''})
        if args.NA_value == "exclude":
            combined_data_frame = combined_data_frame.dropna()
        combined_data_frame.to_csv(output_path, sep='\t', index=False, na_rep=0)

    elif args.mode == "tetranscripts":
        output_path = os.path.join(args.out_dir, "tetranscripts_counts.tsv")
        # Process tetranscripts output
        tetranscripts_output_data_frame.reset_index(inplace=True)
        tetranscripts_output_data_frame.columns.values[0] = ''
        tetranscripts_output_data_frame.columns = tetranscripts_output_data_frame.columns.str.rsplit("-", n=1).str[0]
        if args.NA_value == "exclude":
            tetranscripts_output_data_frame = tetranscripts_output_data_frame.dropna()
        tetranscripts_output_data_frame.to_csv(output_path, sep='\t', index=False, na_rep=0)
    
    elif args.mode == "telocal":
        output_path = os.path.join(args.out_dir, "telocal_counts.tsv")
        annt_output_path = os.path.join(args.out_dir, "telocal_counts_annotated.tsv")
        # Process telocal output
        telocal_output_data_frame.reset_index(inplace=True)
        telocal_output_data_frame.columns.values[0] = ''
        telocal_output_data_frame.columns = telocal_output_data_frame.columns.str.rsplit("-", n=1).str[0]
        if args.NA_value == "exclude":
            telocal_output_data_frame = telocal_output_data_frame.dropna()
        telocal_output_data_frame.to_csv(output_path, sep='\t', index=False, na_rep=0)

        # Annotate telocal outputs
        annotation_file_path = (r'/SMHS/groups/chiappinellilab/genomes/hg38/RepeatMasker/rmsk_hg38_TE_local_locations.locInd')
        telocal_annt_output_data_frame = annotate_telocal_outputs(annotation_file_path, telocal_annt_output_data_frame)
        telocal_output_data_frame.rename(columns={telocal_output_data_frame.columns[0]: "gene/TE"}, inplace=True)
        telocal_output_data_frame[['TE', 'Subfamily', 'Family', 'Class']] = telocal_output_data_frame['gene/TE'].str.split(":", expand=True)
        telocal_output_data_frame.loc[telocal_output_data_frame["gene/TE"].str.contains("ENS", na=False), "TE"] = None
        telocal_output_data_frame = telocal_output_data_frame.merge(telocal_annt_output_data_frame[['TE', 'Chromosome', 'Strand', 'Start', 'Stop']],
                                                                    on='TE', how='left')
        telocal_output_data_frame["Start"] = telocal_output_data_frame["Start"].astype("Int64")
        telocal_output_data_frame["Stop"] = telocal_output_data_frame["Stop"].astype("Int64")
        telocal_output_data_frame.to_csv(annt_output_path, sep='\t', index=False)

    else:
        print("Error! Mode not recognized. Please use tetranscripts, telescope, or telocal.")
        sys.exit()

if __name__ == "__main__":
    main()