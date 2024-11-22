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
    parser.add_argument("-mode", required=True, type=str, choices=['tetranscripts', 'telescope'], help="Generate count table from tetranscripts or telescope data.")
    args = parser.parse_args()

    # Initialize dataframes
    telescope_output_data_frame = pandas.DataFrame()
    tetranscripts_output_data_frame = pandas.DataFrame()
    
    # Process input files
    with open(args.input_files) as file_list:
        for filename in file_list:
            filename = filename.rstrip()
            if filename.endswith('TE_counts.tsv'):
                telescope_output_data_frame = join_telescope_tables(filename, telescope_output_data_frame)
            elif filename.endswith('tetranscripts.cntTable'):
                tetranscripts_output_data_frame = join_tetranscripts_tables(filename, tetranscripts_output_data_frame)

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

    else:
        print("Error! Mode not recognized. Please use tetranscripts or telescope.")
        sys.exit()

if __name__ == "__main__":
    main()