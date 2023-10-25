# Utilities
Scripts and tools used in RNA-seq analyses.
Please add to this readme file whenever new scripts are uploaded.

## _make.telescope.tetranscripts.DESeq2.input.filter.baseMean.10.py_
* This python script was originally created by James McDonald (I think) to generate a single combined count table and a DESeq2.R script from multiple telescope report files.
* This is now updated to generate a single combined count table and DESeq2.R script from Telescope report files or TEtranscripts count tables.
* This script will take the following arguments:
  - Control Samples (.txt file containing full path for each telescope report file and/or tetranscripts count table).
  - Treatment/Experimental Samples (.txt file containing full path for each telescope report file and/or tetranscripts count table).
  - Output directory
  - Optional input for the path to the Telescope annotation file (assumes GTF). This is required when running telescope mode.
  - Optional input will be a flag for output file names.
  - Optinal input will be a flag for handling NA values: zero or exclude. Default is zero.
  - Required input will be a flag for mode: tetranscripts or telescope.
* The input files should provide a txt file containing the absolute filepath for each Telescope report file and TEtranscripts count table.
* The path for each Telescope report file should end in "-telescope_report.tsv".
* The path for each TEtranscripts count table should end in "-tetranscripts.cntTable".
* Please make sure the sample names match before the hyphen in the file name! This is what the script will use to match the samples.
* To run the script:
```
python make.telescope.tetranscripts.DESeq2.input.filter.baseMean.10.py <control samples path> <treat samples path> <output directory> -a <annotation file path> -o <output name> -na <zero, exclude> -mode <tetranscripts, telescope>
```

## New tool
Descrtiption of new tool
