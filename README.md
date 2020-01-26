# BEAVR

BEAVR: A **B**rowser-based tool for the **E**xploration **A**nd **V**isualization of **R**NAseq data. 

BEAVR is a graphical tool to automate analysis and exploration of small and large RNAseq datasets using DESeq2.

# Release

Current stable release is version 0.75.2

# Installation & Requirements

You must have the following components installed in order to run BEAVR:
- R 3.5+	
- library("BiocManager")
- library("shiny")
- library("shinydashboard")

To install the required packages, enter these commands in R:

	install.packages("BiocManager")
	install.packages("shiny")
	install.packages("shinydashboard")

As of version 0.62, all other required packages will be installed automatically. If you run into an error, make sure you have R 3.5+ installed and try reinstalling the above packages. 

# Run BEAVR

Load the required libraries for launch and run the latest version of BEAVR using:

	library(shiny)
	library(shinydashboard)
	runGitHub( "BEAVR", "developerpiru")

If you are using R, a browser window should open automatically showing the app. If you are using RStudio, click "Open in browser" in the window that opens.

# Usage

BEAVR requires two file inputs:
1. Read counts file
2. Column data file

See the **Examples** folder for examples of these two files prepared for the Sehrawat *et al.* (2018) dataset. 

## Preparing the read count table file

The **read counts table file** contains all the raw reads for all the samples in your experiment in a tab-delimited (.txt) or comma-separated (.csv) file type.

The table must be arranged as follows:
1. The first column must contain ENSEMBL IDs for every gene. The heading name for this column must be ```gene_id```.
2. The next ```n``` columns must contain the raw read counts for each ```n``` samples. Label the heading name for each column with a unique sample/replicate identifier.

![Image of read count table](images/readcounttable.jpg)

The **sample treatment matrix file** tells BEAVR which columns in the read count file belong to which treatment groups (ie. Untreated and Treated). The format for this must be as follows:
1. The first column must list the sample\replicate identifiers of each sample you have in your read counts file. For example, for ```n``` samples in the read counts file, you must have ```n``` rows in the column data file. Each row is a unique sample. The heading name for this column can be left blank (it is not used).
	- **Note:** it is critical that the order of the samples here (each row) is in the same order as the samples (each column) in the read counts file.
2. The second column identifies which treatment condition\group the samples belong to. The heading name for this column must be ```condition```. For example, in each row of this column, identify that respective sample as belonging to ```Untreated``` or ```Treated```.
3. In the third column, you can specify additional characteristics for each sample. For example, you can specify different genotype groups or replicates like ```Replicate-A```, ```Replicate-B```, and ```Replicate-C```. The heading name for this column must be ```replicate```.

Once you load these files under the ```Load data``` tab and select your experimental settings under the ```Settings``` tab, click on the ```Gene Table``` tab to begin calculations. You will see a progress bar in the bottom right-hand corner of the window. 

Once complete, you can view the resulting data table, PCA plot, sample clustering/correlation plot, read count plots, heatmap, and volcano plot for your experiment.

The data table can be downloaded and saved using the ```Download``` button in the sidebar. All graphs and plots can be saved by right clicking on the image and saving as an image.

# Known Bugs & Error messages

 Issue | Solution
---------|---------
**Bug:** PCA plot, gene read count plot, and volcano plot won't auto-update after changing the treatment condition. The results table will update correctly. | Edit the parameters in the sidebar for PCA/read count/volcano plots.
**Error message:**```ncol(countData) == nrow(colData) is not TRUE``` | The samples (columns) in your read count table file do not match the samples (rows) in your sample treatment matrix file. Please check to make sure sample names match and that you've selected the correct files (and formats: .csv or .txt) to load into BEAVR.
**Error message:**```None of the keys entered are valid keys for 'ENSEMBL'. Please use the keys method to see a listing of valid arguments``` | This means the ENSEMBL IDs contained in your read count table file cannot be mapped to the reference genome you selected in the Experiment settings tab. Please verify you have selected the correct one.

