# BEAVR

# Introduction
BEAVR: A **B**rowser-based tool for the **E**xploration and **V**isualization of **R**NAseq data. BEAVR is a browser-based graphical tool to automate analysis and exploration of small and large RNAseq datasets using DESeq2.

# Requirements
You must have the following components installed in order to run BEAVR:
- R 3.5+	
- library(shiny)
- library(shinydashboard)
- library(plotly)
- library(scatterD3)
- library(ggplot2)
- library(ggrepel)
- library(data.table)
- library(DT)
- library("DESeq2")
- library("vsn")
- library('apeglm')
- library('org.Hs.eg.db')
- library('EnhancedVolcano')
- library("gridExtra")
- library("ggpubr")
- library("shinyjqui")

# Installation
As of version 0.62, all required packages should be installed automatically. If you run into an error, try relaunching BEAVR using the commands below.
    
# Run BEAVR
Load the required for launch and run the latest version of BEAVR using:

	library(shiny)
	library(shinydashboard)
	runGitHub( "BEAVR", "developerpiru")

If you are using R, a browser window should open automatically showing the app. If you are using RStudio, click "Open in browser" in the window that opens.

# Usage

BEAVR requires two file inputs:
1. Read counts file
2. Column data file

The **read counts file** contains all the reads for all samples in one file (.txt or .csv file). The format must be as follows:
1. The first column must contain ENSEMBL IDs for every gene. The heading name for this column must be ```gene_id```.
2. The next ```n``` columns must contain the raw read counts for each ```n``` samples. Label the heading name for each column with a unique sample\replicate identifier.

The **column data file** tells BEAVR which columns in the read count file belong to which treatment groups (ie. Untreated and Treated). The format for this must be as follows:
1. The first column must list the sample\replicate identifiers of each sample you have in your read counts file. For example, for ```n``` samples in the read counts file, you must have ```n``` rows in the column data file. Each row is a unique sample. The heading name for this column can be left blank (it is not used).
	- **Note:** it is critical that the order of the samples here (each row) is in the same order as the samples (each column) in the read counts file.
2. The second column identifies which treatment condition\group the samples belong to. The heading name for this column must be ```condition```. For example, in each row of this column, identify that respective sample as belonging to ```Untreated``` or ```Treated```.
3. In the third column, you can specify additional characteristics for each sample. For example, you can specify different genotype groups or replicates like ```Replicate-A```, ```Replicate-B```, and ```Replicate-C```. The heading name for this column must be ```replicate```.

Once you load these files under the ```Load data``` tab and select your experimental settings under the ```Select experiment settings``` tab, click on the ```Output - Differentially expressed gene table``` tab to begin calculations. You will see a progress bar in the bottom right-hand corner of the window. 

Once complete, you can view the resulting data table, PCA plot, read count plots, and volcano plot for your experiment.

The data table can be downloaded and saved using the ```Download Table``` button above the table. You can save the PCA plot or read count plots by right clicking and selecting "save image as...". To save the volcano plot, hover over the plot to display the toolbar in the top-right of the plot. Then click the camera icon to save the current view.

# Known Bugs & Error messages

**Bug:** PCA plot, gene read count plot, and volcano plot won't auto-update after changing the treatment condition. The results table will update correctly. 

**Fix:** Just refresh the page and change the treatment condition.


**Error message:**```ncol(countData) == nrow(colData) is not TRUE```

**Cause:** The samples (columns) in your read count table file do not match the samples (rows) in your sample treatment matrix file.

**Fix:** Please check to make sure sample names match and that you've selected the correct files (and formats) to load into BEAVR.
