# VisualRNAseq

# Introduction
A graphic tool to automate analysis of large RNAseq datasets

# Requirements
You must have the following components installed in order to run VisualizeTRACS:
- R 3.5+	
- library(shiny)
- library(shinydashboard)
- library(plotly)
- library(scatterD3)
- library(ggplot2)
- library(data.table)
- library(DT)
- library("DESeq2")
- library("vsn")
- library('apeglm')
- library('org.Hs.eg.db')

# Installation
Download install the latest version of R if you don't have it already from the CRAN project page: https://cran.r-project.org/.

At the R command prompt or in RStudio, run these commands to install the dependencies (listed above) if you don't already have them installed:

	install.packages("shiny")

	install.packages("shinydashboard")
  
	install.packages("plotly")

	install.packages("scatterD3")
  
  	install.packages("ggplot2")
  
  	install.packages("data.table")

	install.packages("DT")
  
	BiocManager::install("DESeq2")
  
	BiocManager::install("vsn")
  
	BiocManager::install("apeglm")
  
	BiocManager::install("org.Hs.eg.db")
  
Note: if you don't have Bioconductor's BiocManager installed, you can install it with:
  
	install.packages("BiocManager")
  
# Run VisualRNAseq
Load the required libraries:

	library(shiny)
	library(shinydashboard)
 	library(plotly)
 	library(scatterD3)
 	library(ggplot2)
 	library(data.table)
	library(DT)
	library("DESeq2")
	library("vsn")
	library('apeglm')
	library('org.Hs.eg.db')

Then run the latest version of VisualRNAseq using:

	runGitHub( "VisualRNAseq", "developerpiru")

If you are using R, a browser window should open automatically showing the app. If you are using RStudio, click "Open in browser" in the window that opens.

# Usage

VisualRNAseq requires two file inputs:
	1. Read counts file
	2. Column data file

The read counts file contains all the reads for all samples in one file (TXT or CSV file). The format must be as follows:
	1. The first column must contain ENSEMBL IDs for every gene. The heading name for this column must be '''gene_id'''.
	2. The next '''n''' columns must contain the raw read counts for each '''n''' samples. Label the heading name for each column with a unique sample\replicate identifier.

The column data file tells VisualRNAseq which columns in the read count file belong to which treatment groups (ie. Control and Knockout, or Untreated and Treated). The format for this must be as follows:
	1. The first column must list the sample\replicate identifiers of each sample you have in your read counts file. For example, for '''n''' samples in the read counts file, you must have '''n''' rows in the column data file. Each row is a unique sample. The heading name for this column can be left blank (it is not used).
	2. The second contain identifies which condition\group the sample belong to. The heading name for this column must be '''condition'''. For example, in each row of this column, identify that respective sample as belonging to '''Control''' or '''Knockout'''.
	3. The third column must specify whether the sample's sequencing was '''single-read''' or '''paired-end'''. The heading name for this column must be '''type'''.

Once you load these files under the '''Load data''' tab and select your experimental settings under the '''Select experiment settings''' tab, click on '''Output - Differentially expressed gene table''' to begin calculations. You will see a progress bar in the bottom right-hand corner of the window. 

Once complete, you can view the resulting data table, PCA plot, gene count plots, and volcano plot for your experiment.

The data table can be downloaded and saved using the '''Download Table''' button above the table. You can save the PCA plot or gene count plot by right clicking and selecting "save image as...". To save the volcano plot, hove over the plot to display the toolbar in the top-right of the plot. Then click the camera icon to save the current view.

