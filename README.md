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
