# BEAVR
# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# See Github for more info & ReadMe: https://github.com/developerpiru/VisualRNAseq
app_version = "0.71.4"

# added:
# +1 to all reads; avoid 0 read count errors
# multiple comparisons
# show >2 conditions on PCA plot
# adjust results based on different conditions
# options for plots to show labels
# select either sample names or replicate names for labels
# boxplot or jitter plot option for read count plot
# use ggrepel for plot labels so labels don't overlap
# automatically install required packages if not already installed
# toggle for biocmanager packages
# fixed volcano plot
# added ability to customize font sizes and point sizes for all graphs/plots
# added ability to plot multiple read count plots at once
# customize legend positions on multiple read count plots
# drag to customize the area of all plots
# option to show y-axis title only on first plot per row
# option to show log10 scale y-axis
# dropped single read plot feature - use multi version with 1x1 grid for a single plot

# bugs"
#### PCA, gene count, volcano plots don't auto-update to new dds dataset after changing treatment condition factor level
#### legend symbols show letter 'a' below symbol on jitter plots
#### alignment of legend is not centered on plots but on canvas

#function to check for required packages and install them if not already installed
installReqs <- function(package_name, bioc){
  if (requireNamespace(package_name, quietly = TRUE) == FALSE) {
    if (bioc == FALSE)
      install.packages(package_name)
    else if (bioc == TRUE) #install using Bioconductor package manager
      BiocManager::install(package_name)
  }
}

#check if required libraries are installed, and install them if needed
installReqs("shiny", bioc = FALSE)
installReqs("shinydashboard", bioc = FALSE)
installReqs("plotly", bioc = FALSE)
installReqs("ggplot2", bioc = FALSE)
installReqs("ggrepel", bioc = FALSE)
installReqs("data.table", bioc = FALSE)
installReqs("DT", bioc = FALSE)
installReqs("BiocManager", bioc = FALSE)
installReqs("DESeq2", bioc = TRUE)
installReqs('apeglm', bioc = TRUE)
installReqs('org.Hs.eg.db', bioc = TRUE)
installReqs('EnhancedVolcano', bioc = TRUE)
installReqs('ggpubr', bioc = FALSE)
installReqs('shinyjqui', bioc = FALSE)

#load required libraries
library("shiny")
library("shinydashboard")
library("plotly")
library("ggplot2")
library("ggrepel")
library("data.table")
library("DT")
library("DESeq2")
library('apeglm')
library('org.Hs.eg.db')
library('EnhancedVolcano')
library("gridExtra")
library("ggpubr")
library("shinyjqui")
library("scales")

ui <- dashboardPage(
  dashboardHeader(title = paste("BEAVR", app_version, sep = " ")),
  
  dashboardSidebar(
    
    #Conditional Panels to show tab-specific settings in sidebar
    
    conditionalPanel("input.navigationTabs == 'loadDataTab'",
                     div(id = 'loadDataTab_sidebar', 
                         h4("Welcome to BEAVR!"),
                         #h3("Please see the GitHub page for help & info."),
                         tags$a(href="https://github.com/developerpiru/VisualRNAseq",target="_blank","Check GitHub for help & info")
                     )),
    
    conditionalPanel("input.navigationTabs == 'expSettingsTab'",
                     div(id = 'expSettingsTab_sidebar' 
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'geneTableTab'",
                     div(id = 'geneTableTab_sidebar',
                         h4("Download Table"),
                         downloadButton("downloadDEGeneTable", ""),
                         
                         h4("Filter results table"),
                         #baseMean
                         uiOutput("baseMean_min"),
                         uiOutput("baseMean_max"),
                         #log2FoldChange
                         uiOutput("log2FC_min"),
                         uiOutput("log2FC_max"),
                         #lfcSE
                         uiOutput("lfcSE_min"),
                         uiOutput("lfcSE_max"),
                         #stat
                         uiOutput("stat_min"),
                         uiOutput("stat_max"),
                         #pvalue
                         uiOutput("pvalue_min"),
                         uiOutput("pvalue_max"),
                         #padj
                         uiOutput("padj_min"),
                         uiOutput("padj_max")
                     )),
    
    conditionalPanel("input.navigationTabs == 'pcaPlotTab'",
                     div(id = 'pcaPlotTab_sidebar', 
                         
                         h4("Appearance"),
                         selectInput("PCAplot_labels", label = "Sample labels",
                                     choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3),
                                     selected = 1),
                         numericInput("pcaPointSize", label = "Point size", value = 3),
                         
                         h4("Font sizes"),
                         numericInput("pcaLabelFontSize", label = "Sample labels", value = 5),
                         numericInput("pcaFontSize_xy_axis", label = "Axis labels", value = 18),
                         numericInput("pcaFontSize_legend_title", label = "Legend title", value = 16),
                         numericInput("pcaFontSize_legend_text", label = "Legend labels", value = 15)
                     )),
    
    conditionalPanel("input.navigationTabs == 'genecountPlotTab'",
                     div(id = 'genecountPlotTab_sidebar',
                         
                         textInput("gene_name", "Enter gene name", value = "KRAS"),
                         
                         h4("Appearance"),
                         selectInput("readcountplot_type", label = "Plot type",
                                     choices = list("Boxplot" = 1, "Jitter plot" = 2), 
                                     selected = 1),
                         numericInput("genecountPointSize", label = "Jitter point size", value = 3),
                         selectInput("readcountplot_labels", label = "Sample labels",
                                     choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3), 
                                     selected = 1),
                         h4("Font sizes"),
                         numericInput("genecountFontSize_plot_title", label = "Gene name", value = 20),
                         numericInput("genecountLabelFontSize", label = "Sample labels", value = 5),
                         numericInput("genecountFontSize_xy_axis", label = "Axis labels", value = 18),
                         numericInput("genecountFontSize_legend_title", label = "Legend title", value = 16),
                         numericInput("genecountFontSize_legend_text", label = "Legend labels", value = 15)
                     )),
    
    conditionalPanel("input.navigationTabs == 'multigenecountPlotTab'",
                     div(id = 'multigenecountPlotTab_sidebar',
                         
                         textInput("multi_gene_name", "Enter gene names separated by a comma", value = "HRAS,KRAS,NRAS"),
                         
                         h4("Appearance"),
                         numericInput("multi_genecountGridRows", label = "Grid rows", value = 2),
                         numericInput("multi_genecountGridColumns", label = "Grid columns", value = 3),
                         selectInput("multi_readcountplot_type", label = "Plot type",
                                     choices = list("Boxplot" = 1, "Jitter plot" = 2),
                                     selected = 1),
                         numericInput("multi_genecountPointSize", label = "Jitter point size", value = 3),
                         selectInput("multi_readcountplot_labels", label = "Sample labels",
                                     choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3), 
                                     selected = 1),
                         checkboxInput("multi_genecountSharedYAxis", label = "Label y-axis on first plot per row", value = TRUE),
                         checkboxInput("multi_log10scale", label = "Log10 scale y-axis", value = FALSE),
                         selectInput("multi_genecountShowLegends", label = "Show legends", 
                                     choices = list("Hide legends" = 1, "Show legends on all plots" = 2, "Show one common legend" = 3), 
                                     selected = 3),
                         selectInput("multi_genecountLegendPosition", label = "Legend position", 
                                     choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                     selected = "Top"),
                         
                         h4("Font sizes"),
                         numericInput("multi_genecountFontSize_plot_title", label = "Gene name", value = 15),
                         numericInput("multi_genecountLabelFontSize", label = "Sample labels", value = 5),
                         numericInput("multi_genecountFontSize_xy_axis", label = "Axis labels", value = 12),
                         numericInput("multi_genecountFontSize_legend_text", label = "Legend labels", value = 12)
                     )),
    
    conditionalPanel("input.navigationTabs == 'volcanoPlotTab'",
                     div(id = 'volcanoPlotTab_sidebar',
                         
                         textInput("FCcutoff", "Log2 fold change cutoff", value = 1, width = NULL,
                                   placeholder = NULL),
                         textInput("padjcutoff", "Adjusted p value cutoff", value = "0.05", width = NULL,
                                   placeholder = NULL),
                         
                         h4("Appearance"),
                         numericInput("volcanoPointSize", label = "Point size", value = 3),
                         checkboxInput("volcanoCutoffLines", label = "Show cutoff lines", value = TRUE),
                         
                         h4("Font sizes"),
                         numericInput("volcanoFontSize_plot_title", label = "Plot title", value = 15),
                         numericInput("volcanoFontSize_label", label = "Point labels", value = 5),
                         numericInput("volcanoFontSize_xy_axis", label = "Axis labels", value = 15),
                         numericInput("volcanoFontSize_legend_title", label = "Legend labels", value = 15)
                         
                     ))
    
  ), #end dashboard Sidebar
  
  dashboardBody(
    
    tabsetPanel(
      
      id = "navigationTabs",
      
      #Load data tab
      tabPanel("Load data", id = "loadDataTab", value= "loadDataTab", fluidRow(
        
        #Load read counts file
        fileInput("rawreadsfile", "Select file with read counts",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        checkboxInput("header1", "Header", TRUE),
        radioButtons("sep1", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        #Load coldata
        fileInput("coldatafile", "Select file with sample treatment matrix",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        checkboxInput("header2", "Header", TRUE),
        radioButtons("sep2", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ",")
        
        # #Input for conditions
        # #selectInput("condslist", "Dataset", c("Complete step 1"))
        # uiOutput("cityControls")
        
      )),
      
      #Experiment settings tab
      tabPanel("Experiment settings", id = "expSettingsTab", value= "expSettingsTab", 
               
               selectInput("ref_genome_organism", label = "Reference organism",
                           choices = list("Human" = 1, "Mouse" = 2), 
                           selected = 1),
               uiOutput("control_condslist"),
               uiOutput("treatment1_condslist"),
               uiOutput("FDR_value"),
               uiOutput("min_reads")
      ),
      
      #DE gane table tab
      tabPanel("DESeq2 output", id = "geneTableTab", value= "geneTableTab",  fluidRow(
        DT::dataTableOutput("calc_res_values")
        
      )),
      
      #PCA plot tab
      tabPanel("PCA plot", id = "pcaPlot", value= "pcaPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("PCA_plot", height = "800", width = "800")
        )
      )),
      
      #Multiple gene read count plots
      tabPanel("Read count plots", id = "multigenecountPlotTab", value= "multigenecountPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("multi_genecount_plot1", height = "800", width="800")
        )
      )),
      
      #Volcano plot
      tabPanel("Volcano plot", id = "volcanoPlotTab", value= "volcanoPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("volcanoPlot", height = "800", width="800")
        )
      ))
      
    )
    
  ), #end dashboardBody
  
  skin = "black"
  
) #end dashboardPage