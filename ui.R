# BEAVR: A Browser-based tool for the Exploration And Visualization of RNA-seq data
# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# See Github for more info & ReadMe: https://github.com/developerpiru/BEAVR

app_version = "0.73.1"

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
# ability to customize font sizes and point sizes for all graphs/plots
# ability to plot multiple read count plots at once
# customize legend positions on multiple read count plots
# drag to customize the area of all plots
# option to show y-axis title only on first plot per row
# option to show log10 scale y-axis
# dropped single read plot feature - use multi version with 1x1 grid for a single plot
# ability to pick human or mouse reference genomes to map ENSEMBL IDs
# ability to filter results table based on any/all columns
# ability to turn filtering on/off
# ability to download filtered results table
# volcano plot now shows filtered results if filtering is enabled
# updated UI colours and other aesthetics
# fixed legend positions for multi read count plots
# sample clustering plot (pearson correlation, euclidean, etc)
# count matrix heatmap to show most significant genes with highest variance between condition and treatment groups
# sample clustering heatmap colors
# fixed colors for all heatmaps
# specify distance and clustering type for count matrix heatmap

# bugs"
#### PCA, gene count, volcano plots don't auto-update to new dds dataset after changing treatment condition factor level
#### legend symbols show letter 'a' below symbol on jitter plots

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
installReqs("BiocManager", bioc = FALSE)
installReqs("shiny", bioc = FALSE)
installReqs("shinydashboard", bioc = FALSE)
installReqs("plotly", bioc = FALSE)
installReqs("ggplot2", bioc = FALSE)
installReqs("ggrepel", bioc = FALSE)
installReqs("data.table", bioc = FALSE)
installReqs("DT", bioc = FALSE)
installReqs("DESeq2", bioc = TRUE)
installReqs("vsn", bioc = TRUE)
installReqs('apeglm', bioc = TRUE)
installReqs('org.Hs.eg.db', bioc = TRUE)
installReqs('org.Mm.eg.db', bioc = TRUE)
installReqs('EnhancedVolcano', bioc = TRUE)
installReqs('gridExtra', bioc = FALSE)
installReqs('ggpubr', bioc = FALSE)
installReqs('shinyjqui', bioc = FALSE)
installReqs('scales', bioc = FALSE)
installReqs('RColorBrewer', bioc = FALSE)
installReqs('pheatmap', bioc = FALSE)
installReqs('colourpicker', bioc = FALSE)

#load required libraries
library("BiocManager")
library("shiny")
library("shinydashboard")
library("plotly")
library("ggplot2")
library("ggrepel")
library("data.table")
library("DT")
library("DESeq2")
library("vsn")
library('apeglm')
library('org.Hs.eg.db')
library('org.Mm.eg.db')
library('EnhancedVolcano')
library("gridExtra")
library("ggpubr")
library("shinyjqui")
library("scales")
library("RColorBrewer")
library("pheatmap")
library("colourpicker")

ui <- dashboardPage(
  dashboardHeader(title = paste("BEAVR", app_version, sep = " "), titleWidth = 300),
  
  dashboardSidebar(width = 300,
    
    tags$style(".skin-purple .sidebar 
               a { 
                  color: #444; 
               }"),
    
    tags$style(".skin-purple .sidebar 
               a.sidebarlink:link, a.sidebarlink:visited { 
                                    color: #FFF;
               }"),
    tags$style(".skin-purple .sidebar
                a.sidebarlink:hover {
                                    color: #777;
               }"),
    
    tags$style(".skin-purple .sidebar
                .center {
                        text-align: center;
               }"),
    
    tags$style(".skin-purple .sidebar
                .borderbox {
                        border: 2px solid #666;
                        padding: 5px 5px 5px 5px;
                        margin: 5px;
               }"),
    
    #Conditional Panels to show tab-specific settings in sidebar
    
    conditionalPanel("input.navigationTabs == 'loadDataTab'",
                     div(id = 'loadDataTab_sidebar', 
                         tags$div('class'="center",
                         h4("Welcome to BEAVR!"),
                         HTML('A <b>B</b>rowser-based tool for the <b>E</b>xploration <b>A</b>nd <b>V</b>isualization of <b>R</b>NAseq data<br><br>'),
                         tags$p(
                         tags$a(href="https://github.com/developerpiru/VisualRNAseq",
                                target="_blank",
                                class ="sidebarlink",
                                "Check GitHub for help & info")
                         ))
                     )),
    
    conditionalPanel("input.navigationTabs == 'expSettingsTab'",
                     div(id = 'expSettingsTab_sidebar' 
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'geneTableTab'",
                     div(id = 'geneTableTab_sidebar',
                         
                         tags$div('class'="center", 
                           #h4("Download Table"),
                           downloadButton("downloadDEGeneTable", "Download Table")
                         ),

                         checkboxInput("filterTableEnabled", label = tags$b("Enable filtering"), value = FALSE),
                         
                         h4("Common filtering options"),
                         tags$div('class'="borderbox",
                           #log2FoldChange
                           tags$b("log2FoldChange"),
                           numericInput("log2FC_min", label = "Min", value = 0),  
                           numericInput("log2FC_max", label = "Max", value = 0),
                           #pvalue
                           tags$b("pvalue"),
                           numericInput("pvalue_min", label = "Min", value = 0),
                           numericInput("pvalue_max", label = "Max", value = 0),
                           #padj
                           tags$b("padj"),
                           numericInput("padj_min", label = "Min", value = 0),
                           numericInput("padj_max", label = "Max", value = 0)
                         ),
                         
                         h4("More filtering options"),
                         tags$div('class'="borderbox",
                           #baseMean
                           tags$b("baseMean"),
                           numericInput("baseMean_min", label = "Min", value = 0),
                           numericInput("baseMean_max", label = "Max", value = 0),
                           #lfcSE
                           tags$b("lfcSE"),
                           numericInput("lfcSE_min", label = "Min", value = 0),
                           numericInput("lfcSE_max", label = "Max", value = 0),
                           #stat
                           tags$b("stat"),
                           numericInput("stat_min", label = "Min", value = 0),
                           numericInput("stat_max", label = "Max", value = 0)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'pcaPlotTab'",
                     div(id = 'pcaPlotTab_sidebar', 
                         
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                           selectInput("PCAplot_labels", label = "Sample labels",
                                       choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3),
                                       selected = 1),
                           numericInput("pcaPointSize", label = "Point size", value = 3)
                         ),
                         
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                           numericInput("pcaLabelFontSize", label = "Sample labels", value = 5),
                           numericInput("pcaFontSize_xy_axis", label = "Axis labels", value = 18),
                           numericInput("pcaFontSize_legend_title", label = "Legend title", value = 16),
                           numericInput("pcaFontSize_legend_text", label = "Legend labels", value = 15)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'sampleClusteringPlotTab'",
                     div(id = 'sampleClusteringPlotTab_sidebar', 
                         
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                                  selectInput("sampleClustering_method", label = "Distance method",
                                              choices = list("Euclidean" = "euclidean", "Pearson" = "correlation", "Maximum" = "maximum",
                                                            "Manhattan" = "manhattan", "Canberra" = "canberra", "Binary" = "binary",
                                                            "Minkowski" = "minkowski"),
                                              selected = 1),
                                  selectInput("sampleClustering_cellNums", label = "Distance values",
                                              choices = list("Hide" = FALSE, "Decimal" = "%.2f", "Exponential" = "%.1e"),
                                              selected = 1),
                                  colourInput("sampleClustering_cellNumsColor", "Distance value color", "black"),
                                  numericInput("sampleClustering_treeHeightRows", label = "Row tree height", value = 50),
                                  numericInput("sampleClustering_treeHeightCols", label = "Column tree width", value = 50),
                                  selectInput("sampleClustering_mapColor", label = "Heatmap color",
                                              choices = list("Blues" = "Blues",
                                                             "Blue-Purple" = "BuPu",
                                                             "Green-Blue" = "GnBu",
                                                             "Greens" = "Greens",
                                                             "Greys" = "Greys",
                                                             "Oranges" = "Oranges", 
                                                             "Orange-Red" = "OrRd",
                                                             "Purple-Blue" = "PuBu",
                                                             "Purple-Blue-Green" = "PuBuGn",
                                                             "Purple-Red" = "PuRd",
                                                             "Purples" = "Purples", 
                                                             "Red-Purple" = "RdPu", 
                                                             "Reds" = "Reds", 
                                                             "Yellow-Green" = "YlGn",
                                                             "Yellow-Green-Blue" = "YlGnBu",
                                                             "Yellow-Orange-Brown" = "YlOrBr",
                                                             "Yellow-Orange-Red" = "YlOrRd"
                                              ), selected = 1),
                                  colourInput("sampleClustering_borderColor", "Border color", "#FFFFFF", allowTransparent = TRUE),
                                  checkboxInput("sampleClustering_legend", label = tags$b("Show legend"), value = TRUE)
                                  
                         ),
                         
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                                  numericInput("sampleClustering_fontsize", label = "Labels", value = 15),
                                  numericInput("sampleClustering_fontsize_cellNums", label = "Cell values", value = 15)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'countMatrixHeatmapTab'",
                     div(id = 'countMatrixHeatmapTab_sidebar', 
                         
                         h4("Clustering"),
                         tags$div('class'="borderbox",
                                  selectInput("heatmap_clustMethod", label = "Clustering method",
                                              choices = list("Ward.D" = "ward.D", 
                                                             "Ward.D2" = "ward.D2", 
                                                             "Single" = "single", 
                                                             "Complete" = "complete", 
                                                             "Average" = "average", 
                                                             "McQuitty" = "mcquitty", 
                                                             "Median" = "median", 
                                                             "Centroid" = "centroid"),
                                              selected = 1),
                                  checkboxInput("heatmap_clustRows", label = tags$b("Cluster rows"), value = TRUE),
                                  checkboxInput("heatmap_clustCols", label = tags$b("Cluster columns"), value = FALSE)
                         
                         ),
                         
                         h4("Heatmap colors"),
                         tags$div('class'="borderbox",
                                  colourInput("heatmap_lowColor", "Low color", "#374AB3", allowTransparent = FALSE),
                                  colourInput("heatmap_midColor", "Mid color", "#FFFFFF", allowTransparent = FALSE),
                                  colourInput("heatmap_highColor", "High color", "#E62412", allowTransparent = FALSE),
                                  colourInput("heatmap_borderColor", "Border color", "#FFFFFF", allowTransparent = TRUE)
                         ),
                         
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                                  selectInput("heatmap_distance", label = "Distance method",
                                              choices = list("Euclidean" = "euclidean", "Pearson" = "correlation", "Maximum" = "maximum",
                                                             "Manhattan" = "manhattan", "Canberra" = "canberra", "Binary" = "binary",
                                                             "Minkowski" = "minkowski"),
                                              selected = 1),
                                  selectInput("heatmap_cellNums", label = "Distance values",
                                              choices = list("Hide" = FALSE, "Decimal" = "%.2f", "Exponential" = "%.1e"),
                                              selected = 1),
                                  colourInput("heatmap_cellNumsColor", "Distance value color", "black"),
                                  
                                  numericInput("heatmap_treeHeightRows", label = "Row tree height", value = 50),
                                  numericInput("heatmap_treeHeightCols", label = "Column tree width", value = 50),
                                  selectInput("heatmap_showGeneNames", label = "Gene names",
                                              choices = list("HGNC symbols" = "HGNC", "ENSEMBL IDs" = "ENSEMBL", "Hide" = "Hide")),
                                  checkboxInput("heatmap_legend", label = tags$b("Show legend"), value = TRUE)
                                  
                         ),

                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                                  numericInput("heatmap_fontsize_geneNames", label = "Gene names", value = 10),
                                  numericInput("heatmap_fontsize_sampleNames", label = "Sample names", value = 10),
                                  numericInput("heatmap_fontsize", label = "Annotations", value = 10),
                                  numericInput("heatmap_fontsize_cellNums", label = "Cell values", value = 10)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'genecountPlotTab'",
                     div(id = 'genecountPlotTab_sidebar',

                         textInput("gene_name", "Enter gene name", value = "KRAS"),

                         h4("Appearance"),
                         tags$div('class'="borderbox",
                           selectInput("readcountplot_type", label = "Plot type",
                                       choices = list("Boxplot" = 1, "Jitter plot" = 2),
                                       selected = 1),
                           numericInput("genecountPointSize", label = "Jitter point size", value = 3),
                           selectInput("readcountplot_labels", label = "Sample labels",
                                       choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3),
                                       selected = 1)
                         ),

                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                           numericInput("genecountFontSize_plot_title", label = "Gene name", value = 20),
                           numericInput("genecountLabelFontSize", label = "Sample labels", value = 5),
                           numericInput("genecountFontSize_xy_axis", label = "Axis labels", value = 18),
                           numericInput("genecountFontSize_legend_title", label = "Legend title", value = 16),
                           numericInput("genecountFontSize_legend_text", label = "Legend labels", value = 15)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'multigenecountPlotTab'",
                     div(id = 'multigenecountPlotTab_sidebar',
                         
                         textInput("multi_gene_name", "Enter gene names separated by a comma", value = "HRAS,KRAS,NRAS"),
                         
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                           numericInput("multi_genecountGridRows", label = "Grid rows", value = 1),
                           numericInput("multi_genecountGridColumns", label = "Grid columns", value = 3),
                           selectInput("multi_readcountplot_type", label = "Plot type",
                                       choices = list("Boxplot" = 1, "Jitter plot" = 2),
                                       selected = 1),
                           numericInput("multi_genecountPointSize", label = "Jitter point size", value = 3),
                           selectInput("multi_readcountplot_labels", label = "Sample labels",
                                       choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3), 
                                       selected = 1),
                           checkboxInput("multi_genecountSharedYAxis", label = tags$b("Label y-axis on first plot per row"), value = TRUE),
                           checkboxInput("multi_log10scale", label = tags$b("Log10 scale y-axis"), value = FALSE),
                           selectInput("multi_genecountShowLegends", label = "Show legends", 
                                       choices = list("Hide legends" = 1, "Show legends on all plots" = 2, "Show one common legend" = 3), 
                                       selected = 3),
                           selectInput("multi_genecountLegendPosition", label = "Legend position", 
                                       choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                       selected = "Top")
                         ),
                         
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                           numericInput("multi_genecountFontSize_plot_title", label = "Gene name", value = 15),
                           numericInput("multi_genecountLabelFontSize", label = "Sample labels", value = 5),
                           numericInput("multi_genecountFontSize_xy_axis", label = "Axis labels", value = 12),
                           numericInput("multi_genecountFontSize_legend_text", label = "Legend labels", value = 12)
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'volcanoPlotTab'",
                     div(id = 'volcanoPlotTab_sidebar',
                         
                         h4("Cutoffs"),
                         tags$div('class'="borderbox",
                           textInput("FCcutoff", "Log2 fold change", value = 1, width = NULL,
                                     placeholder = NULL),
                           textInput("padjcutoff", "Adjusted p value", value = "0.05", width = NULL,
                                     placeholder = NULL)
                         ),
                         
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                           numericInput("volcanoPointSize", label = "Point size", value = 3),
                           checkboxInput("volcanoCutoffLines", label = tags$b("Show cutoff lines"), value = TRUE)
                         ),
                         
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                           numericInput("volcanoFontSize_plot_title", label = "Plot title", value = 15),
                           numericInput("volcanoFontSize_label", label = "Point labels", value = 5),
                           numericInput("volcanoFontSize_xy_axis", label = "Axis labels", value = 15),
                           numericInput("volcanoFontSize_legend_title", label = "Legend labels", value = 15)
                         )
                         
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
                             ".csv", ".txt")),

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
                             ".csv", ".txt")),

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
      tabPanel("Settings", id = "expSettingsTab", value= "expSettingsTab", 
               
               selectInput("ref_genome_organism", label = "Reference organism",
                           choices = list("Human" = 1, "Mouse" = 2), 
                           selected = 1),
               uiOutput("control_condslist"),
               uiOutput("treatment1_condslist"),
               uiOutput("FDR_value"),
               uiOutput("min_reads")
      ),
      
      #DE gane table tab
      tabPanel("Gene Table", id = "geneTableTab", value= "geneTableTab",  fluidRow(
        DT::dataTableOutput("calc_res_values")
        
      )),
      
      #PCA plot tab
      tabPanel("PCA", id = "pcaPlot", value= "pcaPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("PCA_plot", height = "500", width = "500")
        )
      )),
      
      #Sample clustering plot tab
      tabPanel("Sample clustering", id = "sampleClusteringPlot", value= "sampleClusteringPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("sampleClustering_plot", height = "500", width = "500")
        )
      )),
      
      #Count matrix heatmap tab
      tabPanel("Heatmap", id = "countMatrixHeatmap", value= "countMatrixHeatmapTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("countMatrix_heatmap", height = "500", width = "500")
        )
      )),
      
      #Multiple gene read count plots
      tabPanel("Read count plots", id = "multigenecountPlotTab", value= "multigenecountPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("multi_genecount_plot1", height = "400", width="800")
        )
      )),
      
      #Volcano plot
      tabPanel("Volcano plot", id = "volcanoPlotTab", value= "volcanoPlotTab", fluidRow(
        jqui_resizable( #jqui resizable canvas
          plotOutput("volcanoPlot", height = "500", width="800")
        )
      ))
      
    )
    
  ), #end dashboardBody
  
  skin = "purple"
  
  
) #end dashboardPage