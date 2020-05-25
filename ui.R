# BEAVR: A Browser-based tool for the Exploration And Visualization of RNA-seq data
# Developed by Pirunthan Perampalam @ https://github.com/developerpiru/
# See Github for documentation & ReadMe: https://github.com/developerpiru/BEAVR

app_version = "1.0.10"

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
# fix variance transformations for small nsubs (small sample sets)
# volcano plot colors
# enter gene names for heatmap
# fixed color selection for boxplot and jitter plots
# calculate statistics in read count plots
# fixed statistics placements on read count plots
# updated count matrix heatmap function to use ComplexHeatmap
# fixed annoations for count matrix heatmap
# ability to customize colors of count matrix heatmap annotations
# full customization of count matrix heatmap now working
# improve dynamic colorbox rendering
# global num_conditions and num_replicates values added
# fixed heatmap annotation customizations and fonts
# div tags for all plot areas showing plot area boundary
# save all plots as png, jpg, svg, tiff or pdf
# select dpi setting for svg and pdf formats
# cleaned up ui
# pathway enrichment analysis using ReactomePA and enrichplot packages (overenrichment analysis and GSEA maps and plots)
# full customization of pathway and gsea plots/maps
# added results tables for pathway enrichment results and GSEA results
# added shiny.port option to use port 3838
# start info bar containing basic steps
# help tab for basic help/tips info
# fixed heatmap name bug

# bugs"
#### PCA, gene count, volcano plots don't auto-update to new dds dataset after changing treatment condition factor level
#### legend symbols show letter 'a' below symbol on jitter plots

#set port to 3838
options(shiny.port = 3838)

## load required libraries
library("shiny")
library("shinydashboard")
library("shinyWidgets")

#cran packages
library("BiocManager")
library("colourpicker")
library("data.table")
#library("devtools") # to install from github
library("DT")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ggraph")
library("gridExtra")
library("pheatmap")
library("RColorBrewer")
library("scales")
library("shiny")
library("shinydashboard")
library("shinyjqui")
library("shinyWidgets")
library("shinycssloaders")
library("circlize")
library(shinyalert)

# #Bioconductor packages
library("DESeq2")
library("vsn")
library("apeglm")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("ReactomePA") #from bioconductor
library("enrichplot") #from bioconductor

# #GitHub packages
library("EnhancedVolcano")
library("ComplexHeatmap")

spinner_type = 1
spinner_col = "#0C71CF"

ui <- dashboardPage(
  dashboardHeader(title = "BEAVR", titleWidth = 300),
  
  dashboardSidebar(width = 300,
    
    tags$style("

               .skin-blue .sidebar 
               a.sidebarlink:link, a.sidebarlink:visited { 
                        color: #FFFFFF;
               }
               
               .skin-blue .sidebar
               a.sidebarlink:hover {
                        color: #FFFFFF;
               }               

               .skin-blue .sidebar
               .treeview {
                        color: #FFFFFF;
               } 
    
               .skin-blue .sidebar
               .center {
                        text-align: center;
               }

               .skin-blue .sidebar               
               a.shiny-download-link:link, a.shiny-download-link:visited {
                        color: #000000;
               }
               
               .skin-blue .sidebar               
               a.shiny-download-link:hover {
                        color: #000000;
               }

               .skin-blue .sidebar
               .borderbox {
                        padding: 2px 0px 0px 0px;
               }

               .skin-blue .sidebar
               .infobox {
                        padding: 2px 2px 10px 10px;
                        margin: 5px 5px 5px 5px;
                        display:block;
                        clear:both;
                        white-space: normal;
               }

               .content-wrapper, .right-side {
                        background-color: #FFFFFF;
                        overflow-x: auto;
               }

               "),
    
    
    
    sidebarMenu(
      
    #Conditional Panels to show tab-specific settings in sidebar
    conditionalPanel("input.navigationTabs == 'loadDataTab'",

                     div(id = 'loadDataTab_sidebar',
                         tags$div('class'="infobox",
                                  htmlOutput("startinfo")
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'expSettingsTab'",
                     div(id = 'expSettingsTab_sidebar' 
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'geneTableTab'",
                     div(id = 'geneTableTab_sidebar',
                         
                         tags$div('class'="center", 
                           tags$br(),
                           downloadButton("downloadDEGeneTable", "Download Table", class = "btn")
                         ),

                         materialSwitch("filterTableEnabled", label = tags$b("Filter table"), value = FALSE, status = "primary"),
                         
                         menuItem(
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
                             
                           )
                         ),
                         menuItem(
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
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'pcaPlotTab'",
                     div(id = 'pcaPlotTab_sidebar', 
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                             selectInput("PCAplot_labels", label = "Sample labels",
                                         choices = list("No labels" = 1, "Sample names" = 2, "Replicate names" = 3),
                                         selected = 1),
                             numericInput("pcaPointSize", label = "Point size", value = 3)
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    tags$div(id="pcaColorbox")
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                             numericInput("pcaLabelFontSize", label = "Sample labels", value = 5),
                             numericInput("pcaFontSize_xy_axis", label = "Axis labels", value = 18),
                             numericInput("pcaFontSize_legend_title", label = "Legend title", value = 16),
                             numericInput("pcaFontSize_legend_text", label = "Legend labels", value = 15)
                           )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'sampleClusteringPlotTab'",
                     div(id = 'sampleClusteringPlotTab_sidebar', 
                         
                         menuItem(
                         h4("Appearance"),
                         tags$div('class'="borderbox",
                                  selectInput("sampleClustering_method", label = "Distance method",
                                              choices = list("Euclidean" = "euclidean", "Pearson" = "pearson", "Spearman" = "spearman", 
                                                             "Kendall" = "kendall", "Maximum" = "maximum", "Manhattan" = "manhattan", 
                                                             "Canberra" = "canberra", "Binary" = "binary", "Minkowski" = "minkowski"),
                                              selected = "euclidean"),
                                  selectInput("sampleClustering_cellNums", label = "Distance values",
                                              choices = list("Hide" = FALSE, "Decimal" = "%.2f", "Exponential" = "%.1e"),
                                              selected = FALSE),
                                  
                                  numericInput("sampleClusterin_row_dend_width", label = "Row dendrogram width", value = 2),
                                  numericInput("sampleClusterin_col_dend_height", label = "Column dendrogram height", value = 2),
                                  selectInput("sampleClusterin_row_dend_position", label = "Row dendrogram position",
                                              choices = list("Left" = "left", "Right" = "right"),
                                              selected = "left"),
                                  selectInput("sampleClusterin_col_dend_position", label = "Column dendrogram position",
                                              choices = list("Top" = "top", "Bottom" = "bottom"),
                                              selected = "top"),
                                  selectInput("sampleClusterin_rowlabel_position", label = "Row label position",
                                              choices = list("Left" = "left", "Right" = "right"),
                                              selected = "right"),
                                  selectInput("sampleClusterin_collabel_position", label = "Column label position",
                                              choices = list("Top" = "top", "Bottom" = "bottom"),
                                              selected = "bottom"),
                                  numericInput("sampleClusterin_row_rotation", label = "Row label rotation (degrees)", value = 0),
                                  numericInput("sampleClusterin_col_rotation", label = "Column label rotation (degrees)", value = 45)
                                  
                         )
                         ),
                         
                         menuItem(
                           h4("Heatmap colors"),
                           tags$div('class'="borderbox",
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
                                              ), selected = "Blues"),
                                  colourInput("sampleClustering_borderColor", "Border color", "#FFFFFF00", allowTransparent = TRUE)
                           )
                         ),
                         
                         menuItem(
                           h4("Label colors"),
                           tags$div('class'="borderbox",
                                    colourInput("sampleClustering_rowlabelColor", "Row labels", "black"),
                                    colourInput("sampleClustering_collabelColor", "Column labels", "black"),
                                    colourInput("sampleClustering_cellNumsColor", "Distance values", "black"),
                                    colourInput("sampleClustering_legendColor", "Legend labels", "black")
                           )
                         ),
                         
                         menuItem(
                           h4("Legends"),
                           tags$div('class'="borderbox",
                                    selectInput("sampleClustering_main_legend", label = "Legend position",
                                                choices = list("Left" = "left", 
                                                               "Right" = "right", "Top" = "top", "Bottom" = "bottom"),
                                                selected = "right"),
                                    selectInput("sampleClustering_main_legend_dir", label = "Legend direction",
                                                choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                                selected = "horizontal"),
                                    numericInput("sampleClustering_main_legend_size", label = "Legend size", value = 5)
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                                    numericInput("sampleClustering_fontsize_rowNames", label = "Row names", value = 10),
                                    numericInput("sampleClustering_fontsize_colNames", label = "Columns names", value = 10),
                                    numericInput("sampleClustering_fontsize_cellNums", label = "Cell values", value = 10),
                                    numericInput("sampleClustering_fontsize_legends", label = "Legend labels", value = 10)
                           )
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'multigenecountPlotTab'",
                     div(id = 'multigenecountPlotTab_sidebar',
                         
                         menuItem(
                           h4("Genes"),
                           tags$div('class'="borderbox",
                             textAreaInput("multi_gene_name", label = "Enter gene names separated by comma", value = "HRAS,NRAS,KRAS", width = NULL,
                                           height = 100, cols = NULL, rows = NULL, placeholder = NULL,
                                           resize = NULL)
                           )
                         ),
                         
                         menuItem(
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
                             materialSwitch("multi_genecountSharedYAxis", label = tags$b("Label y-axis on first plot per row"), value = TRUE, status = "primary"),
                             materialSwitch("multi_log10scale", label = tags$b("Log10 scale y-axis"), value = FALSE, status = "primary"),
                             materialSwitch("pcaRotateText", label = tags$b("Rotate x-axis labels"), value = TRUE, status = "primary"),
                             selectInput("multi_genecountShowLegends", label = "Show legends", 
                                         choices = list("Hide legends" = 1, "Show legends on all plots" = 2, "Show one common legend" = 3), 
                                         selected = 3),
                             selectInput("multi_genecountLegendPosition", label = "Legend position", 
                                         choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                         selected = "Top")
                           )
                         ),
                         
                         menuItem(
                           h4("Statistics"),
                           tags$div('class'="borderbox",
                             materialSwitch("multi_genecountStats", label = tags$b("Show statistics"), value = TRUE, status = "primary"),
                             numericInput("multi_genecountStatsYcord", label = tags$b("Adjust statistics placement"), value = "1")
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    tags$div(id="multi_genecountColorbox")
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                             numericInput("multi_genecountFontSize_plot_title", label = "Gene name", value = 15),
                             numericInput("multi_genecountLabelFontSize", label = "Sample labels", value = 5),
                             numericInput("multi_genecountFontSize_xy_axis", label = "Axis labels", value = 12),
                             numericInput("multi_genecountFontSize_legend_text", label = "Legend labels", value = 12),
                             numericInput("multi_genecountFontSize_stats_text", label = "Statistics", value = 6)
                           )
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'countMatrixHeatmapTab'",
                     div(id = 'countMatrixHeatmapTab_sidebar', 
                         
                         menuItem(
                         h4("Genes"),
                         tags$div('class'="borderbox",
                                  textAreaInput("heatmap_GeneNames", label = "Enter gene names separated by comma", value = "", 
                                                width = NULL, height = 100, cols = NULL, rows = NULL, placeholder = NULL,
                                                resize = NULL),
                                  materialSwitch("heatmap_pickTopGenes", label = tags$b("Show top genes instead"), value = TRUE, status = "primary"),
                                  numericInput("heatmap_numGenes", label = "Number of top genes to show", value = 5),
                                  selectInput("heatmap_showGeneNames", label = "Gene names",
                                              choices = list("HGNC symbols" = "HGNC", "ENSEMBL IDs" = "ENSEMBL"))

                         )
                         ),
                         
                         menuItem(
                         h4("Clustering"),
                         tags$div('class'="borderbox",
                                  selectInput("heatmap_varlogmethod", label = "Variance stabilization method",
                                              choices = list("Regularized log transformation" = "rlog",
                                                             "Variance stabilization" = "vst"),
                                              selected = "vst"),
                                  selectInput("heatmap_clustMethod", label = "Clustering method",
                                              choices = list("Ward D1 (Ward 1963)" = "ward.D", 
                                                             "Ward D2 (Complete Ward's criteriom)" = "ward.D2", 
                                                             "Single" = "single", 
                                                             "Complete" = "complete", 
                                                             "Average" = "average", 
                                                             "McQuitty" = "mcquitty", 
                                                             "Median" = "median", 
                                                             "Centroid" = "centroid"),
                                              selected = "ward.D"),
                                  selectInput("heatmap_distance", label = "Distance method",
                                              choices = list("Euclidean" = "euclidean", "Pearson" = "pearson", "Spearman" = "spearman", 
                                                             "Kendall" = "kendall", "Maximum" = "maximum", "Manhattan" = "manhattan", 
                                                             "Canberra" = "canberra", "Binary" = "binary", "Minkowski" = "minkowski"),
                                              selected = "euclidean"),
                                  materialSwitch("heatmap_clustRows", label = tags$b("Cluster rows"), value = TRUE, status = "primary"),
                                  materialSwitch("heatmap_clustCols", label = tags$b("Cluster columns"), value = FALSE, status = "primary")
                         
                         )          
                         ),
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                                    sliderInput("heatmap_scale_range", label = "Gene expression scale", 
                                                min = -6, max = 6, value = c(-2,2), step = 0.5),
                                    selectInput("heatmap_annotations", label = "Annotations",
                                                choices = list("None" = "none", "Both" = "both", 
                                                               "Replicate" = "replicate", "Treatment" = "treatment"),
                                                selected = "both"),
                                    selectInput("heatmap_cellNums", label = "Distance values",
                                                choices = list("Hide" = FALSE, "Decimal" = "%.2f", "Exponential" = "%.1e"),
                                                selected = FALSE),
                                    
                                    numericInput("heatmap_row_dend_width", label = "Row dendrogram width", value = 2),
                                    numericInput("heatmap_col_dend_height", label = "Column dendrogram height", value = 2),
                                    selectInput("heatmap_row_dend_position", label = "Row dendrogram position",
                                                choices = list("Left" = "left", "Right" = "right"),
                                                selected = "left"),
                                    selectInput("heatmap_col_dend_position", label = "Column dendrogram position",
                                                choices = list("Top" = "top", "Bottom" = "bottom"),
                                                selected = "top"),
                                    selectInput("heatmap_genelabel_position", label = "Gene label position",
                                                choices = list("Left" = "left", "Right" = "right"),
                                                selected = "right"),
                                    selectInput("heatmap_samplelabel_position", label = "Sample label position",
                                                choices = list("Top" = "top", "Bottom" = "bottom"),
                                                selected = "bottom"),
                                    numericInput("heatmap_row_rotation", label = "Gene name rotation (degrees)", value = 0),
                                    numericInput("heatmap_col_rotation", label = "Sample name rotation (degrees)", value = 45)
                           )
                         ),
                         
                         menuItem(
                         h4("Heatmap colors"),
                         tags$div('class'="borderbox",
                                  colourInput("heatmap_lowColor", "Low color", "#374AB3", allowTransparent = FALSE),
                                  colourInput("heatmap_midColor", "Mid color", "#FFFFFF", allowTransparent = FALSE),
                                  colourInput("heatmap_highColor", "High color", "#E62412", allowTransparent = FALSE),
                                  colourInput("heatmap_borderColor", "Border color", "#FFFFFF00", allowTransparent = TRUE)
                         )
                         ),
                         
                         menuItem(
                         h4("Replicate colors"),
                         tags$div('class'="borderbox",
                                  tags$div(id="heatmap_replicateColorbox")
                         )
                         ),
                         
                         menuItem(
                         h4("Condition/Treatment colors"),
                         tags$div('class'="borderbox",
                                  tags$div(id="heatmap_conditionColorbox")
                         )
                         ),
                         
                         menuItem(
                           h4("Label colors"),
                           tags$div('class'="borderbox",
                                    colourInput("heatmap_rowlabelColor", "Gene labels", "black"),
                                    colourInput("heatmap_collabelColor", "Sample labels", "black"),
                                    colourInput("heatmap_cellNumsColor", "Distance values", "black"),
                                    colourInput("heatmap_legendColor", "Legend labels", "black")
                           )
                         ),
                         
                         menuItem(
                         h4("Legends"),
                         tags$div('class'="borderbox",
                                  selectInput("heatmap_main_legend", label = "Main legend position",
                                              choices = list("Left" = "left", 
                                                             "Right" = "right", "Top" = "top", "Bottom" = "bottom"),
                                              selected = "right"),
                                  selectInput("heatmap_anno_legend", label = "Annotations legend position",
                                              choices = list("Left" = "left", 
                                                             "Right" = "right", "Top" = "top", "Bottom" = "bottom"),
                                              selected = "bottom"),
                                  selectInput("heatmap_main_legend_dir", label = "Main legend direction",
                                              choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                              selected = "vertical"),
                                  selectInput("heatmap_anno_legend_dir", label = "Annotations legend direction",
                                              choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                              selected = "vertical"),
                                  numericInput("heatmap_main_legend_size", label = "Main legend size", value = 5)
                         )
                         ),
                         
                         menuItem(
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                                  numericInput("heatmap_fontsize_geneNames", label = "Gene names", value = 10),
                                  numericInput("heatmap_fontsize_sampleNames", label = "Sample names", value = 10),
                                  numericInput("heatmap_fontsize_annotations", label = "Annotations", value = 10),
                                  numericInput("heatmap_fontsize_cellNums", label = "Cell values", value = 10),
                                  numericInput("heatmap_fontsize_legends", label = "Legend labels", value = 10)
                         )
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'volcanoPlotTab'",
                     div(id = 'volcanoPlotTab_sidebar',
                         
                         menuItem(
                         h4("Cutoffs"),
                         tags$div('class'="borderbox",
                           textInput("volcanoFCcutoff", "Log2 fold change", value = 1, width = NULL,
                                     placeholder = NULL),
                           selectInput("volcanopCutoffType", label = "p value cutoff type", 
                                       choices = list("Unadjusted p value" = "pvalue", "Adjusted p value" = "padj"), 
                                       selected = "unadj"),
                           textInput("volcanopCutoff", "p value or padj value", value = "0.05", width = NULL,
                                    placeholder = NULL)
                         )
                         ),
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                                    numericInput("volcanoPointSize", label = "Point size", value = 3),
                                    materialSwitch("volcanoCutoffLines", label = tags$b("Show cutoff lines"), value = TRUE, status = "primary"),
                                    selectInput("volcanoLegendPosition", label = "Legend position", 
                                                choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                                selected = "bottom")
                           )
                         ),
                         
                         menuItem(
                         h4("Colors"),
                         tags$div('class'="borderbox",
                                  colourInput("volcano_NSColor", "Not significant color", "#5B5B5C", allowTransparent = FALSE),
                                  colourInput("volcano_LFCColor", "LFC only color", "#0FBD32", allowTransparent = FALSE),
                                  colourInput("volcano_pvalColor", "p value only color", "#126FE8", allowTransparent = FALSE),
                                  colourInput("volcano_pvalLFCColor", "LFC and p value color", "#FF0000", allowTransparent = TRUE)
                         )
                         ),
                         
                         menuItem(
                         h4("Font sizes"),
                         tags$div('class'="borderbox",
                           numericInput("volcanoFontSize_plot_title", label = "Plot title", value = 15),
                           numericInput("volcanoFontSize_label", label = "Point labels", value = 5),
                           numericInput("volcanoFontSize_xy_axis", label = "Axis labels", value = 15),
                           numericInput("volcanoFontSize_legend_title", label = "Legend labels", value = 15)
                         )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'enrichmentPlotTab'",
                     div(id = 'enrichmentPlotTab_sidebar',
                         
                         menuItem(
                           h4("Enrichment"),
                           tags$div('class'="borderbox",
                                    numericInput("enrPvalcutoff", "p value cutoff for enrichment", value = "0.05"),
                                    numericInput("enrNumCategories", label = "Max number of categories", value = 10)
                           )
                         ),
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                                    selectInput("enrPlotType", label = "Plot type", 
                                                choices = list("Bar plot" = "bar", "Dot plot" = "dot"), 
                                                selected = "bar")
                           )
                         ),
                         
                         menuItem(
                           h4("Legend"),
                           tags$div('class'="borderbox",
                                    selectInput("enrLegendPosition", label = "Legend position", 
                                                choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                                selected = "right"),
                                    selectInput("enrLegendDirection", label = "Legend direction", 
                                                choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"), 
                                                selected = "vertical")
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    colourInput("enrMostColor", "Most significant color", "#374AB3", allowTransparent = FALSE),
                                    colourInput("enrLeastColor", "Least significant color", "#E62412", allowTransparent = FALSE)
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                                    numericInput("enrFontSize_xy_axis", label = "Axis labels", value = 15),
                                    numericInput("enrFontSize_legend", label = "Legend labels", value = 15)
                           )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'enrichmentMapTab'",
                     div(id = 'enrichmentMapTab_sidebar',
                         
                         menuItem(
                           h4("Enrichment"),
                           tags$div('class'="borderbox",
                                    numericInput("enrMapPvalcutoff", "p value cutoff for enrichment", value = "0.05")
                           )
                         ),
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                                    selectInput("enrMapLegendPosition", label = "Legend position", 
                                                choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                                selected = "right"),
                                    selectInput("enrMapLegendDirection", label = "Legend direction", 
                                                choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"), 
                                                selected = "vertical")
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    colourInput("enrMapMostColor", "Most significant color", "#374AB3", allowTransparent = FALSE),
                                    colourInput("enrMapLeastColor", "Least significant color", "#E62412", allowTransparent = FALSE),
                                    colourInput("enrMaplineColor", "Edge/Line color", "#C9C3C3", allowTransparent = FALSE),
                                    sliderInput("enrMaptransparency", label = "Edge/Line transparency",
                                                min = 0, max = 100, value = 50, step = 1, post = "%")
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                                    numericInput("enrMapFontSize_labels", label = "Labels", value = 5),
                                    numericInput("enrMapFontSize_legend", label = "Legend labels", value = 15)
                           )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'gseaMapTab'",
                     div(id = 'gseaMapTab_sidebar',
                         
                         menuItem(
                           h4("Enrichment"),
                           tags$div('class'="borderbox",
                                    numericInput("gseaMapPvalcutoff", "p value cutoff for enrichment", value = "0.05")
                           )
                         ),
                         
                         menuItem(
                           h4("Legend"),
                           tags$div('class'="borderbox",
                                    selectInput("gseaMapLegendPosition", label = "Legend position", 
                                                choices = list("Top" = "top", "Bottom" = "bottom", "Left" = "left", "Right" = "right"), 
                                                selected = "right"),
                                    selectInput("gseaMapLegendDirection", label = "Legend direction", 
                                                choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"), 
                                                selected = "vertical")
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    colourInput("gseaMapMostColor", "Most significant color", "#374AB3", allowTransparent = FALSE),
                                    colourInput("gseaMapLeastColor", "Least significant color", "#E62412", allowTransparent = FALSE),
                                    colourInput("gseaMaplineColor", "Edge/Line color", "#C9C3C3", allowTransparent = FALSE),
                                    sliderInput("gseaMaptransparency", label = "Edge/Line transparency",
                                                min = 0, max = 100, value = 50, step = 1, post = "%")
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                                    numericInput("gseaMapFontSize_labels", label = "Axis labels", value = 5),
                                    numericInput("gseaMapFontSize_legend", label = "Legend labels", value = 15)
                           )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'gseaPlotTab'",
                     div(id = 'gseaPlotTab_sidebar',
                         
                         menuItem(
                           h4("Enrichment"),
                           tags$div('class'="borderbox",
                                    uiOutput("gseaPlotPathways")
                           )
                         ),
                         
                         menuItem(
                           h4("Appearance"),
                           tags$div('class'="borderbox",
                                    selectInput("gseaPlotShowPvalue", label = "p values",
                                                choices = list("Show" = TRUE, "Hide" = FALSE),
                                                selected = FALSE)
                           )
                         ),
                         
                         menuItem(
                           h4("Colors"),
                           tags$div('class'="borderbox",
                                    colourInput("gseaPlotLineColor", "Enrichment score line color", "#1EFF00", allowTransparent = FALSE)
                           )
                         ),
                         
                         menuItem(
                           h4("Font sizes"),
                           tags$div('class'="borderbox",
                                    numericInput("gseaPlotFontSize", label = "Base font size", value = 15)
                           )
                         )
                         
                     )),
    
    conditionalPanel("input.navigationTabs == 'enrichmentTableTab'",
                     div(id = 'enrichmentTableTab_sidebar',
                         
                         tags$div('class'="center", 
                                  tags$br(),
                                  downloadButton("downloadenrichmentTable", "Download Table", class = "btn")
                         ),
                         
                         materialSwitch("enrichmentTablefilterTableEnabled", label = tags$b("Filter table"), value = FALSE, status = "primary"),
                         
                         menuItem(
                           h4("Filter by statistics"),
                           tags$div('class'="borderbox",
                                    tags$b("Count"),
                                    numericInput("enr_count_min", label = "Min", value = 0),
                                    numericInput("enr_count_max", label = "Max", value = 0),
                                    
                                    #pvalue
                                    tags$b("pvalue"),
                                    numericInput("enr_pvalue_min", label = "Min", value = 0),
                                    numericInput("enr_pvalue_max", label = "Max", value = 0),
                                    
                                    #padj
                                    tags$b("padj"),
                                    numericInput("enr_padj_min", label = "Min", value = 0),
                                    numericInput("enr_padj_max", label = "Max", value = 0),
                                    
                                    #qvalue
                                    tags$b("qvalue"),
                                    numericInput("enr_qvalue_min", label = "Min", value = 0),
                                    numericInput("enr_qvalue_max", label = "Max", value = 0)
                                    
                           )
                         )
                     )),
    
    conditionalPanel("input.navigationTabs == 'gseaTableTab'",
                     div(id = 'gseaTableTab_sidebar',
                         
                         tags$div('class'="center", 
                                  tags$br(),
                                  downloadButton("downloadgseaTable", "Download Table", class = "btn")
                         ),
                         
                         materialSwitch("gseaTablefilterTableEnabled", label = tags$b("Filter table"), value = FALSE, status = "primary"),
                         
                         menuItem(
                           h4("Filter by enrichment scores"),
                           tags$div('class'="borderbox",
                                    #enrichment Score
                                    tags$b("Enrichment score (ES)"),
                                    numericInput("gsea_enrichmentScore_min", label = "Min", value = 0),  
                                    numericInput("gsea_enrichmentScore_max", label = "Max", value = 0),
                                    
                                    #NES
                                    tags$b("Normalized enrichment score (NES)"),
                                    numericInput("gsea_nes_min", label = "Min", value = 0),
                                    numericInput("gsea_nes_max", label = "Max", value = 0)
                                    
                           )
                         ),
                         
                         menuItem(
                           h4("Filter by statistics"),
                           tags$div('class'="borderbox",
                                    #pvalue
                                    tags$b("pvalue"),
                                    numericInput("gsea_pvalue_min", label = "Min", value = 0),
                                    numericInput("gsea_pvalue_max", label = "Max", value = 0),
                                    
                                    #padj
                                    tags$b("padj"),
                                    numericInput("gsea_padj_min", label = "Min", value = 0),
                                    numericInput("gsea_padj_max", label = "Max", value = 0),
                                    
                                    #qvalue
                                    tags$b("qvalue"),
                                    numericInput("gsea_qvalue_min", label = "Min", value = 0),
                                    numericInput("gsea_qvalue_max", label = "Max", value = 0),
                                    
                                    tags$b("Rank"),
                                    numericInput("gsea_rank_min", label = "Min", value = 0),
                                    numericInput("gsea_rank_max", label = "Max", value = 0)
                                    
                           )
                         )
                     ))
    
    )
  ), #end dashboard Sidebar
  
  dashboardBody(
    
    tabsetPanel(
      
      id = "navigationTabs",
      
      #Load data tab
      tabPanel("Load Data", id = "loadDataTab", value= "loadDataTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        #Load read counts file
        fileInput("rawreadsfile", "Select read count table file",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),

        radioButtons("sep1", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        #Load coldata
        fileInput("coldatafile", "Select sample treatment matrix file",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),

        radioButtons("sep2", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        HTML('</p></div>')
        
      )),
      
      #Experiment settings tab
      tabPanel("Settings", id = "expSettingsTab", value= "expSettingsTab", 
        
       #layout HTML tags
       HTML('<div style="padding:10px 10px 10px 10px;"><p>'),

       selectInput("ref_genome_organism", label = "Reference organism",
                   choices = list("Human" = 1, "Mouse" = 2),
                   selected = 1),
       uiOutput("control_condslist"),
       uiOutput("treatment1_condslist"),
       uiOutput("FDR_value"),
       selectInput("shrinkage_method", label = "Shrinkage estimator",
                   choices = list("Approximate posterior estimation (apeglm)" = 1, "Normal (adaptive normal distribution)" = 2),
                   selected = 1),
       uiOutput("min_reads"),

       HTML('</p></div>')
       
      ),
      
      #DE gane table tab
      tabPanel("Gene Table", id = "geneTableTab", value= "geneTableTab",  fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        DT::dataTableOutput("calc_res_values"),
        
        HTML('</p></div>')
        
      )),
      
      #PCA plot tab
      tabPanel("PCA", id = "pcaPlot", value= "pcaPlotTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("PCAplot_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("savePCApng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("savePCAjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("savePCAsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("savePCApdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("savePCAtiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for PCA plot
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("PCA_plot", height = "500", width = "500"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')

      )),
      
      #Sample clustering plot tab
      tabPanel("Sample Clustering", id = "sampleClusteringPlot", value= "sampleClusteringPlotTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("sampleClustering_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveClusteringpng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveClusteringjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveClusteringsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveClusteringpdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveClusteringtiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Sample Clustering plot
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("sampleClustering_plot", height = "500", width = "500"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #Multiple gene read count plots
      tabPanel("Read Count Plots", id = "multigenecountPlotTab", value= "multigenecountPlotTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("ReadCount_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveReadCountpng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveReadCountjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveReadCountsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveReadCountpdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveReadCounttiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Read count plots
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("multi_genecount_plot1", height = "400", width="800"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #Count matrix heatmap tab
      tabPanel("Heatmap", id = "countMatrixHeatmap", value= "countMatrixHeatmapTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("Heatmap_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveHeatmappng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveHeatmapjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveHeatmapsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveHeatmappdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveHeatmaptiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Count Matrix Heatmap
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("countMatrix_heatmap", height = "500", width = "500"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')

      )),
      
      #Volcano plot
      tabPanel("Volcano Plot", id = "volcanoPlotTab", value= "volcanoPlotTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),

        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("Volcanoplot_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveVolcanopng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveVolcanojpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveVolcanosvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveVolcanopdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveVolcanotiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),

        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Count Matrix Heatmap
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("volcanoPlot", height = "500", width="800"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #Pathway Enrichment barplot and dotplot
      tabPanel("Pathway Enrichment Plot", id = "enrichmentPlotTab", value= "enrichmentPlotTab", fluidRow(
        
        useShinyalert(),
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("enrichmentplot_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveEnrichmentplotpng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentplotjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentplotsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentplotpdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentplottiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Enrichment plot
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("enrichmentPlot", height = "500", width="800"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #Pathway Enrichment Map
      tabPanel("Pathway Enrichment Map", id = "enrichmentMapTab", value= "enrichmentMapTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("enrichmentmap_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveEnrichmentmappng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentmapjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentmapsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentmappdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveEnrichmentmaptiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for Enrichment map
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("enrichmentMap", height = "500", width="800"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #Pathway Enrichment results table tab
      tabPanel("Pathway Enrichment Table", id = "enrichmentTableTab", value= "enrichmentTableTab",  fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        DT::dataTableOutput("show_enrichmentTable"),
        
        HTML('</p></div>')
        
      )),
      
      #GSEA Map
      tabPanel("GSEA Map", id = "gseaMapTab", value= "gseaMapTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("gseamap_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveGSEAmappng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAmapjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAmapsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAmappdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAmaptiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for GSEA map
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("gseaMap", height = "500", width="800"),
            HTML('</p></div>')
          )
        ),
        
        HTML('</div>')
        
      )),
      
      #GSEA Plot
      tabPanel("GSEA Plot", id = "gseaPlotTab", value= "gseaPlotTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        # download buttons
        dropdown(label = "Save Plot",
                 selectInput("gseaplot_dpi", label = "Output dpi", width = "100",
                             choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                             selected = 48),
                 downloadButton("saveGSEAplotpng", "PNG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAplotjpg", "JPG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAplotsvg", "SVG"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAplotpdf", "PDF"),
                 HTML('<br><br>'),
                 downloadButton("saveGSEAplottiff", "TIFF"),
                 size = "default",
                 icon = icon("download", class = ""), 
                 up = FALSE
        ),
        HTML('</p>'),
        
        jqui_resizable(#jqui resizable canvass
          tagList(
            # div box for GSEA plot
            HTML('<div style="border:1px solid black;padding:5px 5px 0px 0px;width:500;background-color: #FFFFFF;float:left;display:block;"><p>'),
            plotOutput("gseaPlot", height = "500", width="800"),
            HTML('</p></div>')
          )
        ),
        
        #output box to show GSEA pvalues
        HTML('<div style="background-color:#FFFFFF;clear:both;display:block;font:12px bold Arial, Helvetica, sans-serif;"><p>'),
        textOutput("gseaPlotpvals"),
        HTML('</p></div>'),
        
        HTML('</div>')
        
      )),
      
      #GSEA results table tab
      tabPanel("GSEA Table", id = "gseaTableTab", value= "gseaTableTab",  fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
        
        DT::dataTableOutput("show_gseaTable"),
        
        HTML('</p></div>')
        
      )),
      
      #Help/About tab
      tabPanel("Help", id = "helpTab", value= "helpTab", fluidRow(
        
        #layout HTML tags
        HTML('<div style="padding:10px 10px 10px 10px;"><p>'),

        #output box to show GSEA pvalues
        HTML('<div style="background-color:#FFFFFF;clear:both;display:block;font:12px bold Arial, Helvetica, sans-serif;"><p>'),
        htmlOutput("helpinfo"),
        HTML('</p></div>'),
        
        HTML('</p></div>')
        
      ))
      
    )
    
  ), #end dashboardBody
  
  skin = "blue"
  
  
) #end dashboardPage