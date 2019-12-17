# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# version: 0.66

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
# added font size customization to pca, gene count, and volcano plots + point size customization to volcano plot

# bugs"
# PCA, gene count, volcano plots don't auto-update to new dds after changing treatment condition factor level

ui <- dashboardPage(
  dashboardHeader(title = "VisualRNAseq v0.64"),
  
  dashboardSidebar(
    
    #Conditional Panels to show tab-specific settings in sidebar
    
    conditionalPanel("input.navigationTabs == 'pcaPlotTab'",
                    div(id = 'pcaPlotTab_sidebar', 
                        h4("PCA plot"),
                        h5("Font sizes"),
                        numericInput("pcaFontSize_x_axis", label = "x-axis labels", value = 15),
                        numericInput("pcaFontSize_y_axis", label = "y-axis labels", value = 15),
                        numericInput("pcaFontSize_x_title", label = "x-axis title", value = 15),
                        numericInput("pcaFontSize_y_title", label = "y-axis title", value = 15),
                        numericInput("pcaFontSize_legend_title", label = "Legend title", value = 15),
                        numericInput("pcaFontSize_legend_text", label = "Legend text", value = 15)
                    )),
    
    conditionalPanel("input.navigationTabs == 'genecountPlotTab'",
                    div(id = 'genecountPlotTab_sidebar',
                        h4("Gene count plot"),
                        h5("Font sizes"),
                        numericInput("genecountFontSize_x_axis", label = "x-axis labels", value = 15),
                        numericInput("genecountFontSize_y_axis", label = "y-axis labels", value = 15),
                        numericInput("genecountFontSize_x_title", label = "x-axis title", value = 15),
                        numericInput("genecountFontSize_y_title", label = "y-axis title", value = 15),
                        numericInput("genecountFontSize_legend_title", label = "Legend title", value = 15),
                        numericInput("genecountFontSize_legend_text", label = "Legend text", value = 15)
                    )),
    conditionalPanel("input.navigationTabs == 'volcanoPlotTab'",
                    div(id = 'volcanoPlotTab_sidebar',
                        h4("Volcano plot"),
                        h5("Font sizes"),
                        numericInput("volcanoFontSize", label = "Labels", value = 8),
                        h5("Appearance"),
                        numericInput("volcanoPointSize", label = "Point size", value = 3)
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
        fileInput("coldatafile", "Select file with sample treatment/condition/ info",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),

        checkboxInput("header2", "Header", TRUE),
        radioButtons("sep2", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        #Input for conditions
        #selectInput("condslist", "Dataset", c("Complete step 1"))
        uiOutput("cityControls")
        
      )),
      
      #Experiment settings tab
      tabPanel("Select experiment settings", id = "expSettingsTab", value= "expSettingsTab", 
               
               uiOutput("control_condslist"),
               uiOutput("treatment1_condslist"),
               uiOutput("FDR_value"),
               uiOutput("min_reads")
      ),
      
      #DE gane table tab
      tabPanel("Output - Differentially expressed gene table", id = "geneTableTab", value= "geneTableTab",  fluidRow(
        downloadButton("downloadDEGeneTable", "Download Table"),
        DT::dataTableOutput("calc_res_values")
        
      )),
      
      #PCA plot tab
      tabPanel("PCA Plot", id = "pcaPlot", value= "pcaPlotTab", fluidRow(
        radioButtons("PCAplot_labels", label = "Labels",
                     choices = list("No labels" = 1, "Sample names" = 2, "Replicates" = 3), 
                     selected = 1),
        plotOutput("PCA_plot", height = "500", width="500")
    
      )),
      
      tabPanel("Gene read counts", id = "genecountPlotTab", value= "genecountPlotTab", fluidRow(
        textInput("gene_name", "Enter gene name", value = "KRAS"),
        radioButtons("readcountplot_type", label = "Type of plot",
                     choices = list("Boxplot" = 1, "Jitter plot" = 2), 
                     selected = 1),
        radioButtons("readcountplot_labels", label = "Labels",
                     choices = list("No labels" = 1, "Sample names" = 2, "Replicates" = 3), 
                     selected = 1),
        plotOutput("genecount_plot", height = "500", width="500")
        
      )),
      
      tabPanel("Volcano Plot", id = "volcanoPlotTab", value= "volcanoPlotTab", fluidRow(
        sliderInput("FCcutoff", 
                    "Log2 fold change cutoff", 
                    min = 0,
                    max = 10, 
                    value = c(2)),
        textInput("padjcutoff", "Adjusted p value cutoff", value = "0.05", width = NULL,
                  placeholder = NULL),
        plotOutput("volcanoPlot", height = "800", width="100%")
      ))
    )
    
  ), #end dashboardBody
  
  skin = "black"
  
) #end dashboardPage