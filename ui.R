# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# version: 0.62

# added:
# +1 to all reads; avoid 0 read count errors
# multiple comparisons
# show >2 conditions on PCA plot
# adjust results based on different conditions
# options for plots to show labels
# boxplot or jitter plot option for read count plot
# use ggrepel for plot labels so labels don't overlap
# automatically install required packages if not already installed

# bugs"
# PCA, gene count, volcano plots don't auto-update to new dds after changing treatment condition factor level

ui <- dashboardPage(
  dashboardHeader(title = "VisualRNAseq v0.61"),
  
  dashboardSidebar(
    
    
  ), #end dashboard Sidebar
  
  dashboardBody(
    
    tabsetPanel(
      
      tabPanel("Load data", fluidRow(
        
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
      
      tabPanel("Select experiment settings", 
               uiOutput("control_condslist"),
               uiOutput("treatment1_condslist"),
               uiOutput("FDR_value"),
               uiOutput("min_reads")
      ),
      
      tabPanel("Output - Differentially expressed gene table", fluidRow(
        downloadButton("downloadDEGeneTable", "Download Table"),
        DT::dataTableOutput("calc_res_values")
        
      )),
      
      tabPanel("PCA Plot", fluidRow(
        checkboxInput("PCAplot_show_labels", "Show labels", TRUE),
        plotOutput("PCA_plot", height = "500", width="500")
    
      )),
      
      tabPanel("Gene read counts", fluidRow(
        textInput("gene_name", "Enter Gene name", value = "KRAS"),
        radioButtons("readcountplot_type", label = "Type of plot",
                     choices = list("Boxplot" = 1, "Jitter plot" = 2), 
                     selected = 1),
        checkboxInput("readcountplot_show_labels", "Show labels", TRUE),
        plotOutput("genecount_plot", height = "500", width="500")
        
      )),
      
      tabPanel("Volcano Plot", fluidRow(
        sliderInput("FCcutoff", 
                    "Log2 fold change cutoff", 
                    min = 0,
                    max = 10, 
                    value = c(2)),
        textInput("padjcutoff", "Adjusted p value cutoff", value = "0.05", width = NULL,
                  placeholder = NULL),
        plotlyOutput("volcanoPlot", height = "800", width="100%")
      ))
    )
    
  ), #end dashboardBody
  
  skin = "black"
  
) #end dashboardPage