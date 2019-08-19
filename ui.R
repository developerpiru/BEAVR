#setwd("D:/Google Drive/PhD Stuff/RNAseqGUI-R/RNAseqGUI_v5/")
#runApp("shinyapp", host = "0.0.0.0", port = 80)

ui <- dashboardPage(
  dashboardHeader(title = "VisualRNAseq v0.5"),
  
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
               #uiOutput("treatment2_condslist"),
               #uiOutput("treatment3_condslist"),
               uiOutput("FDR_value")
      ),
      
      tabPanel("Output - Differentially expressed gene table", fluidRow(
        downloadButton("downloadDEGeneTable", "Download Table"),
        DT::dataTableOutput("calc_res_values")
        
      )),
      
      tabPanel("PCA Plot", fluidRow(
        plotOutput("PCA_plot", height = "500", width="500")
        
      )),
      
      tabPanel("Gene counts", fluidRow(
        textInput("gene_name", "Enter Gene name", value = "KRAS"),
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
      )),
      
      tabPanel("raw reads", fluidRow(
        tableOutput("rawreadstable")
        
      )),
      
      tabPanel("coldata", tableOutput("coldatatable")
      ), 

      tabPanel("conds list", tableOutput("condstable")
      )
    )
    
  ), #end dashboardBody
  
  skin = "black"
  
) #end dashboardPage