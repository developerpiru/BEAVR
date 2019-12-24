# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# See Github for more info & ReadMe: https://github.com/developerpiru/VisualRNAseq
app_version = "0.71.3"

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
# filter data table
# use filtered data table for volcano plot

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
installReqs("scatterD3", bioc = FALSE)
installReqs("ggplot2", bioc = FALSE)
installReqs("ggrepel", bioc = FALSE)
installReqs("data.table", bioc = FALSE)
installReqs("DT", bioc = FALSE)
installReqs("BiocManager", bioc = FALSE)
installReqs("DESeq2", bioc = TRUE)
installReqs("vsn", bioc = TRUE)
installReqs('apeglm', bioc = TRUE)
installReqs('org.Hs.eg.db', bioc = TRUE)
installReqs('EnhancedVolcano', bioc = TRUE)
installReqs('ggpubr', bioc = FALSE)
installReqs('shinyjqui', bioc = FALSE)

#load required libraries
library("shiny")
library("shinydashboard")
library("plotly")
library("scatterD3")
library("ggplot2")
library("ggrepel")
library("data.table")
library("DT")
library("DESeq2")
library("vsn")
library('apeglm')
library('org.Hs.eg.db')
library('EnhancedVolcano')
library("gridExtra")
library("ggpubr")
library("shinyjqui")

#increase max file size to 100MB
options(shiny.maxRequestSize = 100*1024^2)

shinyServer(function(input, output, session) {
  
  #---BEGIN DATA INPUT---#
  
  #reactive to get and store raw reads data
  #upload read count file
  cts <- reactive({
    
    req(input$rawreadsfile)
    
    #store in rawreadsdata variable
    rawreadsdata <- read.csv(input$rawreadsfile$datapath,
                    header = input$header1,
                    sep = input$sep1)
    
    rownames(rawreadsdata) <- rawreadsdata$gene_id
    rawreadsdata <- rawreadsdata[,-1]
    
    #increment all reads by 1 to avoid 0 read count errors
    rawreadsdata <- rawreadsdata + 1
    
    return(rawreadsdata)
    
  })
  
  #reactive to get and store coldata table
  #upload column data file
  coldata <- reactive({
    
    req(input$coldatafile)
    
    #store in coldata variable
    coldata <- read.csv(input$coldatafile$datapath,
                        header = input$header2,
                        sep = input$sep2)
    
    return(coldata)
    
  })
  
  #---END DATA INPUT---#
  
  #get control condition from list
  output$control_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("control_condslist", "Choose control condition", temp_condslist)
  })
  
  #get treatment condition from list
  output$treatment1_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("treatment1_condslist", "Choose treatment condition", temp_condslist)
  })
  
  #get false discovery rate from user
  output$FDR_value <- renderUI({
    numericInput("FDRvalue", "False Discovery Rate %",value = 10)
  })
  
  #get minimum read count values to keep from user
  output$min_reads <- renderUI({
    numericInput("min_reads_value", "Drop genes with reads below:",value = 10)
  })
  
  #deprecated function!!
  conditionpicker <- reactive({
    
    temp_coldata <- coldata()
    
    condslist <- unique(temp_coldata[,2])
    
    return(condslist)
    
  })
  
  #DEBUG - function to check coldata against read count table column names
  coldatacompare <- reactive({
    
    #cts <- output$rawreadstable
    #coldata <- output$coldatatable
    
    temp_coldata <- coldata()
    temp_cts <- cts()
    
    #check1 <- all(rownames(temp_coldata) %in% colnames(temp_cts))
    
    check1 <- all(rownames(temp_coldata) %in% colnames(temp_cts))
    check2 <- all(rownames(temp_coldata) == colnames(temp_cts))
    temp_cts <- temp_cts[, rownames(temp_coldata)]
    check3 <- all(rownames(temp_coldata) == colnames(temp_cts))
    
    return(check3)
  })
  
  #DEBUG - check coldata and read count tables for matching row\column names
  output$coldatachecker <- renderText({
    
    coldatacompare()
    
  })
  
  #---BEGIN CALCULATIONS---#
  
  #create DESeq2 data set (dds)
  calc_get_dds <- reactive({
    
    #get values
    temp_cts <- cts()
    temp_coldata <- coldata()
    
    control_factor <- input$control_condslist
    treatment1_factor <- input$treatment1_condslist
    treatment2_factor <- input$treatment2_condslist
    treatment3_factor <- input$treatment3_condslist
    
    #construct a DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = temp_cts, colData = temp_coldata, design = ~ condition)
    #dds
    
    #pre-filter dds table to only keep genes that have at least min_reads_value reads set by user
    keep <- rowSums(counts(dds)) > input$min_reads_value
    dds <- dds[keep,]
    
    #OLD METHOD - Setting the factor level
    #dds$condition <- factor(dds$condition, levels = c(control_factor, treatment1_factor))
    
    #NEW METHOD - Setting the factor level
    dds$condition <- relevel(dds$condition, ref = control_factor)
    
    #Differential expression analysis
    dds <- DESeq(dds)
    
    #write.csv(as.data.frame(dds), file='dds-output.csv')
    
    return(dds)
    
  })
  
  #calculate results from dds
  #calculates LFC and FDR
  calc_res <- reactive({
    
    #Update progress bar
    totalSteps = 8 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    control_factor <- input$control_condslist
    treatment1_factor <- input$treatment1_condslist
    
    #build LFC argument based define experimental conditions
    LFC_coef <- paste("condition_", treatment1_factor, sep="")
    LFC_coef <- paste(LFC_coef, control_factor, sep="_vs_")
    
    FDR_aplha <- (input$FDRvalue)/100
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Performing DESeq2 calculations..."))
    
    #calcate dds values
    dds <<- calc_get_dds()
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Determing differential expression..."))
    
    #OLD WAY - set the contrasts for comparisons
    #res <<- results(dds, contrast=c("condition",control_factor, treatment1_factor))
    
    #new way to set multifactor comparisons - log2FC[final/initial]
    res <<- results(dds, contrast=c("condition", treatment1_factor, control_factor))
    
    #Log fold change shrinkage for visualization and ranking
    #resultsNames(dds)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Using 'apeglm' for LFC shrinkage..."))
    
    #adjust conditions based on contrast
    resLFC <<- lfcShrink(dds, coef=LFC_coef, type="apeglm")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Ordering by p values..."))
    
    #order results by smallest p value
    resOrdered <<- res[order(res$pvalue),]
    
    #get some basic tallies using the summary function.
    #summary(res)
    
    #many adjusted p-values are less than 0.1
    #sum(res$padj < 0.1, na.rm=TRUE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Applying FDR correction..."))
    
    #filter res based on FDR cut off (10% = alpha of 0.1)
    resFDR <<- results(dds, alpha=FDR_aplha)
    #summary(res10)
    #sum(res10$padj < FDR_aplha, na.rm=TRUE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Mapping ENSEMBL names to readable gene symbols..."))
    
    #map ensembl symbols to gene ids
    res$GeneID <<- mapIds(org.Hs.eg.db,keys=rownames(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    resFDR$GeneID <<- mapIds(org.Hs.eg.db,keys=rownames(resFDR),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    resOrdered$GeneID <<- mapIds(org.Hs.eg.db,keys=rownames(resOrdered),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    
    #make a separate table ensembl symbols and gene ids
    listofgenes <<- as.data.frame(res)
    listofgenes <<- subset(listofgenes, select = c(GeneID))
    listofgenes$ENSEMBL <<- rownames(listofgenes)
    
    #write files
    #write.csv(as.data.frame(resOrdered), file='resOrdered-all-genes.csv')
    #write.csv(as.data.frame(resFDR), file='resFDR.csv')
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Preparing to display table..."))
    
    resFDR <<- resFDR[,c(7,1:6)]
    
    return(resFDR)
    #return(resOrdered)
    
  })
  
  #---END CALCULATIONS---#
  
  #function to process filtering of dds table
  filterTable <- reactive({
    
    filteredDDSTable <<- as.data.frame(calc_res())
    
    if (input$tableFilterONOFF == TRUE){
      filteredDDSTable <<- subset(filteredDDSTable, 
                                  filteredDDSTable$log2FoldChange >= input$tableMinLog2FC &
                                    filteredDDSTable$log2FoldChange <= input$tableMaxLog2FC &
                                    filteredDDSTable$pvalue >= input$tableMinPvalue &
                                    filteredDDSTable$pvalue <= input$tableMaxPvalue &
                                    filteredDDSTable$padj >= input$tableMinPadj &
                                    filteredDDSTable$padj <= input$tableMaxPadj)
    }
    
    
    return(filteredDDSTable)
  })
  
  #output calculated dds + FDR table
  #function to show table
  output$calc_res_values <- DT::renderDataTable({
    withProgress(message = 'Performing calculations...', value = 1, min = 1, max = 100, {
      #as.data.frame(calc_res())
      filterTable()
    })
  })
  
  #download DE gene table
  output$downloadDEGeneTable <- downloadHandler(
    filename = function() {
      paste("Differentially Expressed Genes.csv")
    },
    content = function(file) {
      #write.csv(as.data.frame(calc_res()), file, row.names = FALSE)
      write.csv(as.data.frame(filteredDDSTable), file, row.names = FALSE)
    }
  )
  
  #call function to show PCA plot
  output$PCA_plot = renderPlot({
    withProgress(message = 'Generating PCA plot...', value = 1, min = 1, max = 100, {
      do_PCA_plot()
    })
  })
  
  #function to draw PCA plot
  do_PCA_plot <- reactive({
    
    #Update progress bar
    totalSteps = 3 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    vsd <<- vst(dds, blind=FALSE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Plotting..."))
    
    pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    #generate the plot
    p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
      geom_point(size = input$pcaPointSize) +
      labs(shape="Replicate", colour="Condition") +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      theme_classic() +
      theme(axis.text.x = element_text(colour="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.title.y = element_text(colour="black",size=input$pcaFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.title = element_text(colour="black",size=input$pcaFontSize_legend_title,angle=0,hjust=.5,vjust=.5,face="plain"),
            legend.text =  element_text(colour="black",size=input$pcaFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
            legend.text.align = 0,
            text = element_text(size=input$pcaFontSize_x_axis))
      
      
      #theme(text = element_text(size=input$pcaFontSize))
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #show labels for points as determined by user
    if (input$PCAplot_labels == 1){
      #no labels
      p <- p
    } else if (input$PCAplot_labels == 2){
      #sample names as labels
      p <- p + geom_text_repel(size=input$pcaLabelFontSize, nudge_x=0.1, nudge_y=0.1, segment.color=NA, aes(label=rownames(pcaData)))
      aes(shape=rownames(d))
    } else if (input$PCAplot_labels == 3){
      #replicate names as labels
      p <- p + geom_text_repel(size=input$pcaLabelFontSize, nudge_x=0.1, nudge_y=0.1, segment.color=NA, aes(label=replicate))
    }

    #return the plot
    print(p)
    
  })
  
  #make volcano plot highlight genes that have an FDR cutoff and Log2FC cutoff as determined by the user (input$padjcutoff and input$FCcutoff)
  output$volcanoPlot = renderPlot({
    withProgress(message = 'Generating volcano plot...', value = 1, min = 1, max = 100, {
      do_volcano_plot()
    })
  })
  
  #function to plot multiple gene read counts for user defined genes
  do_volcano_plot <- reactive({
    
    #Update progress bar
    totalSteps = 2 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get RNAseq data
    if (input$volcanoFilterONOFF == TRUE){
      RNAseqdatatoplot <<- as.data.frame(filteredDDSTable)
    } else {
      RNAseqdatatoplot <<- as.data.frame(resFDR)
    }
    
    if (input$volcanoPvalueType == TRUE){
      pvaluetype = input$padjcutoff
    } else {
      pvaluetype = input$pvalue
    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    if (input$volcanoCutoffLines == TRUE){
      p <- EnhancedVolcano(RNAseqdatatoplot,
                           title = paste(input$control_condslist, input$treatment1_condslist, sep = " vs. "),
                           subtitle = "",
                           lab = RNAseqdatatoplot$GeneID,
                           x = "log2FoldChange",
                           y = "padj",
                           pCutoff = as.numeric(paste(pvaluetype)),
                           #pCutoff = as.numeric(input$padjcutoff),
                           FCcutoff = as.numeric(input$FCcutoff),
                           titleLabSize = input$volcanoFontSize_plot_title,
                           axisLabSize = input$volcanoFontSize_xy_axis,
                           pointSize = input$volcanoPointSize,
                           labSize  = input$volcanoFontSize_label,
                           legendLabSize = input$volcanoFontSize_legend_title,
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE,
                           legendPosition = "bottom",
                           legendLabels = c('Not Significant', expression(Log[2]~FC~only), "p-value only", expression(p-value~and~log[2]~FC)),
                           cutoffLineType = "longdash",
                           cutoffLineCol = 'black',
                           labCol = 'black',
                           caption = ""
      )

    } else {
      p <- EnhancedVolcano(RNAseqdatatoplot,
                           title = paste(input$control_condslist, input$treatment1_condslist, sep = " vs. "),
                           subtitle = "",
                           lab = RNAseqdatatoplot$GeneID,
                           x = "log2FoldChange",
                           y = "padj",
                           pCutoff = as.numeric(input$padjcutoff),
                           FCcutoff = as.numeric(input$FCcutoff),
                           titleLabSize = input$volcanoFontSize_plot_title,
                           axisLabSize = input$volcanoFontSize_xy_axis,
                           pointSize = input$volcanoPointSize,
                           labSize  = input$volcanoFontSize_label,
                           legendLabSize = input$volcanoFontSize_legend_title,
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE,
                           legendPosition = "bottom",
                           cutoffLineType = "blank",
                           labCol = 'black',
                           caption = ""
      )
    }
    
    
    #return the plot
    print(p)
  })
  
  #multi gene count plot
  output$multi_genecount_plot1 = renderPlot({
    withProgress(message = 'Generating read count plots...', value = 1, min = 1, max = 100, {
     do_multi_genecount_plot()
    })
  })
  
  #function to plot gene counts for user defined genes
  do_multi_genecount_plot <- reactive({
    
    #Update progress bar
    totalSteps = 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get multi gene names from user input
    multi_gene_names <- unlist(strsplit(toupper(input$multi_gene_name), ","))
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Plotting..."))

    #initialize variables to run through and generate all the gene count plots
    p = list()
    i = 0
    
    #loop through and generate the plots for gene names entered
    for (val in multi_gene_names){

      i = i+1
      
      #get data for selected gene from dds data matrix
      d <- plotCounts(dds, gene=listofgenes[which(listofgenes$GeneID==toupper(val)),2], intgroup=c("condition", "replicate"), returnData=TRUE)
      
      #generate each plot
      p[[i]] <- ggplot(d, aes(x=condition, y=count, color=condition)) +
        ggtitle(toupper(val)) +
        geom_point(size = input$genecountPointSize) +
        labs(shape="Replicate", colour="Condition") +
        xlab("") +
        ylab("Normalized read count") +
        theme_classic() +
        theme(axis.text.x = element_text(colour="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
              axis.text.y = element_text(colour="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),  
              axis.title.x = element_text(colour="black",size=0),
              axis.title.y = element_text(colour="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
              plot.title = element_text(colour="black",size=input$multi_genecountFontSize_plot_title,angle=0,hjust=.5,vjust=.5,face="plain"),
              #legend.title = element_text(colour="black",size=input$multi_genecountFontSize_legend_title,angle=0,hjust=.5,vjust=.5,face="plain"),
              legend.title = element_blank(),
              legend.text =  element_text(colour="black",size=input$multi_genecountFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
              legend.text.align = 0)
      
      #change plot type to boxplot or jitter plot based on user selection
      if (input$multi_readcountplot_type == 1){
        p[[i]] <- p[[i]] + geom_boxplot() +
          #hide points
          geom_point(size = -1)
      } else {
        p[[i]] <- p[[i]] + geom_jitter(size=input$multi_genecountPointSize, width=0, height=0) +
          aes(shape=replicate)
      }
      
      #show labels for points as determined by user
      if (input$multi_readcountplot_labels == 1){
        #no labels
        p[[i]] <- p[[i]] + labs(shape="Replicate", colour="Condition")
      } else if (input$multi_readcountplot_labels == 2){
        #sample names as labels
        p[[i]] <- p[[i]] + geom_text_repel(size=input$multi_genecountLabelFontSize, nudge_x=0.1, nudge_y=0.1, segment.color=NA, aes(label=rownames(d))) +
          labs(shape="Replicate", colour="Condition")
        #aes(shape=rownames(d))
      } else if (input$multi_readcountplot_labels == 3){
        #replicate names as labels
        p[[i]] <- p[[i]] + geom_text_repel(size=input$multi_genecountLabelFontSize, nudge_x=0.1, nudge_y=0.1, segment.color=NA, aes(label=replicate)) +
          labs(shape="Replicate", colour="Condition")
      }
      
      #adjust legend based on user selection
      if (input$multi_genecountShowLegends == 1){ # hide legend on all plots
        p[[i]] <- p[[i]] + theme(legend.position = "none")
      } else if (input$multi_genecountShowLegends == 2 | input$multi_genecountShowLegends == 3){ # show legend on all plots or show one common legend
        p[[i]] <- p[[i]] + theme(legend.position = paste(input$multi_genecountLegendPosition))
      } else if (input$multi_genecountShowLegends == 3){ # show one common legend
        
      }
      
      #turn off y-axis unless the plot is the first one in a row, based on user selection
      if (input$multi_genecountSharedYAxis == TRUE){
        if ((i-1)%%input$multi_genecountGridColumns == 0 | i==1){
          p[[i]] <- p[[i]] + theme(axis.title.y = element_text(colour="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"))
        } else {
          p[[i]] <- p[[i]] + theme(axis.title.y = element_blank())
        }
      } else {
        p[[i]] <- p[[i]] + theme(axis.title.y = element_text(colour="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"))
      }
      
      #set y-axis to log scale depending on user selection
      if (input$multi_log10scale == TRUE){
        p[[i]] <- p[[i]] + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                         labels = trans_format("log10", math_format(10^.x)))
      }
        

    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    if (input$multi_genecountShowLegends == 1){
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = FALSE, legend="none"))
    } else if (input$multi_genecountShowLegends == 3){
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = TRUE, legend=paste(input$multi_genecountLegendPosition)))
    } else {
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = FALSE, legend=paste(input$multi_genecountLegendPosition)))
    }

  })

})



