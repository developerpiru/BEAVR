# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# version: 0.61

# added:
# +1 to all reads; avoid 0 read count errors
# multiple comparisons
# show >2 conditions on PCA plot
# adjusdt results based on different conditions
# labels for read count plot

# bugs"
# PCA, gene count, volcano plots don't auto-update to new dds after changing treatment condition factor level

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

#increase max file size to 50MB
options(shiny.maxRequestSize = 50*1024^2)

shinyServer(function(input, output, session) {
  
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
  
  output$control_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("control_condslist", "Choose control condition", temp_condslist)
  })
  
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
  
  conditionpicker <- reactive({
    
    temp_coldata <- coldata()
    
    condslist <- unique(temp_coldata[,2])
    
    return(condslist)
    
  })
  
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
  
  #output calculated dds + FDR table
  output$calc_res_values <- DT::renderDataTable({
    withProgress(message = 'Performing calculations...', value = 1, min = 1, max = 100, {
      as.data.frame(calc_res())
    })
  })
  
  #download all genes table
  output$downloadDEGeneTable <- downloadHandler(
    filename = function() {
      paste("Differentially Expressed Genes.csv")
    },
    content = function(file) {
      write.csv(as.data.frame(calc_res()), file, row.names = FALSE)
    }
  )
  
  #PCA plot
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
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating..."))
    
    pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
    
    print(p)
    
  })
  
  #gene count plot
  output$genecount_plot = renderPlot({
    withProgress(message = 'Generating read count plot...', value = 1, min = 1, max = 100, {
      do_genecount_plot()
    })
  })
  
  #function to plot gene counts for user defined genes
  do_genecount_plot <- reactive({
    
    #Update progress bar
    totalSteps = 2 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get gene name from user
    genename = toupper(input$gene_name)
    
    #get data for selected gene from dds data matrix
    d <- plotCounts(dds, gene=listofgenes[which(listofgenes$GeneID==genename),2], intgroup=c("condition", "replicate"), returnData=TRUE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    p <- ggplot(d, aes(x=condition, y=count, color=condition, shape=replicate)) + 
      geom_point(size=3, position=position_jitter(w=0.2,h=0)) +
      geom_text(aes(label = replicate), position = position_nudge(y = -10)) +
      #scale_y_log10(breaks=c(25,100,400)) +
      ggtitle(genename) +
      xlab("") +
      ylab("Normalized count") +
      theme_bw() +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
            plot.title = element_text(colour="grey20",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
    
    #return the plot
    print(p)
  })
  
  output$volcanoPlot = renderPlotly({
    
    #get RNAseq data
    RNAseqdatatoplot <- as.data.frame(resFDR)
    
    # add a grouping column; default value is "not significant"
    RNAseqdatatoplot["group"] <- "NotSignificant"
    
    # for our plot, we want to highlight 
    # FDR cutoff from text input padjcutoff (input$padjcutoff)
    # Fold Change cutoff from input slider FCcutoff (input$FCcutoff)
    
    # change the grouping for the entries with significance but not a large enough Fold change
    RNAseqdatatoplot[which(RNAseqdatatoplot['padj'] < input$padjcutoff & abs(RNAseqdatatoplot['log2FoldChange']) < input$FCcutoff ),"group"] <- "Significant"
    
    # change the grouping for the entries a large enough Fold change but not a low enough p value
    RNAseqdatatoplot[which(RNAseqdatatoplot['padj'] > input$padjcutoff & abs(RNAseqdatatoplot['log2FoldChange']) > input$FCcutoff ),"group"] <- "FoldChange"
    
    # change the grouping for the entries with both significance and large enough fold change
    RNAseqdatatoplot[which(RNAseqdatatoplot['padj'] < input$padjcutoff & abs(RNAseqdatatoplot['log2FoldChange']) > input$FCcutoff ),"group"] <- "Significant&FoldChange"
    
    # make gene names a URL to genecards.com
    RNAseqdatatoplot$GENEurl <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", RNAseqdatatoplot$GeneID, "'>", RNAseqdatatoplot$GeneID, "</a>")
    
    # Find and label the top peaks..
    top_peaks <- RNAseqdatatoplot[with(RNAseqdatatoplot, order(log2FoldChange, padj)),][1:10,]
    top_peaks <- rbind(top_peaks, RNAseqdatatoplot[with(RNAseqdatatoplot, order(-log2FoldChange, padj)),][1:10,])
    
    # Add gene labels for all of the top genes we found 
    # here we are creating an empty list, and filling it with entries for each row in the dataframe
    # each list entry is another list with named items that will be used by Plot.ly
    a <- list()
    for (i in seq_len(nrow(top_peaks))) {
      m <- top_peaks[i, ]
      a[[i]] <- list(
        x = m[["log2FoldChange"]],
        y = -log10(m[["padj"]]),
        text = m[["GeneID"]],
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40
      )
    }
    
    x <- list(
      title = "Log2 Fold change"
    )
    y <- list(
      title = "FDR"
    )
    
    p <- plot_ly(data = RNAseqdatatoplot, 
                 x = RNAseqdatatoplot$log2FoldChange, 
                 y = -log10(RNAseqdatatoplot$padj), 
                 text = RNAseqdatatoplot$GENEurl, 
                 mode = "markers", 
                 color = RNAseqdatatoplot$group) %>% 
      layout(title ="Volcano Plot", xaxis = x, yaxis = y) %>%
      layout(annotations = a)

    
  })

})



