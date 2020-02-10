# BEAVR: A Browser-based tool for the Exploration And Visualization of RNA-seq data
# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# See Github for more info & ReadMe: https://github.com/developerpiru/BEAVR

app_version = "1.0.1"

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

#increase max file size to 1000MB
options(shiny.maxRequestSize = 1000*1024^2)

shinyServer(function(input, output, session) {
  
  #---BEGIN DATA INPUT---#
  
  #reactive to get and store raw reads data
  #upload read count file
  cts <- reactive({
    
    req(input$rawreadsfile)
    
    #store in rawreadsdata variable
    rawreadsdata <- read.csv(input$rawreadsfile$datapath,
                             header = TRUE,
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
                        header = TRUE,
                        sep = input$sep2)
    
    return(coldata)
    
  })
  
  #---END DATA INPUT---#
  
  #---START DYNAMIC EXPERIMENT SETTINGS---#
  
  #get control condition from list
  output$control_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("control_condslist", "Choose control condition", temp_condslist, selected = 1)
  })
  
  #get treatment condition from list
  output$treatment1_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("treatment1_condslist", "Choose treatment condition", temp_condslist, selected = 2)
  })
  
  #get false discovery rate from user
  output$FDR_value <- renderUI({
    numericInput("FDRvalue", "False Discovery Rate %",value = 10)
  })
  
  #get minimum read count values to keep from user
  output$min_reads <- renderUI({
    numericInput("min_reads_value", "Drop genes with reads below:",value = 10)
  })
  
  #---END DYNAMIC EXPERIMENT SETTINGS---#
  
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
    temp_coldata <<- coldata()
    
    #get number of unique conditions
    num_conditions <<- length(unique(temp_coldata[,2]))
    
    #prepare list of condition names
    conds_names <<- levels(temp_coldata[,2])
    
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
    
    #set flag to FALSE for first run indicator -- used in multi gene count function
    first_run_flag1 <<- TRUE
    first_run_flag2 <<- TRUE
    
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

    #new way to set multifactor comparisons - log2FC[final/initial]
    res <<- results(dds, contrast=c("condition", treatment1_factor, control_factor))
    
    #Log fold change shrinkage for visualization and ranking

    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Performing LFC shrinkage..."))
    
    #adjust conditions based on contrast
    if (input$shrinkage_method == 1){
      resLFC <<- lfcShrink(dds, coef=LFC_coef, type="apeglm")
    } else {
      resLFC <<- lfcShrink(dds, coef=LFC_coef, type="normal")
    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Ordering by p values..."))
    
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
    #check which reference genome was selected by the user and translate ENSEMBL IDs using the correct one
    if (input$ref_genome_organism == 1){
      # 1 = human
      res$GeneID <<- mapIds(org.Hs.eg.db,keys=rownames(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
      resFDR$GeneID <<- mapIds(org.Hs.eg.db,keys=rownames(resFDR),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    } else if (input$ref_genome_organism == 2){
      # 2 = mouse
      res$GeneID <<- mapIds(org.Mm.eg.db,keys=rownames(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
      resFDR$GeneID <<- mapIds(org.Mm.eg.db,keys=rownames(resFDR),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    }

    #make a separate table of ensembl symbols and gene ids
    #required for gene read counts plot function
    #make a copy of res table
    listofgenes <<- as.data.frame(res)
    #keep only the GeneID column; drop everything else
    listofgenes <<- subset(listofgenes, select = c(GeneID))
    #save ENSEMBL IDs in new ENSEMBL column
    listofgenes$ENSEMBL <<- rownames(listofgenes)
    
    #write files
    #write.csv(as.data.frame(resFDR), file='resFDR.csv')
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Preparing to display table..."))
    
    resFDR <<- resFDR[,c(7,1:6)]
    
    #update numericInputs for data filtering with the min/max values from resFDR
    updateNumericInput(session, "log2FC_min", value = min(resFDR$log2FoldChange, na.rm=T))
    updateNumericInput(session, "log2FC_max", value = max(resFDR$log2FoldChange, na.rm=T))
    updateNumericInput(session, "pvalue_min", value = min(resFDR$pvalue, na.rm=T))
    updateNumericInput(session, "pvalue_max", value = max(resFDR$pvalue, na.rm=T))
    updateNumericInput(session, "padj_min", value = min(resFDR$padj, na.rm=T))
    updateNumericInput(session, "padj_max", value = max(resFDR$padj, na.rm=T))
    updateNumericInput(session, "baseMean_min", label = "Min", value = min(resFDR$baseMean, na.rm=T))
    updateNumericInput(session, "baseMean_max", label = "Max", value = max(resFDR$baseMean, na.rm=T))
    updateNumericInput(session, "lfcSE_min", label = "Min", value = min(resFDR$lfcSE, na.rm=T))
    updateNumericInput(session, "lfcSE_max", label = "Max", value = max(resFDR$lfcSE, na.rm=T))
    updateNumericInput(session, "stat_min", value = min(resFDR$stat, na.rm=T))
    updateNumericInput(session, "stat_max", value = max(resFDR$stat, na.rm=T))
    return(resFDR)
    
  })
  
  #---END CALCULATIONS---#
  
  #output calculated dds + FDR table
  #function to show table
  output$calc_res_values <- DT::renderDataTable({
    withProgress(message = 'Performing calculations...', value = 1, min = 1, max = 100, {
      unfilteredTable <<- as.data.frame(calc_res())
      
      #check if the Enable filtering checkbox is checked; if so, enable filtering as below
      #save filtered table to filteredTable global variable
      #if not, save the unfiltered table to filteredTable global variable
      #that way, functions can be simplified and the same variable (filteredTable) can carry both table types
      if (input$filterTableEnabled == TRUE){
        filteredTable <<- subset(unfilteredTable, 
                                 baseMean >= input$baseMean_min &
                                   baseMean <= input$baseMean_max &
                                   log2FoldChange >= input$log2FC_min &
                                   log2FoldChange <= input$log2FC_max &
                                   lfcSE >= input$lfcSE_min &
                                   lfcSE <= input$lfcSE_max &
                                   stat >= input$stat_min &
                                   stat <= input$stat_max &
                                   pvalue >= input$pvalue_min &
                                   pvalue <= input$pvalue_max &
                                   padj >= input$padj_min &
                                   padj <= input$padj_max
        )
      } else {
        filteredTable = unfilteredTable
      }
      
      #show table
      filteredTable
    })
    
  })
  
  #download DE gene table
  output$downloadDEGeneTable <- downloadHandler(
    filename = function() {
      paste("Differentially Expressed Genes.csv")
    },
    content = function(file) {
      #write.csv(as.data.frame(calc_res()), file, row.names = FALSE)
      write.csv(as.data.frame(filteredTable), file, row.names = FALSE)
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
    
    #transform data using variance stabilization method
    #nsub=nrow(dds)
    #vsd <<- vst(dds, blind=FALSE)
    vsd <<- varianceStabilizingTransformation(dds, blind = FALSE)
    #vsd <<- rlogTransformation(dds, blind = FALSE, fitType = "parametric")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Plotting..."))

    #run only if first_run_flag2 boolean is TRUE; meaning ui is being initialized
    if (first_run_flag2 == TRUE){
      #call function to loop through and create more colour widgets for each condition in the experiment
      #widget names are pcaColorX, where X is an integer
      #selector is the target div tag container in ui
      dynamic_colorgen(widget_name = "pcaColor", selector = "#pcaColorbox")
    }
    
    #set first run flag to false so color widgets are no longer made
    first_run_flag2 <<- FALSE
    
    #vector to save colours
    multi_colorslist <- NULL
    
    #loop through and get the colours that the user choses
    for (count in 1:num_conditions){
      multi_colorslist[count] <- input[[paste0('pcaColor', count)]]
    }
    
    #plot transformed data in PCA
    pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    #generate the plot
    p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
      geom_point(size = input$pcaPointSize) +
      scale_color_manual(values = multi_colorslist) + 
      labs(shape="Replicate", color="Condition") +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      theme_classic() +
      theme(axis.text.x = element_text(color="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(color="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(color="black",size=input$pcaFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.title.y = element_text(color="black",size=input$pcaFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.title = element_text(color="black",size=input$pcaFontSize_legend_title,angle=0,hjust=.5,vjust=.5,face="plain"),
            legend.text =  element_text(color="black",size=input$pcaFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
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
  
  #call function to show sample clustering plot
  output$sampleClustering_plot = renderPlot({
    withProgress(message = 'Generating sample clustering heatmap...', value = 1, min = 1, max = 100, {
      do_sampleClustering_plot()
    })
  })
  
  #function to calculate sample clustering
  do_sampleClustering_plot <- reactive({
    
    #Update progress bar
    totalSteps = 3 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get some user defined values
    clust_dist = input$sampleClustering_method
    if (input$sampleClustering_cellNums == FALSE){
      display_cellVals = FALSE
    } else if (input$sampleClustering_cellNums == "%.2f"){
      display_cellVals = TRUE
      cellVals_format = "%.2f"
    } else if (input$sampleClustering_cellNums == "%.1e"){
      display_cellVals = TRUE
      cellVals_format = "%.1e"
    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating..."))
    
    #transform data using variance stabilization method
    #vsd <<- vst(dds, blind=FALSE)
    vsd <<- varianceStabilizingTransformation(dds, blind = FALSE)
    
    #transpose the data
    sampleDists <- dist(t(assay(vsd)))
    
    #generate heatmap for sample clustering
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
    colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
    map_colors <- colorRampPalette( rev(brewer.pal(9, paste(input$sampleClustering_mapColor))) )(255)
    p <- pheatmap(sampleDistMatrix,
                  #base fontsize
                  fontsize = input$sampleClustering_fontsize,
                  #display cell values
                  display_numbers = display_cellVals,
                  #cell value format
                  number_format = cellVals_format,
                  #cell value fontsize
                  fontsize_number = input$sampleClustering_fontsize_cellNums,
                  #cell value font color
                  number_color = input$sampleClustering_cellNumsColor,
                  #show legend
                  legend = input$sampleClustering_legend,
                  #tree draw heights for rows and columns
                  treeheight_row = input$sampleClustering_treeHeightRows,
                  treeheight_col = input$sampleClustering_treeHeightCols,
                  #distance measurement method for rows and columns
                  clustering_distance_rows=clust_dist,
                  clustering_distance_cols=clust_dist,
                  #heatmap color
                  col = map_colors,
                  #border color
                  border_color = input$sampleClustering_borderColor)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))

    #return the plot
    print(p)
    
  })
  
  #call function to show count matrix heatmpa
  output$countMatrix_heatmap = renderPlot({
    withProgress(message = 'Generating heatmap...', value = 1, min = 1, max = 100, {
      do_countMatrix_heatmap()
    })
  })
  
  #function to calculate count matrix heatmap
  do_countMatrix_heatmap <- reactive({
    
    #Update progress bar
    totalSteps = 3 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get some user defined values
    clust_dist = input$heatmap_distance
    if (input$heatmap_cellNums == FALSE){
      display_cellVals = FALSE
    } else if (input$heatmap_cellNums == "%.2f"){
      display_cellVals = TRUE
      cellVals_format = "%.2f"
    } else if (input$heatmap_cellNums == "%.1e"){
      display_cellVals = TRUE
      cellVals_format = "%.1e"
    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating..."))
    
    #transform data using variance stabilization method
    vsd <<- varianceStabilizingTransformation(dds, blind = FALSE)
    rld <<- rlogTransformation(dds, blind = FALSE, fitType = "parametric")
    #ntd <<- normTransform(dds)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Subsetting..."))
    
    #get annotations for labels
    annot <- as.data.frame(colData(dds)[,c("condition","replicate")])
    
    
    
    if (input$heatmap_pickTopGenes == TRUE){
      #get the top X genes as defined by user
      genestokeep <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:input$heatmap_numGenes]
      
    } else {
      #get the list of genes entered by the user
      genestomap_HGNC <<- unlist(strsplit(toupper(input$heatmap_GeneNames), ","))
      
      #convert the user-entered gene symbols to ENSEMBL IDS
      if (input$ref_genome_organism == 1){
        # 1 = human
        genestokeep <<- mapIds(org.Hs.eg.db,keys=genestomap_HGNC,column="ENSEMBL",keytype="SYMBOL",multiVals="first")
      } else if (input$ref_genome_organism == 2){
        # 2 = mouse
        genestokeep <<- mapIds(org.Mm.eg.db,keys=genestomap_HGNC,column="ENSEMBL",keytype="SYMBOL",multiVals="first")
      }
    }

    #subset the required genes from the rld transformation data   
    heatmap_data <- as.data.frame(assay(rld)[genestokeep,])
    
    #get the gene names for the subsetted data
    if (input$ref_genome_organism == 1){
      # 1 = human
      heatmap_data$GeneID <- mapIds(org.Hs.eg.db,keys=rownames(heatmap_data),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    } else if (input$ref_genome_organism == 2){
      # 2 = mouse
      heatmap_data$GeneID <- mapIds(org.Mm.eg.db,keys=rownames(heatmap_data),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    }
    
    # #drop any rows that don't have HGNC symbols (have NA instead)
    # heatmap_data <- na.omit(heatmap_data, cols = c("GeneID"))
    
    #set flag to TRUE by default
    show_geneNames = TRUE
    
    #check which option is set for rownames in the heatmap by the user
    if (input$heatmap_showGeneNames == "HGNC"){
      
      #drop any rows that don't have HGNC symbols (have NA instead)
      heatmap_data <- na.omit(heatmap_data, cols = c("GeneID"))
      
      #set rownames to GeneID
      rownames(heatmap_data) <- heatmap_data$GeneID
      
      #if false, rownames are already set to ENSEMBL IDs so don't change anything
    } else if (input$heatmap_showGeneNames == "Hide") {
      #set flag to false
      show_geneNames = FALSE
    }

    #drop GeneID column before sending to pheatmap, which must be numerical data
    heatmap_data <- subset(heatmap_data, select = -c(GeneID))

    #generate heatmap - no row labels
    p <- pheatmap(heatmap_data, 
                  scale = "row", 
                  cluster_rows = input$heatmap_clustRows, 
                  cluster_cols = input$heatmap_clustCols,
                  #distance measurement method for rows and columns
                  clustering_distance_rows=clust_dist,
                  clustering_distance_cols=clust_dist,
                  #clustering method
                  clustering_method = input$heatmap_clustMethod,
                  #tree draw heights for rows and columns
                  treeheight_row = input$heatmap_treeHeightRows,
                  treeheight_col = input$heatmap_treeHeightCols,
                  #show gene names
                  show_rownames = show_geneNames,
                  #base fontsize
                  fontsize = input$heatmap_fontsize,
                  #display cell values
                  display_numbers = display_cellVals,
                  #cell value format
                  number_format = cellVals_format,
                  #cell value fontsize
                  fontsize_number = input$heatmap_fontsize_cellNums,
                  #gene name fontsize
                  fontsize_row = input$heatmap_fontsize_geneNames,
                  #sample name fontsize
                  fontsize_col = input$heatmap_fontsize_sampleNames,
                  #cell value font color
                  number_color = input$heatmap_cellNumsColor,
                  annotation_col=annot,
                  #border color
                  border_color = input$heatmap_borderColor,
                  color = colorRampPalette(c(input$heatmap_lowColor, input$heatmap_midColor, input$heatmap_highColor))(50))

    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))

    #return the plot
    print(p)
    
  })
  
  #make volcano plot highlight genes that have an FDR cutoff and Log2FC cutoff as determined by the user (input$padjcutoff and input$FCcutoff)
  output$volcanoPlot = renderPlot({
    withProgress(message = 'Generating volcano plot...', value = 1, min = 1, max = 100, {
      do_volcano_plot()
    })
  })
  
  #function to plot gene counts for user defined genes
  do_volcano_plot <- reactive({
    
    #Update progress bar
    totalSteps = 2 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get RNAseq data
    #if filtering is enabled, use filtered data
    if (input$filterTableEnabled == TRUE){
      RNAseqdatatoplot <<- as.data.frame(filteredTable)
    } else {
      #otherwise, use unfiltered data
      RNAseqdatatoplot <<- as.data.frame(unfilteredTable)
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
                           pCutoff = as.numeric(input$padjcutoff),
                           FCcutoff = as.numeric(input$FCcutoff),
                           titleLabSize = input$volcanoFontSize_plot_title,
                           axisLabSize = input$volcanoFontSize_xy_axis,
                           pointSize = input$volcanoPointSize,
                           labSize  = input$volcanoFontSize_label,
                           legendLabSize = input$volcanoFontSize_legend_title,
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE,
                           legendPosition = input$volcanoLegendPosition,
                           legendLabels = c('Not Significant', expression(Log[2]~FC~only), "p-value only", expression(p-value~and~log[2]~FC)),
                           cutoffLineType = "longdash",
                           cutoffLineCol = 'black',
                           labCol = 'black',
                           col = c(input$volcano_NSColor, input$volcano_LFCColor, input$volcano_pvalColor, input$volcano_pvalLFCColor),
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
                           legendPosition = input$volcanoLegendPosition,
                           legendLabels = c('Not Significant', expression(Log[2]~FC~only), "p-value only", expression(p-value~and~log[2]~FC)),
                           cutoffLineType = "blank",
                           labCol = 'black',
                           col = c(input$volcano_NSColor, input$volcano_LFCColor, input$volcano_pvalColor, input$volcano_pvalLFCColor),
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
    #need to check which genome is selected and do upper case or lower case based on that
    if (input$ref_genome_organism == 1){
      # 1 = human
      multi_gene_names <- unlist(strsplit(toupper(input$multi_gene_name), ","))
    } else if (input$ref_genome_organism == 2) {
      # 2 = mouse
      multi_gene_names <- unlist(strsplit(input$multi_gene_name, ","))
    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Plotting..."))
  
    
    #run only if first_run_flag2 boolean is TRUE; meaning ui is being initialized
    # if (first_run_flag1 == TRUE){
    #   #call function to loop through and create more colour widgets for each condition in the experiment
    #   #widget names are pcaColorX, where X is an integer
    #   #selector is the target div tag container in ui
    #   dynamic_colorgen(widget_name = "multi_genecountColor", selector = "#multi_genecountColorbox")
    # }
    # 
    # #set first run flag to false so color widgets are no longer made
    # first_run_flag1 <<- FALSE
    
    #vector to save colours
    #multi_colorslist <- NULL
    
    #loop through and get the colours that the user choses
    # for (count in 1:num_conditions){
    #   multi_colorslist[count] <- input[[paste0('multi_genecountColor', count)]]
    # }
    
    #initialize variables to run through and generate all the gene count plots
    p = list()
    i = 0
    
    #loop through and generate the plots for gene names entered
    for (val in multi_gene_names){
      
      i = i+1
      
      #check which genome was selected
      #if human is selected, make sure gene names are upper case
      if (input$ref_genome_organism == 1){
        # 1 = human
        val = toupper(val)
      }
      
      d <- plotCounts(dds, gene=listofgenes[which(listofgenes$GeneID==val),2], intgroup=c("condition", "replicate"), returnData=TRUE)

      #generate each plot
      # p[[i]] <- ggplot(d, aes(x=condition, y=count, color=condition)) +
      #   #p[[i]] <- ggplot(d, x=condition, y=count) +
      #   ggtitle(val) +
      #   geom_point(size = input$genecountPointSize) +
      #   labs(shape="Replicate", color="Condition") +
      #   xlab("") +
      #   ylab("Normalized read count") +
      #   theme_classic() +
      #   theme(axis.text.x = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
      #         axis.text.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),  
      #         axis.title.x = element_text(color="black",size=0),
      #         axis.title.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
      #         plot.title = element_text(color="black",size=input$multi_genecountFontSize_plot_title,angle=0,hjust=.5,vjust=.5,face="plain"),
      #         legend.title = element_blank(),
      #         legend.text =  element_text(color="black",size=input$multi_genecountFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
      #         legend.text.align = 0)
      # p[[i]] <- p[[i]] + stat_compare_means(label = "p.signif", method = "t.test", paired = TRUE, ref.group = input$control_condslist , label.y = layer_scales(p[[i]])$y$range$range[2], size = 5)
      # p[[i]]$layers[[2]]$aes_params$textsize <- 5   
      
      #change plot type to boxplot or jitter plot based on user selection
      if (input$multi_readcountplot_type == 1){
        #boxplot
        p[[i]] <- ggboxplot(d, x = "condition", y = "count", color = "condition", add = "jitter") +
          ggtitle(val) +
          xlab("") +
          ylab("Normalized read count") +
          theme_classic() +
          theme(axis.text.x = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
                axis.text.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),
                axis.title.x = element_text(color="black",size=0),
                axis.title.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
                plot.title = element_text(color="black",size=input$multi_genecountFontSize_plot_title,angle=0,hjust=.5,vjust=.5,face="plain"),
                legend.title = element_blank(),
                legend.text =  element_text(color="black",size=input$multi_genecountFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
                legend.text.align = 0)
        
        #display statistics on plots if selected
        # if (input$multi_genecountStats == TRUE){
        #   #add statistics
        #   p[[i]] <- p[[i]] + stat_compare_means(label = "p.signif", method = "t.test", paired = TRUE, ref.group = input$control_condslist, label.y = layer_scales(p[[i]])$y$range$range[2]+input$multi_genecountStatsYcord, size = 6)
        #   p[[i]] <- p[[i]] + stat_compare_means(aes(label = paste0("p =", ..p.format..)), method = "t.test", label.y = layer_scales(p[[i]])$y$range$range[2]+input$multi_genecountStatsYcord, size = 6)
        #   p[[i]]$layers[[2]]$aes_params$textsize <- 6
        # }
         
        
      } else {
        #jitter plot
        p[[i]] <- ggplot(d, aes(x=condition, y=count, color=condition)) +
          ggtitle(val) +
          #geom_point(size = input$genecountPointSize) +
          geom_jitter(aes(size=input$multi_genecountPointSize)) +
          xlab("") +
          ylab("Normalized read count") +
          theme_classic() +
          theme(axis.text.x = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=.5,vjust=.5,face="plain"),
                axis.text.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=0,hjust=1,vjust=0,face="plain"),
                axis.title.x = element_text(color="black",size=0),
                axis.title.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"),
                plot.title = element_text(color="black",size=input$multi_genecountFontSize_plot_title,angle=0,hjust=.5,vjust=.5,face="plain"),
                legend.title = element_blank(),
                legend.text =  element_text(color="black",size=input$multi_genecountFontSize_legend_text,angle=0,hjust=.5,vjust=.5,face="plain"),
                legend.text.align = 0)
          
          #display statistics on plots if selected
          # if (input$multi_genecountStats == TRUE){
          #   #add statistics
          #   p[[i]] <- p[[i]] + stat_compare_means(label = "p.signif", method = "t.test", paired = TRUE, ref.group = input$control_condslist, label.y = layer_scales(p[[i]])$y$range$range[2]+input$multi_genecountStatsYcord, size = 6)
          #   p[[i]] <- p[[i]] + stat_compare_means(aes(label = paste0("p =", ..p.format..)), method = "t.test", label.y = layer_scales(p[[i]])$y$range$range[2]+input$multi_genecountStatsYcord, size = 6)
          # }
        
      }
      
      #show labels for points as determined by user
      if (input$multi_readcountplot_labels == 1){
        #no labels
        p[[i]] <- p[[i]]
      } else if (input$multi_readcountplot_labels == 2){
        #sample names as labels
        p[[i]] <- p[[i]] + geom_text_repel(size=input$multi_genecountLabelFontSize, 
                                           nudge_x=0.1, nudge_y=0.1, 
                                           segment.color=NA, aes(label=rownames(d))) 
          # + labs(shape="Replicate", color="Condition")
          #aes(shape=rownames(d))
      } else if (input$multi_readcountplot_labels == 3){
        #replicate names as labels
        p[[i]] <- p[[i]] + geom_text_repel(size=input$multi_genecountLabelFontSize, 
                                           nudge_x=0.1, nudge_y=0.1, 
                                           segment.color=NA, aes(label=replicate))
          # + labs(shape="Replicate", color="Condition")
      }
      # 
      # ###### this code is redundant now ######
      # #adjust legend based on user selection
      # if (input$multi_genecountShowLegends == 1){ # hide legend on all plots
      #   p[[i]] <- p[[i]] + theme(legend.position = "none")
      # } else if (input$multi_genecountShowLegends == 2 | input$multi_genecountShowLegends == 3){ # show legend on all plots or show one common legend
      #   p[[i]] <- p[[i]] + theme(legend.position = paste(input$multi_genecountLegendPosition))
      # } else if (input$multi_genecountShowLegends == 3){ # show one common legend
      #   
      # }
      # 
      #turn off y-axis unless the plot is the first one in a row, based on user selection
      if (input$multi_genecountSharedYAxis == TRUE){
        if ((i-1) %% input$multi_genecountGridColumns == 0 | i == 1){
          p[[i]] <- p[[i]] + theme(axis.title.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"))
        } else {
          p[[i]] <- p[[i]] + theme(axis.title.y = element_blank())
        }
      } else {
        p[[i]] <- p[[i]] + theme(axis.title.y = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=90,hjust=.5,vjust=.5,face="plain"))
      }
      
      #set y-axis to log scale depending on user selection
      if (input$multi_log10scale == TRUE){
        p[[i]] <- p[[i]] + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                         labels = trans_format("log10", math_format(10^.x)))
      }
      
      #rotate x-axis text based on user selection
      if (input$pcaRotateText == TRUE){
        p[[i]] <- p[[i]] + theme(axis.text.x = element_text(color="black",size=input$multi_genecountFontSize_xy_axis,angle=45,hjust=.5,vjust=.5,face="plain"))
      }

    }
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    #specify legend position
    if (input$multi_genecountShowLegends == 1){
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = FALSE, legend="none"))
    } else if (input$multi_genecountShowLegends == 3){
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = TRUE, legend=paste(input$multi_genecountLegendPosition)))
    } else {
      do.call(ggarrange, c(p, nrow=input$multi_genecountGridRows, ncol=input$multi_genecountGridColumns, common.legend = FALSE, legend=paste(input$multi_genecountLegendPosition)))
    }
    
  })
  
  get_unique_conds <- reactive({
    
    #use a random plot table sample
    temp <- plotCounts(dds, gene=listofgenes[1,2], intgroup=c("condition", "replicate"), returnData=TRUE)
    #determine number of rows
    nrow(temp)
    
    #empty vector to hold unique conditions
    newlist <- NULL
    
    newtemp <- as.matrix(temp)
    
    newlist[1] <- newtemp[1,2]
    
    counter = 2
    for (i in 2:nrow(newtemp)){
      
      prev = newtemp[i-1,2]
      curr = newtemp[i,2]
      
      if (curr != prev){
        newlist[counter] <- curr
        counter = counter+1
      }
    }
    
    return(newlist)
    
  })  
  
  dynamic_colorgen <- function(widget_name, selector)({
    
    #get list of unique condition names
    unique_conds <- get_unique_conds()
    
    #initialize vector to hold dynamically generated names for color widgets
    colorwidget_names <- NULL
    
    for (count in 1:num_conditions){
      
      #prepare unique identifier name for current widget in the loop
      colorwidget_names[count] <- paste0(widget_name, count)
      
      insertUI(
        selector = selector,
        ui = colourInput(colorwidget_names[count], unique_conds[count], "#EB4141", allowTransparent = FALSE)
      )
    }
    
  })
  
})


