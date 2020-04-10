# BEAVR: A Browser-based tool for the Exploration And Visualization of RNA-seq data
# GUI to analyze RNAseq data using DESeq2
# input: transcript read counts (ie. from STAR aligner or HTseq), and column data matrix file containing sample info
# See Github for more info & ReadMe: https://github.com/developerpiru/BEAVR

app_version = "1.0.7"

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
# installReqs("BiocManager", bioc = FALSE)
# installReqs("shiny", bioc = FALSE)
# installReqs("shinydashboard", bioc = FALSE)
# installReqs("plotly", bioc = FALSE)
# installReqs("ggplot2", bioc = FALSE)
# installReqs("ggrepel", bioc = FALSE)
# installReqs("data.table", bioc = FALSE)
# installReqs("DT", bioc = FALSE)
# installReqs("DESeq2", bioc = TRUE)
# installReqs("vsn", bioc = TRUE)
# installReqs('apeglm', bioc = TRUE)
# installReqs('org.Hs.eg.db', bioc = TRUE)
# installReqs('org.Mm.eg.db', bioc = TRUE)
# installReqs('EnhancedVolcano', bioc = TRUE)
# installReqs('gridExtra', bioc = FALSE)
# installReqs('ggpubr', bioc = FALSE)
# installReqs('shinyjqui', bioc = FALSE)
# installReqs('scales', bioc = FALSE)
# installReqs('RColorBrewer', bioc = FALSE)
# installReqs('pheatmap', bioc = FALSE)
# installReqs('colourpicker', bioc = FALSE)

#load required libraries
#install_github("jokergoo/ComplexHeatmap") # install ComplexHeatmap package
# require("BiocManager")
# require("colourpicker")
# require("ComplexHeatmap")
# require("data.table")
# require("DESeq2")
# require("devtools") # to install from github
# require("DT")
# require("ggplot2")
# require("ggpubr")
# require("ggrepel")
# require("gridExtra")
# require("pheatmap")
# require("RColorBrewer")
# require("scales")
# require("shiny")
# require("shinydashboard")
# require("shinyjqui")
# require("shinyWidgets")
# require("vsn")
# require("apeglm")
# require("circlize")
# require("EnhancedVolcano")
# require("enrichplot") #from bioconductor
# require("ggraph") #from cran, for enrichment maps
# require("org.Hs.eg.db")
# require("org.Mm.eg.db")
# require("ReactomePA") #from bioconductor
# require("shinycssloaders")

#set port to 3838
options(shiny.port = 3838)

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
    
    selectInput("control_condslist", "Choose control condition", temp_condslist, selected = temp_condslist[1])
  })
  
  #get treatment condition from list
  output$treatment1_condslist <- renderUI({
    temp_coldata <- coldata()
    temp_condslist <- unique(temp_coldata[,2])
    
    selectInput("treatment1_condslist", "Choose treatment condition", temp_condslist, selected = temp_condslist[2])
  })
  
  #get false discovery rate from user
  output$FDR_value <- renderUI({
    numericInput("FDRvalue", "False Discovery Rate %",value = 10)
  })
  
  #get minimum read count values to keep from user
  output$min_reads <- renderUI({
    numericInput("min_reads_value", "Drop genes with reads below:",value = 10)
  })
  
  #function to dynamically update the dropdown box for detected GSEA pathways
  output$gseaPlotPathways <- renderUI({
    #gsea_plot_data[paste0(input$gseaPlotPathway), "Description"]
    
    selectInput("gseaPlotPathways", "Detected pathways (from most to least significant)", 
                gsea_plot_data$Description, 
                selected = gsea_plot_data[1, "Description"])
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
    num_replicates <<- length(unique(temp_coldata[,3]))
    
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
    first_run_flag3 <<- TRUE
    
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
        filteredTable <<- unfilteredTable
        
        #reset all values to max values
        #update numericInputs for data filtering with the min/max values from resFDR
        updateNumericInput(session, "log2FC_min", value = min(filteredTable$log2FoldChange, na.rm=T))
        updateNumericInput(session, "log2FC_max", value = max(filteredTable$log2FoldChange, na.rm=T))
        updateNumericInput(session, "pvalue_min", value = min(filteredTable$pvalue, na.rm=T))
        updateNumericInput(session, "pvalue_max", value = max(filteredTable$pvalue, na.rm=T))
        updateNumericInput(session, "padj_min", value = min(filteredTable$padj, na.rm=T))
        updateNumericInput(session, "padj_max", value = max(filteredTable$padj, na.rm=T))
        updateNumericInput(session, "baseMean_min", label = "Min", value = min(filteredTable$baseMean, na.rm=T))
        updateNumericInput(session, "baseMean_max", label = "Max", value = max(filteredTable$baseMean, na.rm=T))
        updateNumericInput(session, "lfcSE_min", label = "Min", value = min(filteredTable$lfcSE, na.rm=T))
        updateNumericInput(session, "lfcSE_max", label = "Max", value = max(filteredTable$lfcSE, na.rm=T))
        updateNumericInput(session, "stat_min", value = min(filteredTable$stat, na.rm=T))
        updateNumericInput(session, "stat_max", value = max(filteredTable$stat, na.rm=T))
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
      #write to file
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
    vsd <<- varianceStabilizingTransformation(dds, blind = FALSE)

    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Plotting..."))

    #run only if first_run_flag2 boolean is TRUE; meaning ui is being initialized
    if (first_run_flag2 == TRUE){
      #call function to loop through and create more colour widgets for each condition in the experiment
      #widget names are pcaColorX, where X is an integer
      #selector is the target div tag container in ui
      dynamic_colorgen(widget_name = "pcaColor", selector = "#pcaColorbox", sourcelist = "condition")
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
      do_sampleClustering_plot2()
    })
  })
  
  #function to calculate sample clustering using ComplexHeatmap
  do_sampleClustering_plot2 <- reactive({
    
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
    vsd <<- varianceStabilizingTransformation(dds, blind = FALSE)
    
    #transpose the data
    sampleDists <- dist(t(assay(vsd)))
    
    #generate heatmap for sample clustering
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
    colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
    map_colors <- colorRampPalette( rev(brewer.pal(9, paste(input$sampleClustering_mapColor))) )(255)
    
    
    #generate heatmap
    p = Heatmap(sampleDistMatrix,
                #name
                name = " ",
                #colour
                col = map_colors,
                #row distance method
                clustering_distance_rows = input$sampleClustering_method,
                #column distance method
                clustering_distance_columns = input$sampleClustering_method,
                #position of gene names
                row_names_side = input$sampleClusterin_rowlabel_position,
                #position of sample names
                column_names_side = input$sampleClusterin_collabel_position,
                #rotate row names
                row_names_rot = input$sampleClusterin_row_rotation,
                #rotate column names
                column_names_rot = input$sampleClusterin_col_rotation,
                #position of row dendrogam
                row_dend_side = input$sampleClusterin_row_dend_position, 
                #position of column dendrogram
                column_dend_side = input$sampleClusterin_col_dend_position,
                #width of row dendrogram
                row_dend_width = unit(input$sampleClusterin_row_dend_width, "cm"),
                #height of column dendrogram
                column_dend_height = unit(input$sampleClusterin_col_dend_height, "cm"),
                #cell border settings
                rect_gp = gpar(col = input$sampleClustering_borderColor, lwd = 1),
                #size and color of row names
                row_names_gp = gpar(fontsize = input$sampleClustering_fontsize_rowNames, col = input$sampleClustering_rowlabelColor),
                #size and color of column names
                column_names_gp = gpar(fontsize = input$sampleClustering_fontsize_colNames, col = input$sampleClustering_collabelColor),
                #legend direction
                heatmap_legend_param = list(direction = input$sampleClustering_main_legend_dir,
                                            labels_gp = gpar(fontsize = input$sampleClustering_fontsize_legends, col = input$sampleClustering_legendColor)),
                #show cell values
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(input$sampleClustering_cellNums == "%.2f"){
                    grid.text(sprintf("%.1f", sampleDistMatrix[i, j]), x, y, gp = gpar(fontsize = input$sampleClustering_fontsize_cellNums, col = input$sampleClustering_cellNumsColor))
                  } else if (input$sampleClustering_cellNums == "%.1e"){
                    grid.text(sprintf("%.1e", sampleDistMatrix[i, j]), x, y, gp = gpar(fontsize = input$sampleClustering_fontsize_cellNums, col = input$sampleClustering_cellNumsColor))
                  }
                }
                
    )

    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #show the heatmap
    draw(p, padding = unit(c(10, 10, 10, 10), "mm"),
         merge_legend = TRUE, 
         heatmap_legend_side = input$sampleClustering_main_legend)
    
  })
  
  #call function to show count matrix heatmap
  output$countMatrix_heatmap = renderPlot({
    withProgress(message = 'Generating heatmap...', value = 1, min = 1, max = 100, {
      do_countMatrix_heatmap2()
    })
  })
  
  #new function to draw count matrix heatmap using ComplexHeatmap
  do_countMatrix_heatmap2 <- reactive({
    
    #Update progress bar
    totalSteps = 3 + 3
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get some user defined values
    #clust_dist = input$heatmap_distance
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
    if (input$heatmap_varlogmethod == "vst")
      transform_data <<- vst(dds, blind=FALSE)
    else
      transform_data <<- rlog(dds, blind=FALSE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Subsetting..."))
    
    #get annotations for labels
    annot <- as.data.frame(colData(dds)[,c("condition","replicate")])
    
    if (input$heatmap_pickTopGenes == TRUE){
      #get the top X genes as defined by user
      #genestokeep <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:input$heatmap_numGenes]
      genestokeep <- order(rowVars(assay(transform_data)), decreasing = TRUE)[1:input$heatmap_numGenes]
      
      #old code from version 0.73.2
      #genestokeep <- head(order(rowVars(assay(transform_data)), decreasing=TRUE), 50)
      #heatmap_data <- as.data.frame(assay(transform_data)[genestokeep,])
      
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
    
    #subset the required genes from the transformed data   
    heatmap_data <- as.data.frame(assay(transform_data)[genestokeep,])
    
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
    
    #scale the data across all genes
    heatmap_data2 = t(apply(heatmap_data, 1, function(x) {
      scale(x)
    }))
    
    #set sample names (column names)
    colnames(heatmap_data2) = colnames(heatmap_data)
    
    #run only if first_run_flag3 boolean is TRUE; meaning ui is being initialized
    if (first_run_flag3 == TRUE){
      #call function to loop through and create more colour widgets for each replicate and condition in the experiment
      #widget names are heatmap_replicateColorX and heatmap_conditionColorX, where X is an integer
      #selector is the target div tag container in ui
      dynamic_colorgen(widget_name = "heatmap_replicateColor", selector = "#heatmap_replicateColorbox", sourcelist = "replicate")
      dynamic_colorgen(widget_name = "heatmap_conditionColor", selector = "#heatmap_conditionColorbox", sourcelist = "condition")
    }
    
    #set first run flag to false so color widgets are no longer made
    first_run_flag3 <<- FALSE
    
    #vector to save colours
    replicate_colors = 0
    condition_colors = 0

    #loop through and get the colours that the user chooses for replicate annotations
    for (count in 1:num_replicates){
      replicate_colors[count] = input[[paste0('heatmap_replicateColor', count)]]
    }
    
    #loop through and get the colours that the user chooses for condition annotations
    for (count in 1:num_conditions){
      condition_colors[count] = input[[paste0('heatmap_conditionColor', count)]]
    }
    
    #set the names of colour vectors to replicate names or condition names, respectively
    names(replicate_colors) = unique(temp_coldata[,3]) # column 3 is replicates
    names(condition_colors) = unique(temp_coldata[,2]) # column 2 is conditions

    if (input$heatmap_anno_legend_dir == "vertical")
      anno_horizontal_flip = num_replicates
    else if (input$heatmap_anno_legend_dir == "horizontal")
      anno_horizontal_flip = 1
    
    #define the replicate annotation
    if (input$heatmap_annotations == "replicate"){
      heatmap_anno = HeatmapAnnotation(Replicates = as.matrix(colData(dds)[,c("replicate")]), 
                                       col = list(Replicates = replicate_colors),
                                       annotation_name_gp = gpar(fontsize = input$heatmap_fontsize_annotations),
                                       annotation_legend_param = list(Replicates = list(nrow = anno_horizontal_flip, 
                                                                                        title_gp = gpar(fontsize = input$heatmap_fontsize_legends, fontface = "bold"), 
                                                                                        labels_gp = gpar(fontsize = input$heatmap_fontsize_legends))))
    } else if (input$heatmap_annotations == "treatment") {
      heatmap_anno = HeatmapAnnotation(Condition = as.matrix(colData(dds)[,c("condition")]), 
                                       col = list(Condition = condition_colors),
                                       annotation_name_gp = gpar(fontsize = input$heatmap_fontsize_annotations),
                                       annotation_legend_param = list(Condition = list(nrow = anno_horizontal_flip, 
                                                                                       title_gp = gpar(fontsize = input$heatmap_fontsize_legends, fontface = "bold"), 
                                                                                       labels_gp = gpar(fontsize = input$heatmap_fontsize_legends))))
    } else if (input$heatmap_annotations == "both") {
      heatmap_anno = HeatmapAnnotation(Replicates = as.matrix(colData(dds)[,c("replicate")]), 
                                       Condition = as.matrix(colData(dds)[,c("condition")]), 
                                       col = list(Replicates = replicate_colors, Condition = condition_colors),
                                       annotation_name_gp = gpar(fontsize = input$heatmap_fontsize_annotations),
                                       annotation_legend_param = list(Replicates = list(nrow = anno_horizontal_flip, 
                                                                                        title_gp = gpar(fontsize = input$heatmap_fontsize_legends, fontface = "bold"), 
                                                                                        labels_gp = gpar(fontsize = input$heatmap_fontsize_legends)),
                                                                      Condition = list(nrow = anno_horizontal_flip, 
                                                                                       title_gp = gpar(fontsize = input$heatmap_fontsize_legends, fontface = "bold"), 
                                                                                       labels_gp = gpar(fontsize = input$heatmap_fontsize_legends))))
    } else if (input$heatmap_annotations == "none") {  
      heatmap_anno = NULL
    }
    
    #generate heatmap
    p = Heatmap(heatmap_data2,
                #name
                name = "Expression",
                #sample annotation,
                top_annotation = heatmap_anno,
                #colour
                col = colorRamp2(c(input$heatmap_scale_range[1], 0, input$heatmap_scale_range[2]), c(input$heatmap_lowColor, input$heatmap_midColor, input$heatmap_highColor)),
                #cluster rows?
                cluster_rows = input$heatmap_clustRows,
                #cluster columns?
                cluster_columns = input$heatmap_clustCols,
                #row distance method
                clustering_distance_rows = input$heatmap_distance,
                #column distance method
                clustering_distance_columns = input$heatmap_distance,
                #row cluster method
                clustering_method_rows = input$heatmap_clustMethod,
                #column cluster method
                clustering_method_columns = input$heatmap_clustMethod,
                #position of gene names
                row_names_side = input$heatmap_genelabel_position,
                #position of sample names
                column_names_side = input$heatmap_samplelabel_position,
                #rotate row names
                row_names_rot = input$heatmap_row_rotation,
                #rotate column names
                column_names_rot = input$heatmap_col_rotation,
                #position of row dendrogam
                row_dend_side = input$heatmap_row_dend_position, 
                #position of column dendrogram
                column_dend_side = input$heatmap_col_dend_position,
                #width of row dendrogram
                row_dend_width = unit(input$heatmap_row_dend_width, "cm"),
                #height of column dendrogram
                column_dend_height = unit(input$heatmap_col_dend_height, "cm"),
                #cell border settings
                rect_gp = gpar(col = input$heatmap_borderColor, lwd = 1),
                #size and color of row names
                row_names_gp = gpar(fontsize = input$heatmap_fontsize_geneNames, col = input$heatmap_rowlabelColor),
                #size and color of column names
                column_names_gp = gpar(fontsize = input$heatmap_fontsize_sampleNames, col = input$heatmap_collabelColor),
                #legend direction, size, color
                heatmap_legend_param = list(direction = input$heatmap_main_legend_dir,
                                            title_gp = gpar(fontsize = input$heatmap_fontsize_legends, fontface = "bold"),
                                            labels_gp = gpar(fontsize = input$heatmap_fontsize_legends, 
                                                             col = input$heatmap_legendColor)),
                #show cell values
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(input$heatmap_cellNums == "%.2f"){
                    grid.text(sprintf("%.1f", heatmap_data2[i, j]), x, y, gp = gpar(fontsize = input$heatmap_fontsize_cellNums, col = input$heatmap_cellNumsColor))
                  } else if (input$heatmap_cellNums == "%.1e"){
                    grid.text(sprintf("%.1e", heatmap_data2[i, j]), x, y, gp = gpar(fontsize = input$heatmap_fontsize_cellNums, col = input$heatmap_cellNumsColor))
                  }
                }
                
                )
    
    # #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #show the heatmap
    draw(p, padding = unit(c(10, 10, 10, 10), "mm"),
         #merge_legend = TRUE, 
         heatmap_legend_side = input$heatmap_main_legend, 
         annotation_legend_side = input$heatmap_anno_legend)
    
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
      RNAseqdatatoplot <- as.data.frame(filteredTable)
    } else {
      #otherwise, use unfiltered data
      RNAseqdatatoplot <- as.data.frame(unfilteredTable)
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
  
    
    #run only if first_run_flag1 boolean is TRUE; meaning ui is being initialized
    if (first_run_flag1 == TRUE){
      #call function to loop through and create more colour widgets for each condition in the experiment
      #widget names are pcaColorX, where X is an integer
      #selector is the target div tag container in ui
      dynamic_colorgen(widget_name = "multi_genecountColor", selector = "#multi_genecountColorbox", sourcelist = "condition")
    }

    #set first run flag to false so color widgets are no longer made
    first_run_flag1 <<- FALSE
    
    #vector to save colours
    multi_colorslist <- NULL

    #loop through and get the colours that the user choses
    for (count in 1:num_conditions){
      multi_colorslist[count] <- input[[paste0('multi_genecountColor', count)]]
    }
    
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
      
      #get read counts
      d <<- plotCounts(dds, gene=listofgenes[which(listofgenes$GeneID==val),2], intgroup=c("condition", "replicate"), returnData=TRUE)

      #change plot type to boxplot or jitter plot based on user selection
      if (input$multi_readcountplot_type == 1){
        #boxplot
        p[[i]] <- ggboxplot(d, x = "condition", y = "count", fill = "condition", palette = multi_colorslist) +
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

      } else {
        #jitter plot
        p[[i]] <- ggplot(d, aes(x = condition, y = count)) +
          ggtitle(val) +
          geom_jitter(aes(color = condition), size = input$multi_genecountPointSize) +
          scale_color_manual(values = multi_colorslist) + # color palette for jitter plot
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
      } else if (input$multi_readcountplot_labels == 3){
        #replicate names as labels
        p[[i]] <- p[[i]] + geom_text_repel(size=input$multi_genecountLabelFontSize, 
                                           nudge_x=0.1, nudge_y=0.1, 
                                           segment.color=NA, aes(label=replicate))
      }

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
      
      #display statistics on plots if selected
      if (input$multi_genecountStats == TRUE){
        #vector to save colours
        stat_y_locations <- NULL
        
        #for loop to cycle through conditions and determine max values in each so stats can be plotted above it
        for (cond.i in length(unique(d$condition))){
          #determine max
          if (input$multi_log10scale == TRUE){
            stat_y_locations[cond.i] <- log10(max(subset(d, condition == unique(d$condition)[cond.i])[,1]) * (1 + input$multi_genecountStatsYcord/100))
          } else {
            stat_y_locations[cond.i] <- max(subset(d, condition == unique(d$condition)[cond.i])[,1]) * (1 + input$multi_genecountStatsYcord/100)
          }
        }
        
        #add statistics, asterisks
        p[[i]] <- p[[i]] + stat_compare_means(label = "p.signif", method = "t.test", paired = FALSE, ref.group = input$control_condslist, label.y = stat_y_locations, size = input$multi_genecountFontSize_stats_text)
        #add p value
        #p[[i]] <- p[[i]] + stat_compare_means(aes(label = paste0("p =", ..p.format..)), method = "t.test", label.y = layer_scales(p[[i]])$y$range$range[2]+input$multi_genecountStatsYcord, size = 6)
        #p[[i]]$layers[[2]]$aes_params$textsize <- input$multi_genecountFontSize_stats_text
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
  
  #generate enrichment pathway barplot or dotplot
  output$enrichmentPlot = renderPlot({
    withProgress(message = 'Generating pathway enrichment plot...', value = 1, min = 1, max = 100, {
      do_enrichment_plot()
    })
  })
  
  #function to plot enrichment pathway barplot or dotplot
  do_enrichment_plot <- reactive({
    
    #Update progress bar
    totalSteps = 3 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating enrichment results..."))
    
    #call function to calculate pathway enrichment results
    enrichment_data <- calc_enrichment_results()
    
    #convert to data frame
    enrichment_data2 <- as.data.frame(enrichment_data)
    colnames(enrichment_data2) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "padjust", "qvalue", "geneID", "Count")
    
    #sort and take only the top x (enrNumCategories) categories
    enrichment_data2 <- enrichment_data2[order(enrichment_data2$padjust),]
    enrichment_data2 <- enrichment_data2[1:input$enrNumCategories,]
    enrichment_data2 <- enrichment_data2[order(enrichment_data2$Count),]
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Generating enrichment plot..."))
    
    #check which type of plot is selected: bar plot or dot plot
    if (input$enrPlotType == "bar") {
      #draw bar plot
      p <- ggplot(enrichment_data2, aes(x = reorder(Description, -padjust), y = Count, fill = padjust)) +
        coord_flip() +
        xlab("") +
        ylab("Gene count") +
        geom_bar(stat="identity") +
        scale_fill_continuous(low = input$enrLeastColor, high = input$enrMostColor, trans = "reverse", labels = to_scientific) +
        theme_classic() +
        theme(axis.text.y = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.title.y = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.text.x = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.title.x = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              legend.position = input$enrLegendPosition,
              legend.direction = input$enrLegendDirection,
              legend.title = element_text(color="black",size=input$enrFontSize_legend,angle=0,hjust=0,vjust=.5,face="bold"),
              legend.text =  element_text(color="black",size=input$enrFontSize_legend,angle=0,hjust=0,vjust=.5,face="plain"),
              legend.text.align = 0)
      
      p <- p + guides(colour = "colorbar")
      p <- p + guides(fill = guide_colourbar(title = "P adjusted"))
    
    } else if (input$enrPlotType == "test"){
      
      p <- dotplot(enrichment_data, showCategory = input$enrNumCategories) +
        scale_color_gradient(low = "#56B1F7", high = "#132B43")
        
    } else {
      #draw dot plot
      p <- dotplot(enrichment_data, showCategory = input$enrNumCategories) +
        scale_color_continuous(low = input$enrLeastColor, high = input$enrMostColor, trans = "reverse", labels = to_scientific) +
        theme_classic() +
        theme(axis.text.y = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.title.y = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.text.x = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              axis.title.x = element_text(color="black",size=input$enrFontSize_xy_axis,angle=0,face="plain"),
              legend.position = input$enrLegendPosition,
              legend.direction = input$enrLegendDirection,
              legend.title = element_text(color="black",size=input$enrFontSize_legend,angle=0,hjust=0,vjust=.5,face="bold"),
              legend.text =  element_text(color="black",size=input$enrFontSize_legend,angle=0,hjust=0,vjust=.5,face="plain"),
              legend.text.align = 0)
      
      p <- p + guides(colour = "colorbar", size = "legend")
      p <- p + guides(size = guide_legend(title = "Gene count"))
      p <- p + guides(color = guide_colorbar(title = "P adjusted"))
      
    }

    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    print(p)
  })
  
  calc_enrichment_results <- reactive ({
    
    #check if filtering is enabled for the data table
    if (input$filterTableEnabled == FALSE){
      print("Warning: you have not enabled filtering of the data table!")
    }
    
    #Update progress bar
    totalSteps = 3 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get gene input list from the filtered gene data table
    #filtered data is stored in the global variable filteredTable
    input_data <- filteredTable
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Getting Entrez IDs..."))
    
    #convert from ENSEMBL to ENTREZID
    input_data$EntrezID <- mapIds(org.Hs.eg.db,keys=rownames(input_data),column="ENTREZID",keytype="ENSEMBL",multiVals="first")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating enrichment..."))
    
    #calculate pathway enrichment
    enrichment_data <- enrichPathway(gene = input_data$EntrezID, pvalueCutoff = input$enrPvalcutoff)#, readable=T)
    
    #return the data
    return(enrichment_data)
  })
  
  #generate enrichment pathway map
  output$enrichmentMap = renderPlot({
    withProgress(message = 'Generating pathway enrichment plot...', value = 1, min = 1, max = 100, {
      do_enrichment_map()
    })
  })
  
  #function to generate enrichment pathway map
  do_enrichment_map <- reactive({
    
    #check if filtering is enabled for the data table
    if (input$filterTableEnabled == FALSE){
      print("Warning: you have not enabled filtering of the data table!")
    }
    
    #Update progress bar
    totalSteps = 6 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get gene input list from the filtered gene data table
    #filtered data is stored in the global variable filteredTable
    input_data <- filteredTable
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Getting Entrez IDs..."))
    
    #convert from ENSEMBL to ENTREZID
    input_data$EntrezID <- mapIds(org.Hs.eg.db,keys=rownames(input_data),column="ENTREZID",keytype="ENSEMBL",multiVals="first")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating enrichment..."))
    
    #calculate pathway enrichment
    enrichment_data <- enrichPathway(gene = input_data$EntrezID, pvalueCutoff = input$enrMapPvalcutoff, readable = TRUE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Generating enrichment map..."))
    
    #generate enrichment map
    p <- emapplot(enrichment_data) +
      scale_color_continuous(low = input$enrMapLeastColor, high = input$enrMapMostColor, trans = "reverse", labels = to_scientific) +
      theme_void() +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = input$enrMapLegendPosition,
            legend.direction = input$enrMapLegendDirection,
            legend.title = element_text(color="black",size=input$enrMapFontSize_legend,angle=0,hjust=0,vjust=.5,face="bold"),
            legend.text =  element_text(color="black",size=input$enrMapFontSize_legend,angle=0,hjust=0,vjust=.5,face="plain"),
            legend.text.align = 0)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Adjusting aesthetics..."))
    
    #modify legend titles
    p <- p + guides(colour = "colorbar", size = "legend")
    p <- p + guides(size = guide_legend(title = "Gene count"))
    p <- p + guides(color = guide_colorbar(title = "P adjusted"))
    
    #set layer 3 in p to null because it contains the default geom_node_text labels
    p$layers[[3]] <- NULL
    
    #create a new geom_node_text layer with customizable node labels
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE, size = input$enrMapFontSize_labels)
    
    #create a new layer with new geom_edge_links to customize edge/line colors (layer 4)
    p <- p + geom_edge_link(alpha=(100-input$enrMaptransparency)/100, aes_(width=~I(width)), colour=input$enrMaplineColor)
    
    #set layer 1 to equal layer 4 so the new edges are drawn first, then set layer 4 to NULL
    p$layers[[1]] <- p$layers[[4]]
    p$layers[[4]] <- NULL
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    print(p)
  })
  
  #generate GSEA map
  output$gseaMap = renderPlot({
    withProgress(message = 'Generating pathway enrichment plot...', value = 1, min = 1, max = 100, {
      do_gsea_map()
    })
  })
  
  #function to generate GSEA map
  do_gsea_map <- reactive({
    
    #check if filtering is enabled for the data table
    if (input$filterTableEnabled == FALSE){
      print("Warning: you have not enabled filtering of the data table!")
    }
    
    #Update progress bar
    totalSteps = 4 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #calculate function to GSEA
    enrichment_data <- calc_gsea()
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Generating enrichment map..."))
    
    #generate enrichment map
    p <- emapplot(enrichment_data, color="pvalue") +
      scale_color_continuous(low = input$gseaMapLeastColor, high = input$gseaMapMostColor, trans = "reverse", labels = to_scientific) +
      theme_void() +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = input$gseaMapLegendPosition,
            legend.direction = input$gseaMapLegendDirection,
            legend.title = element_text(color="black",size=input$gseaMapFontSize_legend,angle=0,hjust=0,vjust=.5,face="bold"),
            legend.text =  element_text(color="black",size=input$gseaMapFontSize_legend,angle=0,hjust=0,vjust=.5,face="plain"),
            legend.text.align = 0)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Adjusting aesthetics..."))
    
    #modify legend titles
    p <- p + guides(colour = "colorbar", size = "legend")
    p <- p + guides(size = guide_legend(title = "Gene count"))
    p <- p + guides(color = guide_colorbar(title = "P adjusted"))
    
    #set layer 3 in p to null because it contains the default geom_node_text labels
    p$layers[[3]] <- NULL
    
    #create a new geom_node_text layer with customizable node labels
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE, size = input$gseaMapFontSize_labels)
    
    #create a new layer with new geom_edge_links to customize edge/line colors (layer 4)
    p <- p + geom_edge_link(alpha=(100-input$gseaMaptransparency)/100, aes_(width=~I(width)), colour=input$gseaMaplineColor)
    
    #set layer 1 to equal layer 4 so the new edges are drawn first, then set layer 4 to NULL
    p$layers[[1]] <- p$layers[[4]]
    p$layers[[4]] <- NULL
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    print(p)
  })
  
  #function to calculate gsea results
  calc_gsea <- reactive({
    
    #check if filtering is enabled for the data table
    if (input$filterTableEnabled == FALSE){
      print("Warning: you have not enabled filtering of the data table!")
    }
    
    #Update progress bar
    totalSteps = 3 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get gene input list from the filtered gene data table
    #filtered data is stored in the global variable filteredTable
    input_data <- filteredTable
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Getting Entrez IDs..."))
    
    #convert from SYMBOL to ENTREZID
    input_data$EntrezID <- mapIds(org.Hs.eg.db,keys=input_data$GeneID,column="ENTREZID",keytype="SYMBOL",multiVals="first")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating GSEA results..."))
    
    #Format data into geneList compatible format
    #get the LFC values
    geneList = input_data$log2FoldChange
    
    #get the Entrez IDs
    names(geneList) = as.character(input_data$EntrezID)
    
    #sort by decreasing order
    geneList = sort(geneList, decreasing = TRUE)
    
    #calculate GSEA
    enrichment_data <- gsePathway(geneList, nPerm = 10000, 
                                  pvalueCutoff = input$gseaMapPvalcutoff, pAdjustMethod = "BH", verbose = FALSE)
    
    #return data
    return(enrichment_data)
  })  
  
  #generate GSEA plot for a single pathway
  output$gseaPlot = renderPlot({
    withProgress(message = 'Generating pathway enrichment plot...', value = 1, min = 1, max = 100, {
      do_gsea_plot()
    })
  })
  
  #function to update p values below GSEA plot
  output$gseaPlotpvals = renderText({
    
    if (is.null(input$gseaPlotPathways) == FALSE){
      if (is.null(gsea_plot_data) == FALSE){
        #pull p value and p adjusted values from gsea_plot_data once calculated
        paste("p value: ", format(gsea_plot_data[which(gsea_plot_data$Description == input$gseaPlotPathways),"pvalue"], digits = 3, scientific = TRUE),
              " and p adjusted: ", format(gsea_plot_data[which(gsea_plot_data$Description == input$gseaPlotPathways),"p.adjust"], digits = 3, scientific = TRUE))
        }
    }
    
  })

  #function to generate GSEA plot for a single pathway
  do_gsea_plot <- reactive({
    
    #check if filtering is enabled for the data table
    if (input$filterTableEnabled == FALSE){
      print("Warning: you have not enabled filtering of the data table!")
    }
    
    #Update progress bar
    totalSteps = 5 + 2
    currentStep = 1
    incProgress(currentStep/totalSteps*100, detail = paste("Initializing..."))
    
    #get gene input list from the filtered gene data table
    #filtered data is stored in the global variable filteredTable
    input_data <- filteredTable
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Getting Entrez IDs..."))
    
    #convert from SYMBOL to ENTREZID
    input_data$EntrezID <- mapIds(org.Hs.eg.db,keys=input_data$GeneID,
                                  column="ENTREZID",keytype="SYMBOL",multiVals="first")
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Calculating GSEA..."))
    
    #Format data into geneList compatible format
    #get the LFC values
    geneList = input_data$log2FoldChange
    
    #get the Entrez IDs
    names(geneList) = as.character(input_data$EntrezID)
    
    #sort by decreasing order
    geneList = sort(geneList, decreasing = TRUE)
    
    #calculate GSEA
    gsea_plot_data <<- gsePathway(geneList, nPerm = 10000,
                                  pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Generating enrichment map..."))
    
    #set the name of the pathway to plot
    if (is.null(input$gseaPlotPathways)){
      #if input$gseaPlotPathways is NULL, then set pathway_id to default value
      pathway_id = gsea_plot_data[1, "ID"]
    } else {
      #otherwise use determine the Reactome ID for the pathway selected in input$gseaPlotPathways
      pathway_id = gsea_plot_data[which(gsea_plot_data$Description == input$gseaPlotPathways),"ID"]
    }
    
    p0 <<- gseaplot2(gsea_plot_data, geneSetID = pathway_id,
                     title = gsea_plot_data[paste0(pathway_id), "Description"],
                     pvalue_table = input$gseaPlotShowPvalue,
                     base_size = input$gseaPlotFontSize,
                     color = input$gseaPlotLineColor)
    
    #Update progress bar
    currentStep = currentStep + 1
    incProgress(currentStep/totalSteps*100, detail = paste("Finalizing..."))
    
    #return the plot
    print(p0)
  })
  
  #function to show GSEA results table
  output$show_gseaTable <- DT::renderDataTable({
    
    withProgress(message = 'Performing calculations...', value = 1, min = 1, max = 100, {
      gsea_results <- as.data.frame(calc_gsea())
      rownames(gsea_results) <- NULL
      
      #rename columns
      colnames(gsea_results) <- c("ID", "Description", "setSize", "ES", "NES", "pvalue", "padj", "qvalue", "rank", "leadEdge", "coreEnrichment")
      
      if (input$gseaTablefilterTableEnabled == TRUE) {
        #filter the table based on user defined settings
        filtered_gsea_results <<- subset(gsea_results, 
                                 ES >= input$gsea_enrichmentScore_min &
                                  ES <= input$gsea_enrichmentScore_max &
                                  NES >= input$gsea_nes_min &
                                  NES <= input$gsea_nes_max &
                                  pvalue >= input$gsea_pvalue_min &
                                  pvalue <= input$gsea_pvalue_max &
                                  padj >= input$gsea_padj_min &
                                  padj <= input$gsea_padj_max &
                                  qvalue >= input$gsea_qvalue_min &
                                  qvalue <= input$gsea_qvalue_max &
                                  rank >= input$gsea_rank_min &
                                  rank <= input$gsea_rank_max
                                )
        
      } else {
        #no filtering enabled, show all values and update filtering options
        filtered_gsea_results <<-gsea_results
        
        #update numericInputs for data filtering
        updateNumericInput(session, "gsea_enrichmentScore_min", value = min(gsea_results$ES, na.rm=T))
        updateNumericInput(session, "gsea_enrichmentScore_max", value = max(gsea_results$ES, na.rm=T))
        updateNumericInput(session, "gsea_nes_min", value = min(gsea_results$NES, na.rm=T))
        updateNumericInput(session, "gsea_nes_max", value = max(gsea_results$NES, na.rm=T))
        updateNumericInput(session, "gsea_pvalue_min", value = min(gsea_results$pvalue, na.rm=T))
        updateNumericInput(session, "gsea_pvalue_max", value = max(gsea_results$pvalue, na.rm=T))
        updateNumericInput(session, "gsea_padj_min", value = min(gsea_results$padj, na.rm=T))
        updateNumericInput(session, "gsea_padj_max", value = max(gsea_results$padj, na.rm=T))
        updateNumericInput(session, "gsea_qvalue_min", value = min(gsea_results$qvalue, na.rm=T))
        updateNumericInput(session, "gsea_qvalue_max", value = max(gsea_results$qvalue, na.rm=T))
        updateNumericInput(session, "gsea_rank_min", value = min(gsea_results$rank, na.rm=T))
        updateNumericInput(session, "gsea_rank_max", value = max(gsea_results$rank, na.rm=T))
        
      }
      
      #show table
      filtered_gsea_results
    })
    
  })
  
  #function to show Pathway Enrichment results table
  output$show_enrichmentTable <- DT::renderDataTable({
    
    withProgress(message = 'Performing calculations...', value = 1, min = 1, max = 100, {
      enrichment_results <- as.data.frame(calc_enrichment_results())
      rownames(enrichment_results) <- NULL
      
      #rename columns
      colnames(enrichment_results) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "padj", "qvalue", "GeneID", "Count")
      #rearrange table for display
      enrichment_results <- enrichment_results[,c(1:7,9,8)]
      
      if (input$enrichmentTablefilterTableEnabled == TRUE) {
        #filter the table based on user defined settings
        filtered_enrichment_results <<- subset(enrichment_results,
                                          pvalue >= input$enr_pvalue_min &
                                          pvalue <= input$enr_pvalue_max &
                                          padj >= input$enr_padj_min &
                                          padj <= input$enr_padj_max &
                                          qvalue >= input$enr_qvalue_min &
                                          qvalue <= input$enr_qvalue_max &
                                          Count >= input$enr_count_min &
                                          Count <= input$enr_count_max
        )

      } else {
        #no filtering enabled, show all values and update filtering options
         filtered_enrichment_results <<- enrichment_results

        #update numericInputs for data filtering
        updateNumericInput(session, "enr_pvalue_min", value = min(enrichment_results$pvalue, na.rm=T))
        updateNumericInput(session, "enr_pvalue_max", value = max(enrichment_results$pvalue, na.rm=T))
        updateNumericInput(session, "enr_padj_min", value = min(enrichment_results$padj, na.rm=T))
        updateNumericInput(session, "enr_padj_max", value = max(enrichment_results$padj, na.rm=T))
        updateNumericInput(session, "enr_qvalue_min", value = min(enrichment_results$qvalue, na.rm=T))
        updateNumericInput(session, "enr_qvalue_max", value = max(enrichment_results$qvalue, na.rm=T))
        updateNumericInput(session, "enr_count_min", value = min(enrichment_results$Count, na.rm=T))
        updateNumericInput(session, "enr_count_max", value = max(enrichment_results$Count, na.rm=T))

      }
      
      #show table
      filtered_enrichment_results
    })
    
  })
  
  #download Pathway enrichment results table
  output$downloadenrichmentTable <- downloadHandler(
    filename = function() {
      paste("Pathway enrichment results.csv")
    },
    content = function(file) {
      #write to file
      write.csv(filtered_enrichment_results, file, row.names = FALSE)
    }
  )
  
  #download GSEA results table
  output$downloadgseaTable <- downloadHandler(
    filename = function() {
      paste("GSEA results.csv")
    },
    content = function(file) {
      #write to file
      write.csv(filtered_gsea_results, file, row.names = FALSE)
    }
  )
  
  #function to load html file for start info tab
  output$startinfo <- renderUI ({
    return(includeHTML("startinfo.html"))
  })
  
  #function to load html file for help info tab
  output$helpinfo <- renderUI ({
    return(includeHTML("helpinfo.html"))
  })
  
  ### HELPER FUNCTIONS ###
  #function to get unique condition/treatment names
  get_unique_conds <- function(target_col)({
    
    #use a random plot table sample
    temp <- plotCounts(dds, gene=listofgenes[1,2], intgroup=c("condition", "replicate"), returnData=TRUE)

    #empty vector to hold unique conditions
    newlist <- NULL
    
    newtemp <- as.matrix(temp)
    
    newlist[1] <- newtemp[1,target_col]
    
    counter = 2
    for (i in 2:nrow(newtemp)){
      
      prev = newtemp[i-1,target_col]
      curr = newtemp[i,target_col]
      
      if (curr != prev){
        newlist[counter] <- curr
        counter = counter+1
      }
    }
    
    return(newlist)
    
  })  
  
  #function to dynamically draw colorboxes for n number of replicates
  dynamic_colorgen <- function(widget_name, selector, sourcelist)({
    
    #get list of unique condition names
    unique_IDs <- get_unique_conds(target_col = sourcelist)
    
    #initialize vector to hold dynamically generated names for color widgets
    colorwidget_names <- NULL
    
    if (sourcelist == "condition")
      num = num_conditions
    else if (sourcelist == "replicate")
      num = num_replicates
    
    for (count in 1:num){
      
      #prepare unique identifier name for current widget in the loop
      colorwidget_names[count] <- paste0(widget_name, count)
      
      insertUI(
        selector = selector,
        ui = colourInput(colorwidget_names[count], unique_IDs[count], "#EB4141", allowTransparent = FALSE)
      )
    }
    
  })
  
  #function to transform R exponent notation to human readable scientific notation
  to_scientific <- function(num)({
    #convert to scientific notation string
    num <- format(num, digits=3, scientific = TRUE)
    
    #make sure a zero value prints at "0"
    num <- gsub("0e\\+00","0",num)
    
    #quote the part before the exponent to keep all the digits
    num <- gsub("^(.*)e", "'\\1'e", num)
    
    #replace exponent "e" with proper "10^" notation
    num <- gsub("e", "%*%10^", num)
    
    # return string
    parse(text=num)
  })

  ### START PCA PLOT SAVE AS FUNCTIONS ### 
  # png
  output$savePCApng <- downloadHandler(
    filename = function() { paste0("PCA_plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_PCA_plot_width"]], 
          height = session$clientData[["output_PCA_plot_height"]])
      print(do_PCA_plot())
      dev.off()
    })
  
  # jpg
  output$savePCAjpg <- downloadHandler(
    filename = function() { paste0("PCA_plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_PCA_plot_width"]], 
           height = session$clientData[["output_PCA_plot_height"]])
      print(do_PCA_plot())
      dev.off()
    })
  
  # svg
  output$savePCAsvg <- downloadHandler(
    filename = function() { paste0("PCA_plot", '.svg') },
    content = function(file) {
      svg(file, 
          width = session$clientData[["output_PCA_plot_width"]]/as.integer(input$PCAplot_dpi), 
          height = session$clientData[["output_PCA_plot_height"]]/as.integer(input$PCAplot_dpi))
      print(do_PCA_plot())
      dev.off()
    })
  
  # pdf
  output$savePCApdf <- downloadHandler(
    filename = function() { paste0("PCA_plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_PCA_plot_width"]]/as.integer(input$PCAplot_dpi), 
          height = session$clientData[["output_PCA_plot_height"]]/as.integer(input$PCAplot_dpi))
      print(do_PCA_plot())
      dev.off()
    })
  
  # tiff
  output$savePCAtiff <- downloadHandler(

    filename = function() { paste0("PCA_plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_PCA_plot_width"]], 
           height = session$clientData[["output_PCA_plot_height"]])
      print(do_PCA_plot())
      dev.off()
    })
  ### END PCA PLOT SAVE AS FUNCTIONS ###
  
  ### START SAMPLE CLUSTERING PLOT SAVE AS FUNCTIONS ### 
  # png
  output$saveClusteringpng <- downloadHandler(
    filename = function() { paste0("SampleClustering_plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_sampleClustering_plot_width"]], 
          height = session$clientData[["output_sampleClustering_plot_height"]])
      print(do_sampleClustering_plot2())
      dev.off()
    })
  
  # jpg
  output$saveClusteringjpg <- downloadHandler(
    filename = function() { paste0("SampleClustering_plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_sampleClustering_plot_width"]], 
           height = session$clientData[["output_sampleClustering_plot_height"]])
      print(do_sampleClustering_plot2())
      dev.off()
    })
  
  # svg
  output$saveClusteringsvg <- downloadHandler(
    filename = function() { paste0("SampleClustering_plot", '.svg') },
    content = function(file) {
      svg(file, 
          width = session$clientData[["output_sampleClustering_plot_width"]]/as.integer(input$sampleClustering_dpi), 
          height = session$clientData[["output_sampleClustering_plot_height"]]/as.integer(input$sampleClustering_dpi))
      print(do_sampleClustering_plot2())
      dev.off()
    })
  
  # pdf
  output$saveClusteringpdf <- downloadHandler(
    filename = function() { paste0("SampleClustering_plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_sampleClustering_plot_width"]]/as.integer(input$sampleClustering_dpi), 
          height = session$clientData[["output_sampleClustering_plot_height"]]/as.integer(input$sampleClustering_dpi))
      print(do_sampleClustering_plot2())
      dev.off()
    })
  
  # tiff
  output$saveClusteringtiff <- downloadHandler(
    
    filename = function() { paste0("SampleClustering_plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_sampleClustering_plot_width"]], 
           height = session$clientData[["output_sampleClustering_plot_height"]])
      print(do_sampleClustering_plot2())
      dev.off()
    })
  ### END SAMPLE CLUSTERING PLOT SAVE AS FUNCTIONS ###
  
  ### START READ COUNT PLOTS SAVE AS FUNCTIONS ### 
  # png
  output$saveReadCountpng <- downloadHandler(
    filename = function() { paste0("ReadCounts_plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_multi_genecount_plot1_width"]], 
          height = session$clientData[["output_multi_genecount_plot1_height"]])
      print(do_multi_genecount_plot())
      dev.off()
    })
  
  # jpg
  output$saveReadCountjpg <- downloadHandler(
    filename = function() { paste0("ReadCounts_plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_multi_genecount_plot1_width"]], 
           height = session$clientData[["output_multi_genecount_plot1_height"]])
      print(do_multi_genecount_plot())
      dev.off()
    })
  
  # svg
  output$saveReadCountsvg <- downloadHandler(
    filename = function() { paste0("ReadCounts_plot", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_multi_genecount_plot1_width"]]/as.integer(input$ReadCount_dpi), 
          height = session$clientData[["output_multi_genecount_plot1_height"]]/as.integer(input$ReadCount_dpi))
      print(do_multi_genecount_plot())
      dev.off()
    })
  
  # pdf
  output$saveReadCountpdf <- downloadHandler(
    filename = function() { paste0("ReadCounts_plot", '.pdf') },
    content = function(file) {
      pdf(file,
          width = session$clientData[["output_multi_genecount_plot1_width"]]/input$ReadCount_dpi, 
          height = session$clientData[["output_multi_genecount_plot1_height"]]/input$ReadCount_dpi)
      print(do_multi_genecount_plot())
      dev.off()
    })
  
  # tiff
  output$saveReadCounttiff <- downloadHandler(
    
    filename = function() { paste0("ReadCounts_plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_multi_genecount_plot1_width"]], 
           height = session$clientData[["output_multi_genecount_plot1_height"]])
      print(do_multi_genecount_plot())
      dev.off()
    })
  ### END READ COUNT PLOTS SAVE AS FUNCTIONS ###
  
  ### START HEATMAP SAVE AS FUNCTIONS ### 
  # png
  output$saveHeatmappng <- downloadHandler(
    filename = function() { paste0("Heatmap_plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_countMatrix_heatmap_width"]], 
          height = session$clientData[["output_countMatrix_heatmap_height"]])
      print(do_countMatrix_heatmap2())
      dev.off()
    })
  
  # jpg
  output$saveHeatmapjpg <- downloadHandler(
    filename = function() { paste0("Heatmap_plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_countMatrix_heatmap_width"]], 
           height = session$clientData[["output_countMatrix_heatmap_height"]])
      print(do_countMatrix_heatmap2())
      dev.off()
    })
  
  # svg
  output$saveHeatmapsvg <- downloadHandler(
    filename = function() { paste0("Heatmap_plot", '.svg') },
    content = function(file) {
      svg(file, 
          width = session$clientData[["output_countMatrix_heatmap_width"]]/as.integer(input$Heatmap_dpi), 
          height = session$clientData[["output_countMatrix_heatmap_height"]]/as.integer(input$Heatmap_dpi))
      print(do_countMatrix_heatmap2())
      dev.off()
    })
  
  # pdf
  output$saveHeatmappdf <- downloadHandler(
    filename = function() { paste0("Heatmap_plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_countMatrix_heatmap_width"]]/as.integer(input$Heatmap_dpi), 
          height = session$clientData[["output_countMatrix_heatmap_height"]]/as.integer(input$Heatmap_dpi))
      print(do_countMatrix_heatmap2())
      dev.off()
    })
  
  # tiff
  output$saveHeatmaptiff <- downloadHandler(
    
    filename = function() { paste0("Heatmap_plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_countMatrix_heatmap_width"]], 
           height = session$clientData[["output_countMatrix_heatmap_height"]])
      print(do_countMatrix_heatmap2())
      dev.off()
    })
  ### END HEATMAP SAVE AS FUNCTIONS ###
  
  ### START VOLCANO PLOT SAVE AS FUNCTIONS ### 
  # png
  output$saveVolcanopng <- downloadHandler(
    filename = function() { paste0("Volcano_plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_volcanoPlot_width"]], 
          height = session$clientData[["output_volcanoPlot_height"]])
      print(do_volcano_plot())
      dev.off()
    })
  
  # jpg
  output$saveVolcanojpg <- downloadHandler(
    filename = function() { paste0("Volcano_plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_volcanoPlot_width"]], 
           height = session$clientData[["output_volcanoPlot_height"]])
      print(do_volcano_plot())
      dev.off()
    })
  
  # svg
  output$saveVolcanosvg <- downloadHandler(
    filename = function() { paste0("Volcano_plot", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_volcanoPlot_width"]]/as.integer(input$Volcanoplot_dpi),
          height = session$clientData[["output_volcanoPlot_height"]]/as.integer(input$Volcanoplot_dpi))
      print(do_volcano_plot())
      dev.off()
    })
  
  # pdf
  output$saveVolcanopdf <- downloadHandler(
    filename = function() { paste0("Volcano_plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_volcanoPlot_width"]]/as.integer(input$Volcanoplot_dpi), 
          height = session$clientData[["output_volcanoPlot_height"]]/as.integer(input$Volcanoplot_dpi))
      print(do_volcano_plot())
      dev.off()
    })
  
  # tiff
  output$saveVolcanotiff <- downloadHandler(
    
    filename = function() { paste0("Volcano_plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_volcanoPlot_width"]], 
           height = session$clientData[["output_volcanoPlot_height"]])
      print(do_volcano_plot())
      dev.off()
    })
  ### END VOLCANO PLOT SAVE AS FUNCTIONS ###
  
  ### START ENRICHMENT PLOT SAVE AS FUNCTIONS ### 
  # png
  output$saveEnrichmentplotpng <- downloadHandler(
    filename = function() { paste0("Enrichment plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_enrichmentPlot_width"]], 
          height = session$clientData[["output_enrichmentPlot_height"]])
      print(do_enrichment_plot())
      dev.off()
    })
  
  # jpg
  output$saveEnrichmentplotjpg <- downloadHandler(
    filename = function() { paste0("Enrichment plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_enrichmentPlot_width"]], 
           height = session$clientData[["output_enrichmentPlot_height"]])
      print(do_enrichment_plot())
      dev.off()
    })
  
  # svg
  output$saveEnrichmentplotsvg <- downloadHandler(
    filename = function() { paste0("Enrichment plot", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_enrichmentPlot_width"]]/as.integer(input$enrichmentplot_dpi),
          height = session$clientData[["output_enrichmentPlot_height"]]/as.integer(input$enrichmentplot_dpi))
      print(do_enrichment_plot())
      dev.off()
    })
  
  # pdf
  output$saveEnrichmentplotpdf <- downloadHandler(
    filename = function() { paste0("Enrichment plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_enrichmentPlot_width"]]/as.integer(input$enrichmentplot_dpi), 
          height = session$clientData[["output_enrichmentPlot_height"]]/as.integer(input$enrichmentplot_dpi))
      print(do_enrichment_plot())
      dev.off()
    })
  
  # tiff
  output$saveEnrichmentplottiff <- downloadHandler(
    
    filename = function() { paste0("Enrichment plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_enrichmentPlot_width"]], 
           height = session$clientData[["output_enrichmentPlot_height"]])
      print(do_enrichment_plot())
      dev.off()
    })
  ### END ENRICHMENT PLOT SAVE AS FUNCTIONS ###
  
  ### START ENRICHMENT MAP SAVE AS FUNCTIONS ### 
  # png
  output$saveEnrichmentmappng <- downloadHandler(
    filename = function() { paste0("Enrichment map", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_enrichmentMap_width"]], 
          height = session$clientData[["output_enrichmentMap_height"]])
      print(do_enrichment_map())
      dev.off()
    })
  
  # jpg
  output$saveEnrichmentmapjpg <- downloadHandler(
    filename = function() { paste0("Enrichment map", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_enrichmentMap_width"]], 
           height = session$clientData[["output_enrichmentMap_height"]])
      print(do_enrichment_map())
      dev.off()
    })
  
  # svg
  output$saveEnrichmentmapsvg <- downloadHandler(
    filename = function() { paste0("Enrichment map", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_enrichmentMap_width"]]/as.integer(input$enrichmentmap_dpi),
          height = session$clientData[["output_enrichmentMap_height"]]/as.integer(input$enrichmentmap_dpi))
      print(do_enrichment_map())
      dev.off()
    })
  
  # pdf
  output$saveEnrichmentmappdf <- downloadHandler(
    filename = function() { paste0("Enrichment map", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_enrichmentMap_width"]]/as.integer(input$enrichmentmap_dpi), 
          height = session$clientData[["output_enrichmentMap_height"]]/as.integer(input$enrichmentmap_dpi))
      print(do_enrichment_map())
      dev.off()
    })
  
  # tiff
  output$saveEnrichmentmaptiff <- downloadHandler(
    
    filename = function() { paste0("Enrichment map", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_enrichmentMap_width"]], 
           height = session$clientData[["output_enrichmentMap_height"]])
      print(do_enrichment_map())
      dev.off()
    })
  ### END ENRICHMENT MAP SAVE AS FUNCTIONS ###
  
  ### START GSEA MAP SAVE AS FUNCTIONS ### 
  # png
  output$saveGSEAmappng <- downloadHandler(
    filename = function() { paste0("GSEA map", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_gseaMap_width"]], 
          height = session$clientData[["output_gseaMap_height"]])
      print(do_gsea_map())
      dev.off()
    })
  
  # jpg
  output$saveGSEAmapjpg <- downloadHandler(
    filename = function() { paste0("GSEA map", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_gseaMap_width"]], 
           height = session$clientData[["output_gseaMap_height"]])
      print(do_gsea_map())
      dev.off()
    })
  
  # svg
  output$saveGSEAmapsvg <- downloadHandler(
    filename = function() { paste0("GSEA map", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_gseaMap_width"]]/as.integer(input$gseamap_dpi),
          height = session$clientData[["output_gseaMap_height"]]/as.integer(input$gseamap_dpi))
      print(do_gsea_map())
      dev.off()
    })
  
  # pdf
  output$saveGSEAmappdf <- downloadHandler(
    filename = function() { paste0("GSEA map", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_gseaMap_width"]]/as.integer(input$gseamap_dpi), 
          height = session$clientData[["output_gseaMap_height"]]/as.integer(input$gseamap_dpi))
      print(do_gsea_map())
      dev.off()
    })
  
  # tiff
  output$saveGSEAmaptiff <- downloadHandler(
    
    filename = function() { paste0("GSEA map", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_gseaMap_width"]], 
           height = session$clientData[["output_gseaMap_height"]])
      print(do_gsea_map())
      dev.off()
    })
  ### END GSEA MAP SAVE AS FUNCTIONS ###
  
  ### START GSEA PLOT SAVE AS FUNCTIONS ### 
  # png
  output$saveGSEAplotpng <- downloadHandler(
    filename = function() { paste0("GSEA plot", '.png') },
    content = function(file) {
      png(file, 
          width = session$clientData[["output_gseaPlot_width"]], 
          height = session$clientData[["output_gseaPlot_height"]])
      print(do_gsea_plot())
      dev.off()
    })
  
  # jpg
  output$saveGSEAplotjpg <- downloadHandler(
    filename = function() { paste0("GSEA plot", '.jpg') },
    content = function(file) {
      jpeg(file, 
           width = session$clientData[["output_gseaPlot_width"]], 
           height = session$clientData[["output_gseaPlot_height"]])
      print(do_gsea_plot())
      dev.off()
    })
  
  # svg
  output$saveGSEAplotsvg <- downloadHandler(
    filename = function() { paste0("GSEA plot", '.svg') },
    content = function(file) {
      svg(file,
          width = session$clientData[["output_gseaPlot_width"]]/as.integer(input$gseaplot_dpi),
          height = session$clientData[["output_gseaPlot_height"]]/as.integer(input$gseaplot_dpi))
      print(do_gsea_plot())
      dev.off()
    })
  
  # pdf
  output$saveGSEAplotpdf <- downloadHandler(
    filename = function() { paste0("GSEA plot", '.pdf') },
    content = function(file) {
      pdf(file, 
          width = session$clientData[["output_gseaPlot_width"]]/as.integer(input$gseaplot_dpi), 
          height = session$clientData[["output_gseaPlot_height"]]/as.integer(input$gseaplot_dpi))
      print(do_gsea_plot())
      dev.off()
    })
  
  # tiff
  output$saveGSEAplottiff <- downloadHandler(
    
    filename = function() { paste0("GSEA plot", '.tiff') },
    content = function(file) {
      tiff(file, 
           width = session$clientData[["output_gseaPlot_width"]], 
           height = session$clientData[["output_gseaPlot_height"]])
      print(do_gsea_plot())
      dev.off()
    })
  ### END GSEA PLOT SAVE AS FUNCTIONS ###

})


