installPkgs <- function(package_name, bioc){
  if (bioc == FALSE)
    install.packages(package_name)
  else if (bioc == TRUE) #install using Bioconductor package manager
    BiocManager::install(package_name)
}

## load required libraries
installPkgs("shiny")
installPkgs("shinydashboard")
installPkgs("shinyWidgets")

#cran packages
installPkgs("BiocManager")
installPkgs("colourpicker")
installPkgs("data.table")
installPkgs("devtools") # to install from github
installPkgs("DT")
installPkgs("ggplot2")
installPkgs("ggpubr")
installPkgs("ggrepel")
installPkgs("ggraph")
installPkgs("gridExtra")
installPkgs("pheatmap")
installPkgs("RColorBrewer")
installPkgs("scales")
installPkgs("shiny")
installPkgs("shinydashboard")
installPkgs("shinyjqui")
installPkgs("shinyWidgets")
installPkgs("shinycssloaders")
installPkgs("circlize")

# #Bioconductor packages
installPkgs("DESeq2", bioc == TRUE)
installPkgs("vsn", bioc == TRUE)
installPkgs("apeglm", bioc == TRUE)
installPkgs("org.Hs.eg.db", bioc == TRUE)
installPkgs("org.Mm.eg.db", bioc == TRUE)
installPkgs("ReactomePA", bioc == TRUE)
installPkgs("enrichplot", bioc == TRUE)

# #GitHub packages
devtools::install_github("kevinblighe/EnhancedVolcano")
devtools::install_github("jokergoo/ComplexHeatmap")