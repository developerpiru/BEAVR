#install all required R packages for BEAVR

installCran <- function(package_name){
  install.packages(package_name)
}

installBioc <- function(package_name){
  BiocManager::install(package_name)
}

#may need to run this
options(repos = "https://cloud.r-project.org")
options(install.packages.compile.from.source="interactive")

## load required libraries
installCran("shiny")
installCran("shinydashboard")
installCran("shinyWidgets")

#cran packages
installCran("BiocManager")
installCran("colourpicker")
installCran("data.table")
installCran("devtools") # to install from github
installCran("DT")
installCran("ggplot2")
installCran("ggpubr")
installCran("ggrepel")
installCran("ggraph")
installCran("gridExtra")
installCran("pheatmap")
installCran("RColorBrewer")
installCran("scales")
installCran("shiny")
installCran("shinydashboard")
installCran("shinyjqui")
installCran("shinyWidgets")
installCran("shinycssloaders")
installCran("circlize")

# #Bioconductor packages
installBioc("DESeq2")
installBioc("vsn")
installBioc("apeglm")
installBioc("org.Hs.eg.db")
installBioc("org.Mm.eg.db")
installBioc("ReactomePA")
installBioc("enrichplot")

# #GitHub packages
devtools::install_github("kevinblighe/EnhancedVolcano")
devtools::install_github("jokergoo/ComplexHeatmap")