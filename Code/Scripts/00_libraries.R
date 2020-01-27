
package_list = c('knitr','bookdown','styler',
                 'tidyverse','magrittr','stringr',
                 'readr','feather',
                 'plyr','dplyr',
                 'data.table',
                 'SDMTools',
                 'ggplot2','ggthemes','hexbin','mapdata','gridExtra','quantreg','grid',
                 'plotly','cowplot',
                 'scales','reshape',
                 'HMM',
                 'chron','lubridate',
                 'foreach',
                 'parallel','doParallel',
                 'tictoc',
                 'mhsmm',
                 'igraph',
                 'shiny'
)

pckgs = installed.packages()

for(package in package_list){
  cat(package,"\n")
  pckgs = installed.packages()
  need.to.install = (!(package %in% pckgs[,1]))
  if(need.to.install){cat('installing package ',package,'\n');install.packages(package,repos = "https://cloud.r-project.org",dependencies = TRUE)}
  library(package = package, character.only = TRUE)
}

rm(package_list, pckgs,need.to.install,package)


