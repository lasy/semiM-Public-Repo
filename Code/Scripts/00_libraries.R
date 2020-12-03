
package_list = c('data.table',
                 'cowplot',
                 'ggpubr',
                 'tictoc',
                 'igraph',
                 'shiny',
                 'knitr','bookdown','styler', 
                 'kableExtra',
                 'feather',
                 'tidyverse',
                 'ggthemes',
                 'magrittr',
                 'sn'
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

library(HiddenSemiMarkov)


