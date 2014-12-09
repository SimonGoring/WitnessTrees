
library(devtools)
install_github('rmarkdown', 'rstudio')

library(rmarkdown)
render('GoringetalMarkdown.Rmd', 'word_document')
render('GoringetalMarkdown.Rmd', 'pdf_document')
