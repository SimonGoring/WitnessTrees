all: markdown word pdf

pdf: Neotoma_paper.Rmd
	Rscript -e "library('rmarkdown'); render('GoringetalMarkdown.Rmd', 'pdf_document')"

word: Neotoma_paper.Rmd
	Rscript -e "library('rmarkdown'); render('GoringetalMarkdown.Rmd', 'word_document')"

markdown: Neotoma_paper.Rmd
	Rscript -e "library('rmarkdown'); render('GoringetalMarkdown.Rmd', 'md_document')"
