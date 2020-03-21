#!/usr/bin/Rscript

#rmarkdown::render ("gwas-markdown.Rmd", output_file="out.pdf", output_format="pdf_document", 
#				   params=list (workingDir=getwd(), reportTitle="Naive report"))

rmarkdown::render ("gwas-markdown.Rmd", output_file="out.html", output_format="html_document", 
				   params=list (workingDir=getwd(), reportTitle="Naive report", nBest=5))

