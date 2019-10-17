#!/usr/local/packages/r-3.4.0/bin Rscript
library(rmarkdown)
library(knitr)
args = commandArgs(trailingOnly=TRUE)
arg1=args[1]
arg2=args[2]
arg3=args[3]
arg4=args[4]
stopifnot(length(args) ==4)
rmarkdown::render("/local/projects/RNASEQ/Report_Generation/Scripts/v1.0/FastQC_Report.Rmd", params=list( projectname=arg1, pathd=arg2, infof=arg3), output_dir = arg4, output_file = paste(arg1,"_", Sys.Date(), '_FastQC_Report.pdf', sep=''))
#rmarkdown::render("/local/projects/RNASEQ/Report_Generation/Scripts/v1.0/FastQC_Report.Rmd", params=list( projectname=arg1, pathd=arg2, infof=arg3), output_dir = arg4)
