#!/usr/local/packages/r-3.4.0/bin Rscript
library(rmarkdown)
library(knitr)
args = commandArgs(trailingOnly=TRUE)
arg0=args[1]
arg1=args[2]
arg2=args[3]
arg3=args[4]
#stopifnot(length(args) ==4)
arg4=args[5]
stopifnot(length(args) ==5)
#print(arg2)
rmarkdown::render(arg0, params=list( projectname=arg1, pathd=arg2, infof=arg3), output_dir = arg4, output_file = paste(arg1,"_", Sys.Date(), '_Alignment_Report.pdf', sep=''))
