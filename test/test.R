
library(GetoptLong)
library(GenomicRanges)
library(rjson)

source("GTF.R")
gtf = new("GTF")
gtf$read("D:\\personal_project\\GTF\\inst\\extdata\\test.gtf")
gtf$availableTypes(type = "gene")
gtf$availableTypes(type = "transcript")

gtf$genes()
gtf$transcripts()
gtf$exons()
