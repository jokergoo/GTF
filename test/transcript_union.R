###############################################################
##
## transfrom gtf to json
## read raw gencode
## merge multiple transcripts in single gene
## save origin and merged gencode as RData
## write merged gencode as gtf (for alignment and counting)
## generate gene length
## generate gene.bed and exon.bed
##
###############################################################

library(GetoptLong)
setwd(dirname(get_scriptname()))
setwd("../")
source("ROO/Class/gencode.R")

gencode = Gencode$new()
gencode$readGTF("gencode/gencode.v17_broad.gtf")
save(gencode, file = "RData/gencode_broad.RData")

qqcat("merge multiple transcripts for single gene\n")
gencode$mergeTranscripts()
save(gencode, file = "RData/gencode_broad_merged.RData")


gencode$writeGTF("gencode/gencode_broad_merged.gtf")


qqcat("calculate sum of length of merged exons in genes\n")
gene_length = gencode$geneLength("exon")
df = data.frame(id = names(gene_length), gene_length = gene_length)
write.table(df, file = "gencode/gencode_gene_length.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# gene.bed and exon.bed are already sorted by chr and start positions
gencode$toBed(file = "bed/gene.bed", type = "gene")
gencode$toBed(file = "bed/exon.bed", type = "exon")
gencode$toBed(type = "5UTR")
gencode$toBed(type = "3UTR")
