### Simple class to process GTF data

First initialize a `GTF` class object:

```r
library(GTF)
gencode = new("GTF")
```

Read the Gencode GTF file:

```
gencode$read(paste(system.file(package = "GTF"), "/extdata/test.gtf", sep=""))
```

Normally, a GTF file will provide alternative transcripts for one gene. 
Users can choose to merge all transcripts by taking the union of all exons:

```r
gencode$mergeTranscripts()
gencode$mergeTranscripts(mc.cores = 2)
```

Write back to GTF file:

```r
gencode$write("output.gtf")
```

Get the gene/exon length of each gene:

```r
gene_length = gencode$geneLength("exon")
```

Return or write back to BED file:

```r
gene = gencode$toBed(type = "gene")
exon = gencode$toBed(type = "exon")
gencode$toBed(file = "gene.bed", type = "gene")
gencode$toBed(file = "exon.bed", type = "exon")
```

Query in Gencode data:

```r
gi = gencode$getGeneID(1:10)
chr = gencode$getValueByGeneID(gi, type = "chr")
tr = gencode$getTranscriptsByGeneID(gi)
```

For other gene annotation which is not from Gencode, users can use http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format to generate files in GTF format (e.g. refGene.gtf)
