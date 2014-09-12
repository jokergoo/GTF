######################################
# 
# Gencode class
#
######################################


#' Reference class to manipulate GTF data
#'
#' @name GTF-class
#' @rdname GTF-class
#' @import methods
#' @import rjson
#' @import IRanges
#' @import GetoptLong
#' @import parallel
#' @exportClass GTF
#' @field GTF complicated structure of list
#' @field sorted data is sorted or not
#' @examples
#' \dontrun{
#' library(GTF)
#' gencode = new("GTF")
#' gencode$read(paste(system.file(package = "GTF"), "/extdata/test.gtf", sep=""))
#' gencode$mergeTranscripts()
#' gencode$mergeTranscripts(mc.cores = 2)
#' gencode$write("output.gtf")
#' gene_length = gencode$geneLength("exon")
#' gene = gencode$toBed(type = "gene")
#' exon = gencode$toBed(type = "exon")
#' gencode$toBed(file = "gene.bed", type = "gene")
#' gencode$toBed(file = "exon.bed", type = "exon")
#' gi = gencode$getGeneID(1:10)
#' chr = gencode$getValueByGeneID(gi, type = "chr")
#' tr = gencode$getTranscriptsByGeneID(gi)
#' }
GTF = setRefClass("GTF",
	fields = list(
		gtf = "list",
		sorted  = "logical")
)

GTF$methods(initialize = function() {
	gtf <<- vector("list", length = 0)
	sorted <<- FALSE
})

# common methods
GTF$methods(show = function() {
	"show method"
	if(length(.self$gtf) == 0) {
		qqcat("It contains nothing.\n")
	} else {
		qqcat("It contains @{length(.self$gtf)} genes\n")
		if(.self$sorted) {
			qqcat("Has been sorted.\n")
		}
	}
})

GTF$methods(read = function(gtf) {
	"read from GTF file"
	
	perlScript = qq("@{system.file(package = 'GTF')}/extdata/gtf_to_json.pl")
	cmd = qq("perl \"@{perlScript}\" \"@{normalizePath(gtf)}\"")
	p = pipe(cmd)
	gtf <<- fromJSON(file = p)
	close(p)
	.self$sort()
	
	return(invisible(.self))
})


GTF$methods(sort = function() {
	"sort GTF"
	
	chr_all = sapply(.self$gtf, function(x) x$chr)
	
	g2 = list()
	for(chr in sort_chr(unique(chr_all))) {

		gene_list = .self$gtf[chr_all == chr]
		#qqcat("@{chr} contains @{length(gene_list)} genes ...\n")
		gene_list = sort_by_start(gene_list)

		for(i in seq_along(gene_list)) {
			gene = gene_list[[i]]
			transcript_list = gene$transcript
			transcript_list = sort_by_start(transcript_list)

			for(j in seq_along(transcript_list)) {
				transcript = transcript_list[[j]]
				exon_list = transcript$exon
				exon_list = sort_by_start(exon_list)
				CDS_list = transcript$CDS
				CDS_list = sort_by_start(CDS_list)
				transcript_list[[j]]$exon = exon_list
				transcript_list[[j]]$CDS = CDS_list
			}

			gene_list[[i]]$transcript = transcript_list
		}

		g2 = c(g2, gene_list)
	}
	gtf <<- g2
	sorted <<- TRUE
	
	return(invisible(.self))
})


GTF$methods(toBed = function(file = NULL, type = c("gene", "exon"), chromosome = paste("chr", c(1:22, "X", "Y"), sep = "")) {
	"convert and write to BED file"
	type = match.arg(type)

	if(type == "gene") {
		chr = sapply(.self$gtf, function(x) x$chr)
		start = sapply(.self$gtf, function(x) x$start)
		end = sapply(.self$gtf, function(x) x$end)
		id = names(.self$gtf)
		value = rep(0, length(id))
		strand = sapply(.self$gtf, function(x) x$strand)

	} else if(type == "exon") {

		n_exon = sum(sapply(.self$gtf, function(x) length(x$transcript[[1]]$exon)))
		qqcat("totally @{n_exon} exons\n")
		chr = character(n_exon)
		start = integer(n_exon)
		end = integer(n_exon)
		id = character(n_exon)
		value = integer(n_exon)
		strand = character(n_exon)
		gi = names(.self$gtf)
		n = 0
		for(i in seq_along(gi)) {
			if(length(.self$gtf[[i]]$transcript) > 1) {
				stop(qq("find more than one transcripts in @{gi[i]}, please merge your `gtf` first.\n"))
			}

			exon = .self$gtf[[i]]$transcript[[1]]$exon
			l_exon = length(exon)
			chr[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$chr, l_exon)
			start[n + seq_len(l_exon)] = sapply(exon, function(x) x$start)
			end[n + seq_len(l_exon)] = sapply(exon, function(x) x$end)
			id[n + seq_len(l_exon)] = paste(gi[i], seq_along(exon), sep = "_")
			strand[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$strand, l_exon)
			n = n + l_exon

			if(i %% 500 == 0) {
				qqcat("@{i}/@{length(gi)} genes finished\n")
			}
		}

		value = rep(0, length(id))
		
	} else {
		stop("Currently only support 'gene' and 'exon'\n")
	}

	l = chr %in% chromosome
	chr = chr[l]
	start = start[l]
	end = end[l]
	id = id[l]
	value = value[l]
	strand = strand[l]
	
	start = as.integer(start)
	end = as.integer(end)
	df = data.frame(chr = chr, start = start, end = end, id = id, value = value, strand = strand)
	df = df[order.bed(df), ]
	if(is.null(file)) {
		return(df)
	} else {
		write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		return(invisible(.self))
	}
})


GTF$methods(mergeTranscripts = function(mc.cores = 1) {
	"merge alternative transcripts into one union gene"
	
	gn = names(.self$gtf)
	gcd = mclapply(seq_along(gn), function(i) union_transcript(.self$gtf[[i]], gn[i]), mc.cores = mc.cores)
	names(gcd) = gn
	gtf <<- gcd
	
	.self$sort()
	
	return(invisible(.self))
})


GTF$methods(write = function(file) {
	"write to GTF file"
	out = file(file, "w")
	
	gn = names(.self$gtf)
	
	for(i in seq_along(.self$gtf)) {
		gene = .self$gtf[[i]]
		tn = names(gene$transcript)
		type = sapply(gene$transcript, function(t) t$type)
		last_field = paste("gene_id \"", gn[i], "\"; ", "transcript_id \"", tn[1], "\"; ", "transcript_type \"", type[1], "\";", sep = "")
		writeLines(paste(gene$chr, ".", "gene", gene$start, gene$end, ".", gene$strand, ".", last_field, sep = "\t"), con = out)
		for(j in seq_along(gene$transcript)) {
			last_field = paste("gene_id \"", gn[i], "\"; ", "transcript_id \"", tn[j], "\"; ", "transcript_type \"", type[j], "\";", sep = "")
			writeLines(paste(gene$chr, ".", "transcript", gene$transcript[[j]]$start, gene$transcript[[j]]$end, ".", gene$strand, ".", last_field, sep = "\t"), con = out)
			for(k in seq_along(gene$transcript[[j]]$exon)) {
				exon = gene$transcript[[j]]$exon[[k]]
				last_field = paste("gene_id \"", gn[i], "\"; ", "transcript_id \"", tn[j], "\"; ", "transcript_type \"", type[j], "\"; ", "exon_number ", exon$i, ";", sep = "")
				writeLines(paste(gene$chr, ".", "exon", exon$start, exon$end, ".", gene$strand, ".", last_field, sep = "\t"), con = out)
			}
		}

		if(i %% 1000 == 0) {
			qqcat("@{i} genes finished\n")
		}
	}
	
	close(out)
	
	return(invisible(.self))
})


GTF$methods(geneLength = function(type = c("gene", "exon")) {
	"get gene/exon length"
	type = match.arg(type)
	
	if(type == "gene") {
		return(sapply(.self$gtf, function(x) x$end - x$start + 1))
	} else if(type == "exon") {
		n_tr = sapply(.self$gtf, function(x) length(x$transcript))
		if(any(n_tr > 1)) {
			stop("If you choose type == exon, multiple transcripts for a same gene should be merged first.")
		}
		return(sapply(.self$gtf, function(x) sum(sapply(x$transcript[[1]]$exon, function(y) y$end - y$start + 1))))
	}
})


GTF$methods(getGeneID = function(index = NULL) {
	"get Gene ID (ensembl ID)"
	if(is.null(index)) {
		names(.self$gtf)
	} else {
		names(.self$gtf)[index]
	}
})


GTF$methods(getValueByGeneID = function(geneID, type = c("chr", "name", "strand", "type", "end", "start")) {
	"get values corresponding to gene IDs"
	type = match.arg(type)
	sapply(.self$gtf[geneID], function(x) x[[type]])
})


GTF$methods(getTranscriptsByGeneID = function(geneID) {
	"get transcripts by gene IDs"
	.self$gtf[geneID]
})


# private method
sort_by_start = function(x) {
	if(is.null(x) || length(x) == 0) {
		return(x)
	} else {
		start = sapply(x, function(y) y$start)
		o = order(start)
		return(x[o])
	}
}

sort_chr = function(x) {
	y = gsub("^chr(\\d)$", "chr0\\1", x)
	y = gsub("^chr(\\d)_", "chr0\\1_", y)
	x[order(y)]
}

union_transcript = function(gene, gn = NULL) {
	
	if(length(gene$transcript) == 1) return(gene)
	
	exon = gene$transcript[[1]]$exon
	start = sapply(exon, function(x) x$start)
	end = sapply(exon, function(x) x$end)
	type = gene$type
	
	ir = IRanges(start, end)
	
	for(i in seq_along(gene$transcript)[-1]) {
		exon = gene$transcript[[i]]$exon
		start = sapply(exon, function(x) x$start)
		end = sapply(exon, function(x) x$end)
		
		ir2 = IRanges(start, end)
		ir = union(ir, ir2)
	}
	
	start = start(ir)
	end = end(ir)
	exon = vector("list", length(start))
	names(exon) = paste("exon", seq_along(exon), sep = "_")
	for(i in seq_along(exon)) {
		exon[[i]] = list(start = start[i], end = end[i], i = i)
	}
	
	s = min(start)
	e = max(end)
	
	CDS = gene$transcript[[1]]$CDS
	start = sapply(CDS, function(x) x$start)
	if(length(start) == 0) start = NULL
	end = sapply(CDS, function(x) x$end)
	if(length(end) == 0) end = NULL
	
	ir = IRanges(start, end)
	
	for(i in seq_along(gene$transcript)[-1]) {
		CDS = gene$transcript[[i]]$CDS
		start = sapply(CDS, function(x) x$start)
		end = sapply(CDS, function(x) x$end)
		
		start = sapply(CDS, function(x) x$start)
		if(length(start) == 0) start = NULL
		end = sapply(CDS, function(x) x$end)
		if(length(end) == 0) end = NULL
		
		ir2 = IRanges(start, end)
		ir = union(ir, ir2)
	}
	
	start = start(ir)
	end = end(ir)
	CDS = vector("list", length(start))
	if(length(CDS)) {
		names(CDS) = paste("CDS", seq_along(CDS), sep = "_")
		for(i in seq_along(CDS)) {
			CDS[[i]] = list(start = start[i], end = end[i], i = i)
		}
	}
	
	gene2 = gene
	if(length(CDS)) {
		gene2$transcript = list(transcript_1 = list(exon = exon, CDS = CDS, start = s, end = e, type = type))
	} else {
		gene2$transcript = list(transcript_1 = list(exon = exon, start = s, end = e, type = type))
	}
	
	if(!is.null(gn)) {
		names(gene2$transcript) = gn
	}
	return(gene2)
}

order.bed = function(bed) {
	chr = bed[, 1]
	start = bed[, 2]
	chr = gsub("^chr", "", chr, ignore.case = TRUE)
	unique_chr = unique(chr)
	l = grepl("^\\d+$", unique_chr)
	letter_chr = unique_chr[!l]
	letter_chr_n = as.character(rank(letter_chr) + 1000)
	
	for(i in seq_along(letter_chr)) {
		chr[chr == letter_chr[i]] = letter_chr_n[i]
	}
	chr = as.numeric(chr)
	order(chr, start)
}
