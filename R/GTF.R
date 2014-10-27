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
#' gene = gencode$toBed(category = "gene")
#' exon = gencode$toBed(category = "exon")
#' gencode$toBed(file = "gene.bed", category = "gene")
#' gencode$toBed(file = "exon.bed", category = "exon")
#' gi = gencode$getGeneID(1:10)
#' chr = gencode$getValueByGeneID(gi, type = "chr")
#' tr = gencode$getTranscriptsByGeneID(gi)
#' gencode_subset = gencode$subset(1:10)
#' gencode_subset = gencode$subset(gi[1:5])
#' gencode$subset(1:10, copy = FALSE)
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
		qqcat("It contains nothing.\n", cat_prefix = "")
	} else {
		qqcat("It contains @{length(.self$gtf)} genes\n", cat_prefix = "")
		if(.self$sorted) {
			qqcat("Has been sorted.\n", cat_prefix = "")
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


GTF$methods(toBed = function(file = NULL, category = c("gene", "exon", "transcript", "tss", "upstream", "downstream"), 
	extend = 2000, type = NULL) {
	"convert and write to BED file. If gtf is merged, or if gtf is not merged and category == 'gene', possible values for 'type' should be in 'pseudogene', 'lincRNA', 'protein_coding', 'antisense', 'processed_transcript', 'snRNA', 'sense_intronic', 'miRNA', 'misc_RNA', 'snoRNA', 'rRNA', 'polymorphic_pseudogene', 'sense_overlapping', '3prime_overlapping_ncrna', 'IG_V_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_pseudogene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_V_pseudogene', 'IG_C_pseudogene', 'TR_D_gene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'IG_D_gene', 'Mt_tRNA', 'Mt_rRNA', if 'gtf' is not merged, and 'category' is not equal to `gene`, possible values for type should be in 'processed_transcript', 'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene', 'lincRNA', 'protein_coding', 'processed_pseudogene', 'antisense', 'snRNA', 'pseudogene', 'retained_intron', 'nonsense_mediated_decay', 'sense_intronic', 'miRNA', 'misc_RNA', 'transcribed_processed_pseudogene', 'snoRNA', 'non_stop_decay', 'rRNA', 'unitary_pseudogene', 'polymorphic_pseudogene', 'sense_overlapping', 'TEC', '3prime_overlapping_ncrna', 'IG_V_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_pseudogene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_V_pseudogene', 'IG_C_pseudogene', 'TR_D_gene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'IG_D_gene', 'Mt_tRNA', 'Mt_rRNA'"
	
	category = match.arg(category)

	if(category == "gene") {
		chr = sapply(.self$gtf, function(x) x$chr)
		start = sapply(.self$gtf, function(x) x$start)
		end = sapply(.self$gtf, function(x) x$end)
		id = names(.self$gtf)
		value = rep(0, length(id))
		strand = sapply(.self$gtf, function(x) x$strand)
		gt = sapply(.self$gtf, function(x) x$type)

	} else if(category == "exon") {

		n_exon = sum(sapply(.self$gtf, function(x) {
			sum(sapply(x$transcript, function(tr) length(tr$exon)))
		}))
		qqcat("totally @{n_exon} exons\n")
		chr = character(n_exon)
		start = integer(n_exon)
		end = integer(n_exon)
		id = character(n_exon)
		value = integer(n_exon)
		strand = character(n_exon)
		gt = character(n_exon)
		gi = names(.self$gtf)
		n = 0
		for(i in seq_along(gi)) {
			for(k in seq_along(.self$gtf[[i]]$transcript)) {
				exon = .self$gtf[[i]]$transcript[[k]]$exon
				l_exon = length(exon)
				chr[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$chr, l_exon)
				start[n + seq_len(l_exon)] = sapply(exon, function(x) x$start)
				end[n + seq_len(l_exon)] = sapply(exon, function(x) x$end)
				id[n + seq_len(l_exon)] = paste(gi[i], k, seq_along(exon), sep = "_")
				strand[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$strand, l_exon)
				gt[n+seq_len(l_exon)] = rep(.self$gtf[[i]]$type, l_exon)
				n = n + l_exon
			}
			
			if(i %% 500 == 0) {
				qqcat("@{i}/@{length(gi)} genes finished\n")
			}
		}

		value = rep(0, length(id))
		
	} else if(category %in% c("transcript", "tss", "upstream", "downstream")) {

		n_tss = sum(sapply(.self$gtf, function(x) {
			length(x$transcript)
		}))

		qqcat("totally @{n_tss} @{category}s\n")
		chr = character(n_tss)
		start = integer(n_tss)
		end = integer(n_tss)
		id = character(n_tss)
		value = integer(n_tss)
		strand = character(n_tss)
		gt = character(n_tss)
		gi = names(.self$gtf)
		n = 0
		for(i in seq_along(gi)) {
			for(k in seq_along(.self$gtf[[i]]$transcript)) {
				transcript_start = .self$gtf[[i]]$transcript[[k]]$start
				transcript_end = .self$gtf[[i]]$transcript[[k]]$end
				
				if(.self$gtf[[i]]$strand == "+") {
					if(category == "tss") {
						pos_start = transcript_start
						pos_end = transcript_start
					} else if(category == "upstream") {
						pos_start = transcript_start - 1 - extend + 1
						pos_end = transcript_start - 1
					} else if(category == "downstream") {
						pos_start = transcript_end + 1
						pos_end = transcript_end + 1 + extend - 1
					} else {
						pos_start = transcript_start
						pos_end = transcript_end
					}
				} else {
					if(category == "tss") {
						pos_start = transcript_end
						pos_end = transcript_end
					} else if(category == "upstream") {
						pos_start = transcript_end + 1
						pos_end = transcript_end + 1 + extend - 1
					} else if(category == "downstream") {
						pos_start = transcript_start - 1 - extend + 1
						pos_end = transcript_start - 1
					} else {
						pos_start = transcript_start
						pos_end = transcript_end
					}
				}
				chr[n + 1] = .self$gtf[[i]]$chr
				start[n + 1] = pos_start
				end[n + 1] = pos_end
				id[n + 1] = paste(gi[i], k, sep = "_")
				strand[n + 1] = .self$gtf[[i]]$strand
				gt[n + 1] = .self$gtf[[i]]$type
				n = n + 1
			}
			
			if(i %% 500 == 0) {
				qqcat("\r")
				qqcat("@{i}/@{length(gi)} genes finished           \n")
				flush.console()
			}
		}

		value = rep(0, length(id))
		
	} 

	start = as.integer(start)
	end = as.integer(end)
	df = data.frame(chr = chr, start = start, end = end, name = id, value = value, strand = strand, type = gt, stringsAsFactors = FALSE)
	df = subset(df, start <= end)
	if(!is.null(type)) {
		df = subset(df, type %in% type)
	}
	df = df[order.bed(df), ]

	rownames(df) = NULL

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

GTF$methods(subset = function(index, copy = TRUE) {
	"subset"
	if(copy) {
		new_gtf = .self$copy()
		new_gtf$gtf = new_gtf$gtf[index]
		new_gtf$sorted = FALSE
		return(new_gtf)
	} else {
		.self$gtf = .self$gtf[index]
		.self$sorted = FALSE
		return(.self)
	}
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
