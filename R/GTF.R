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
		sorted  = "logical",
		merged = "logical")
)

GTF$methods(initialize = function() {
	gtf <<- vector("list", length = 0)
	sorted <<- FALSE
	merged <<- FALSE
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
		if(.self$merged) {
			qqcat("Has been merged.\n", cat_prefix = "")
		}
	}
})

#
GTF$methods(read = function(gtf, gene_fields = c(), transcript_fields = c()) {
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

GTF$methods(availableTypes = function() {
	unique(.self$getValueByGeneID(type = "type"))
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
			ti = names(.self$gtf[[i]]$transcript)
			for(k in seq_along(.self$gtf[[i]]$transcript)) {
				exon = .self$gtf[[i]]$transcript[[k]]$exon
				l_exon = length(exon)
				chr[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$chr, l_exon)
				start[n + seq_len(l_exon)] = sapply(exon, function(x) x$start)
				end[n + seq_len(l_exon)] = sapply(exon, function(x) x$end)
				id[n + seq_len(l_exon)] = paste(ti[i], k, seq_along(exon), sep = "_")
				strand[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$strand, l_exon)
				gt[n+seq_len(l_exon)] = rep(.self$gtf[[i]]$type, l_exon)
				n = n + l_exon
			}
			
			if(i %% 500 == 0) {
				qqcat("@{i}/@{length(gi)} genes finished\n")
			}
		}

		value = rep(0, length(id))
		
	} else if(category == "intron") {

		n_intron = sum(sapply(.self$gtf, function(x) {
			sum(sapply(x$transcript, function(tr) length(tr$exon)+1))
		}))
		qqcat("around @{n_intron} introns\n")
		chr = character(n_intron)
		start = integer(n_intron)
		end = integer(n_intron)
		id = character(n_intron)
		value = integer(n_intron)
		strand = character(n_intron)
		gt = character(n_intron)
		gi = names(.self$gtf)
		n = 0
		for(i in seq_along(gi)) {
			ti = names(.self$gtf[[i]]$transcript)
			for(k in seq_along(.self$gtf[[i]]$transcript)) {
				exon = .self$gtf[[i]]$transcript[[k]]$exon
				ir_exon = IRanges(sapply(exon, function(x) x$start), sapply(exon, function(x) x$end))
				ir_tx = IRanges(.self$gtf[[i]]$transcript[[k]]$start, .self$gtf[[i]]$transcript[[k]]$end)
				ir_intron = setdiff(ir_tx, ir_exon)
				intron = as.data.frame(ir_intron)
				
				l_intron = length(ir_intron)
				if(l_intron) {
					chr[n + seq_len(l_intron)] = rep(.self$gtf[[i]]$chr, l_intron)
					start[n + seq_len(l_intron)] = intron[, 1]
					end[n + seq_len(l_intron)] = intron[, 2]
					id[n + seq_len(l_intron)] = paste(ti[i], k, seq_len(nrow(intron)), sep = "_")
					strand[n + seq_len(l_intron)] = rep(.self$gtf[[i]]$strand, l_intron)
					gt[n+seq_len(l_intron)] = rep(.self$gtf[[i]]$type, l_intron)
					n = n + l_intron
				}
			}
			
			if(i %% 500 == 0) {
				qqcat("@{i}/@{length(gi)} genes finished\n")
			}
		}

		value = rep(0, length(id))
		
	} else if(category %in% c("3utr", "5utr")) {

		n_utr = sum(sapply(.self$gtf, function(x) {
			length(x$transcript)
		}))

		qqcat("around @{n_utr} @{category}s\n")
		chr = character(n_utr)
		start = integer(n_utr)
		end = integer(n_utr)
		id = character(n_utr)
		value = integer(n_utr)
		strand = character(n_utr)
		gt = character(n_utr)
		gi = names(.self$gtf)
		n = 0
		for(i in seq_along(gi)) {
			ti = names(.self$gtf[[i]]$transcript)
			for(k in seq_along(.self$gtf[[i]]$transcript)) {
				if(length(.self$gtf[[i]]$transcript[[k]]$CDS)) {
					
					exon = .self$gtf[[i]]$transcript[[k]]$exon
					CDS = .self$gtf[[i]]$transcript[[k]]$CDS
					ir_exon = sort(IRanges(sapply(exon, function(x) x$start), sapply(exon, function(x) x$end)))
					ir_CDS = sort(IRanges(sapply(CDS, function(x) x$start), sapply(exon, function(x) x$end)))
					df_exon = as.data.frame(ir_exon)
					
					ir_utr = setdiff(ir_exon, ir_CDS)
					utr = as.data.frame(ir_utr)
					
					flag = 0
					if(.self$gtf[[i]]$strand == "+") {
						if(category == "3utr" && utr[nrow(utr), 2] == df_exon[nrow(df_exon), 2]) {
							pos_start = utr[nrow(utr), 1]
							pos_end = utr[nrow(utr), 2]
							flag = 1
						} else if(category == "5utr" && utr[1, 1] == df_exon[1, 1]) {
							pos_start = utr[1, 1]
							pos_end = utr[1, 2]
							flag = 1
						}
					} else {
						if(category == "3utr" && utr[1, 1] == df_exon[1, 1]) {
							pos_start = utr[1, 1]
							pos_end = utr[1, 2]
							flag = 1
						} else if(category == "5utr" && utr[nrow(utr), 2] == df_exon[nrow(df_exon), 2]) {
							pos_start = utr[nrow(utr), 1]
							pos_end = utr[nrow(utr), 2]
							flag = 1
						}
					}
					
					if(flag) {
						chr[n + 1] = .self$gtf[[i]]$chr
						start[n + 1] = pos_start
						end[n + 1] = pos_end
						id[n + 1] = ti[k]
						strand[n + 1] = .self$gtf[[i]]$strand
						gt[n + 1] = .self$gtf[[i]]$type
						n = n + 1
					}
				}
			}
			
			if(i %% 500 == 0) {
				qqcat("\r")
				qqcat("@{i}/@{length(gi)} genes finished           \n")
				flush.console()
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
			ti = names(.self$gtf[[i]]$transcript)
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
				id[n + 1] = ti[k]
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
	
	chr = chr[1:n]
	start = start[1:n]
	end = end[1:n]
	id = id[1:n]
	value = value[1:n]
	strand = strand[1:n]
	gt = gt[1:n]

	start = as.integer(start) - 1  # bed file is 0-based
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
	merged <<- TRUE
	
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


GTF$methods(getValueByGeneID = function(geneID = NULL, type = c("chr", "name", "strand", "type", "end", "start")) {
	"get values corresponding to gene IDs"
	type = match.arg(type)
	if(is.null(geneID)) {
		sapply(.self$gtf, function(x) x[[type]])
	} else {
		sapply(.self$gtf[geneID], function(x) x[[type]])
	}
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

# the name for transcript will be named `gn`
union_transcript = function(gene, gn = NULL) {
	
	if(length(gene$transcript) == 1) {
		if(!is.null(gn)) {
			names(gene$transcript) = gn
		}
		return(gene)
	}
	
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
	
	s = min(sapply(gene$transcript, function(x) x$start))
	e = max(sapply(gene$transcript, function(x) x$end))
	gene2 = gene
	gene2$transcript = list(merged_transcript = list(exon = exon, CDS = CDS, start = s, end = e, type = type))
	
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
