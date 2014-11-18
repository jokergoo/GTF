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

GTF$methods(read = function(gtf) {
	"read from GTF file"
	
	#perlScript = qq("@{system.file(package = 'GTF')}/extdata/gtf_to_json.pl")
	perlScript = "D:\\personal_project\\GTF\\inst\\extdata\\gtf_to_json.pl"
	tmp = tempfile()
	cmd = qq("perl \"@{perlScript}\" --input \"@{normalizePath(gtf)}\" --output \"@{tmp}\"")
	error = try(system(cmd))
	if(class(error) == "try-error") {
		file.remove(tmp)
		stop(error, appendLF = FALSE)
	}
	gtf <<- fromJSON(file = tmp)
	file.remove(tmp)
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

GTF$methods(availableTypes = function(type = c("gene", "transcript")) {
	type = match.arg(type)[1]
	if(type == "gene") {
		unique(sapply(.self$gtf, function(x) x$type))
	} else {
		unique(unlist(lapply(.self$gtf, function(gene) unique(sapply(gene$transcript, function(tr) tr$type)))))
	}
})

GTF$methods(genes = function(type = NULL) {
	
	chr = sapply(.self$gtf, function(x) x$chr)
	start = sapply(.self$gtf, function(x) x$start)
	end = sapply(.self$gtf, function(x) x$end)
	gene_id = names(.self$gtf)
	gene_name = sapply(.self$gtf, function(x) x$name)
	strand = sapply(.self$gtf, function(x) x$strand)
	gene_type = sapply(.self$gtf, function(x) x$type)
	gr = GRanges(seqnames = Rle(chr),
	        ranges = IRanges(start = start,
			                 end = end),
			strand = Rle(strand),
			gene_id = gene_id, 
			gene_type = gene_type, 
			gene_name = gene_name)
	names(gr) = gene_id
	
	if(!is.null(type)) {
		l = gt %in% type
		gr = gr[l]
	}
	return(gr)
})

GTF$methods(transcripts = function(type = NULL) {

	n_tx = sum(sapply(.self$gtf, function(x) {
		length(x$transcript)
	}))

	chr = character(n_tx)
	start = integer(n_tx)
	end = integer(n_tx)
	id = character(n_tx)
	value = integer(n_tx)
	strand = character(n_tx)
	tt = character(n_tx)
	gi = character(n_tx)
	gn = character(n_tx)
	
	all_gi = names(.self$gtf)
	all_gn = sapply(.self$gtf, function(x) x$name)
	n = 0
	for(i in seq_along(all_gi)) {
		tx = .self$gtf[[i]]$transcript
		if(!is.null(type)) {
			l = sapply(tx, function(x) x$type) %in% type
			tx = tx[l]
		}
		k = length(tx)
		
		chr[n + seq_len(k)] = rep(.self$gtf[[i]]$chr, k)
		start[n + seq_len(k)] = sapply(tx, function(x) x$start)
		end[n + seq_len(k)] = sapply(tx, function(x) x$end)
		id[n + seq_len(k)] = names(tx)
		strand[n + seq_len(k)] = rep(.self$gtf[[i]]$strand, k)
		tt[n + seq_len(k)] = sapply(tx, function(x) x$type)
		gi[n + seq_len(k)] = rep(all_gi[i], k)
		gn[n + seq_len(k)] = rep(all_gn[i], k)
		
		n = n + k
	}
	
	gr = GRanges(seqnames = Rle(chr),
	        ranges = IRanges(start = start,
			                 end = end),
			strand = Rle(strand),
			transcript_id = id, 
			transcript_type = tt, 
			gene_id = gi,
			gene_name = gn)
	names(gr) = id
	return(gr)
})

GTF$methods(exons = function() {
	
	n_tx = sum(sapply(.self$gtf, function(x) {
		length(x$transcript)
	}))
	
	n_exon = sum(sapply(.self$gtf, function(x) {
		sum(sapply(x$transcript, function(tr) length(tr$exon)))
	}))
		
	chr = character(n_exon)
	start = integer(n_exon)
	end = integer(n_exon)
	id = character(n_exon)
	value = integer(n_exon)
	strand = character(n_exon)
	
	gi = character(n_exon)
	gn = character(n_exon)
	
	ti = character(n_exon)
	tt = character(n_exon)
	
	all_gi = names(.self$gtf)
	all_gn = sapply(.self$gtf, function(x) x$name)
	
	n = 0
	for(i in seq_along(all_gi)) {
		tx = .self$gtf[[i]]$transcript
		tii = names(tx)
		
		for(k in seq_along(tx)) {
			exon = tx[[k]]$exon
			l_exon = length(exon)
			if(l_exon) {
				chr[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$chr, l_exon)
				start[n + seq_len(l_exon)] = sapply(exon, function(x) x$start)
				end[n + seq_len(l_exon)] = sapply(exon, function(x) x$end)
				if(.self$gtf[[i]]$strand == "+") {
					id[n + seq_len(l_exon)] = paste(tii[k], k, seq_along(exon), sep = "_")
				} else {
					id[n + seq_len(l_exon)] = paste(tii[k], k, rev(seq_along(exon)), sep = "_")
				}
				strand[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$strand, l_exon)
				ti[n+seq_len(l_exon)] = rep(tii[k], l_exon)
				tt[n+seq_len(l_exon)] = rep(tx[[k]]$type, l_exon)
				gi[n+seq_len(l_exon)] = rep(all_gi[i], l_exon)
				gn[n+seq_len(l_exon)] = rep(all_gn[i], l_exon)
				
				n = n + l_exon
			}
		}
	}
	
	gr = GRanges(seqnames = Rle(chr),
	        ranges = IRanges(start = start,
			                 end = end),
			strand = Rle(strand),
			exon_id = id, 
			transcript_id = ti,
			transcript_type = tt,
			gene_id = gi,
			gene_name = gn)
	names(gr) = id
	return(gr)
})

GTF$methods(CDS = function() {
	
	n_tx = sum(sapply(.self$gtf, function(x) {
		length(x$transcript)
	}))
	
	n_exon = sum(sapply(.self$gtf, function(x) {
		sum(sapply(x$transcript, function(tr) length(tr$exon)))
	}))
		
	chr = character(n_exon)
	start = integer(n_exon)
	end = integer(n_exon)
	id = character(n_exon)
	value = integer(n_exon)
	strand = character(n_exon)
	
	gi = character(n_exon)
	gn = character(n_exon)
	
	ti = character(n_exon)
	tt = character(n_exon)
	
	all_gi = names(.self$gtf)
	all_gn = sapply(.self$gtf, function(x) x$name)
	
	n = 0
	for(i in seq_along(all_gi)) {
		tx = .self$gtf[[i]]$transcript
		tii = names(tx)
		
		for(k in seq_along(tx)) {
			exon = tx[[k]]$CDS
			l_exon = length(exon)
			if(l_exon) {
				chr[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$chr, l_exon)
				start[n + seq_len(l_exon)] = sapply(exon, function(x) x$start)
				end[n + seq_len(l_exon)] = sapply(exon, function(x) x$end)
				if(.self$gtf[[i]]$strand == "+") {
					id[n + seq_len(l_exon)] = paste(tii[k], k, seq_len(l_exon), sep = "_")
				} else {
					id[n + seq_len(l_exon)] = paste(tii[k], k, rev(seq_len(l_exon)), sep = "_")
				}
				strand[n + seq_len(l_exon)] = rep(.self$gtf[[i]]$strand, l_exon)
				ti[n+seq_len(l_exon)] = rep(tii[k], l_exon)
				tt[n+seq_len(l_exon)] = rep(tx[[k]]$type, l_exon)
				gi[n+seq_len(l_exon)] = rep(all_gi[i], l_exon)
				gn[n+seq_len(l_exon)] = rep(all_gn[i], l_exon)
				
				n = n + l_exon
			}
		}
	}
	
	gr = GRanges(seqnames = Rle(chr),
	        ranges = IRanges(start = start,
			                 end = end),
			strand = Rle(strand),
			CDS_id = id, 
			transcript_id = ti,
			transcript_type = tt,
			gene_id = gi,
			gene_name = gn)
	names(gr) = id
	return(gr)
})

# introns = transcript - exons
GTF$methods(introns = function() {
	
	n_tx = sum(sapply(.self$gtf, function(x) {
		length(x$transcript)
	}))
	
	n_exon = sum(sapply(.self$gtf, function(x) {
		sum(sapply(x$transcript, function(tr) length(tr$exon)))
	}))
	
	n_intron = n_tx + n_exon
		
	chr = character(n_intron)
	start = integer(n_intron)
	end = integer(n_intron)
	id = character(n_intron)
	value = integer(n_intron)
	strand = character(n_intron)
	
	gi = character(n_intron)
	gn = character(n_intron)
	
	ti = character(n_intron)
	tt = character(n_intron)
	
	all_gi = names(.self$gtf)
	all_gn = sapply(.self$gtf, function(x) x$name)
	
	n = 0
	for(i in seq_along(all_gi)) {
		tx = .self$gtf[[i]]$transcript
		tii = names(tx)
		
		for(k in seq_along(tx)) {
			tx_ir = IRangs(start = tx[[k]]$start, end = tx[[k]]$end)
			exon = tx[[k]]$exon
			l_exon = length(exon)
			if(l_exon) {
				exon_ir = IRanges(start = sapply(exon, function(x) x$start),
				                  end = sapply(exon, function(x) x$end))
				intron_ir = setdiff(tx_ir, exon_ir)
				
			} else {
				intron_ir = tx_ir
			}
			
			l_intron = length(intron_ir)
			if(l_intron) {
				chr[n + seq_len(l_intron)] = rep(.self$gtf[[i]]$chr, l_exon)
				start[n + seq_len(l_intron)] = start(intron_ir)
				end[n + seq_len(l_intron)] = end(intron_ir)
				if(.self$gtf[[i]]$strand == "+") {
					id[n + seq_len(l_intron)] = paste(tii[k], k, seq_len(l_intron), sep = "_")
				} else {
					id[n + seq_len(l_intron)] = paste(tii[k], k, rev(seq_len(l_intron)), sep = "_")
				}
				strand[n + seq_len(l_intron)] = rep(.self$gtf[[i]]$strand, l_exon)
				ti[n+seq_len(l_intron)] = rep(tii[k], l_intron)
				tt[n+seq_len(l_intron)] = rep(tx[[k]]$type, l_intron)
				gi[n+seq_len(l_intron)] = rep(all_gi[i], l_intron)
				gn[n+seq_len(l_intron)] = rep(all_gn[i], l_intron)
				
				n = n + l_intron
			}
		}
	}
	
	gr = GRanges(seqnames = Rle(chr),
	        ranges = IRanges(start = start,
			                 end = end),
			strand = Rle(strand),
			intron_id = id, 
			transcript_id = ti,
			transcript_type = tt,
			gene_id = gi,
			gene_name = gn)
	names(gr) = id
	return(gr)
})

# exon - CDS, last one
GTF$methods(threeUTRs = function() {

})

# exon - CDS, first one
GTF$methods(fiveUTRs = function() {

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
