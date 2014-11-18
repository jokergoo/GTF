
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
