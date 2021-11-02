# d = read.csv('GCF_000001405.37_GRCh38.p11_feature_table.txt', header = T, sep = '\t')
# d = d[which(d[,'seq_type'] == 'chromosome'),]
# d = d[which(as.vector(d[,'product_accession']) != ''),]
# g = d[,c('X..feature', 'product_accession', 'symbol', 'GeneID')]
# g$refseq = sapply(as.vector(g$product_accession), function(x){strsplit(x, '\\.')[[1]][1]})
# write.csv(g, 'hg38_refseq_annot.txt', quote =F)

convert = function(geneID=NULL, refID=NULL, entrzID=NULL, sp='hg'){
	if (sp == 'hg'){
		db = read.csv('/data5/home/rongbin/software/genome/refseq/hg38_refseq_annot.csv')
	}else if (sp == 'mm'){
		db = read.csv('/data5/home/rongbin/software/genome/refseq/mm10_refseq_annot.csv')
	}else{
		cat('wrong specias, it should be either hg or mm')
	}
	if (!is.null(refID) & is.null(geneID)){
		#convert refseq id to gene id
		refID = sapply(refID, function(x){strsplit(x, '\\.')[[1]][1]})
		intersect_refID = intersect(refID, as.vector(db$refseq))
		restID = refID[which(!refID %in% intersect_refID)]
		# deal with known ids
		if (length(intersect_refID) != 0){
			get = merge(db, as.matrix(intersect_refID), by.x = 'refseq', by.y = 1)
			get_vector = as.vector(get$symbol); names(get_vector) = as.vector(get$refseq)
		}else{
			get = NULL
			get_vector = NULL
		}
		#deal with beyond ids
		rest_get = data.frame('refseq' = restID)
		rest_get$X..feature = rep(NA, length(restID))
		rest_get$product_accession = rep(NA, length(restID))
		rest_get$symbol = rep(NA, length(restID))
		rest_get$GeneID = rep(NA, length(restID))
		rest_vector = as.vector(rest_get$symbol); names(rest_vector) = as.vector(rest_get$refseq)
		get = rbind(get, rest_get)
		get_vector = c(get_vector, rest_vector)
		return(list('res_df' = get, 'res_vect' = get_vector))
	}else if (!is.null(geneID) & is.null(refID)){
		intersect_geneID = intersect(geneID, as.vector(db$symbol))
		restID = refID[which(!geneID %in% intersect_geneID)]
		# deal with known ids
		if (length(intersect_geneID) != 0){
			get = merge(db, as.matrix(intersect_geneID), by.x = 'symbol', by.y = 1)
			get_vector = as.vector(get$refseq); names(get_vector) = as.vector(get$symbol)
		}else{
			get = NULL
			get_vector = NULL
		}
		#deal with beyond ids
		rest_get = data.frame('symbol' = restID)
		rest_get$X..feature = rep(NA, length(restID))
		rest_get$product_accession = rep(NA, length(restID))
		rest_get$refseq = rep(NA, length(restID))
		rest_get$GeneID = rep(NA, length(restID))
		rest_vector = as.vector(rest_get$refseq); names(rest_vector) = as.vector(rest_get$symbol)
		get = rbind(get, rest_get)
		get_vector = c(get_vector, rest_vector)
		return(list('res_df' = get, 'res_vect' = get_vector))
	}else if (!is.null(entrzID) &){

	}else{
		cat('please input one of geneID vector or refseqID vector.')
	}
}