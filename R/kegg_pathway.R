kegg.pathway <- function(x, pathway, ...) {
	run.stats <- run.perl(cmd = "keggProcess.pl", perl = NULL, pathway)
	output.cmap <- paste("colorMap_", pathway, ".txt", sep = "")
	output.gkey <- paste("keggKey_", pathway, ".txt", sep = "")
	if (file.exists(output.cmap) & (file.exists(output.gkey))) {
		# array data: Symbol and logFC
		x0 <- read.table(x, sep = "\t", header = TRUE)
		keggKey <- read.table(output.gkey, header = TRUE, comment.char = "")
		keggMap <- read.table(output.cmap, header = TRUE, comment.char = "")
		
		p.symbols <- keggKey$Symbol[match(keggMap$X.hsa, keggKey$X.hsa)]
		p.fc <- x0$logFC[match(p.symbols, x0$Symbol)]
		keggMap$FG <- kegg.fcmap(p.fc, ...)
		colnames(keggMap)[1] <- "#hsa"
		write.table(keggMap, file = output.cmap, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
		return(output.cmap)
	}
	stop("Error encountered during processing KEGG data for: ", pathway)
	return(NA)
}

kegg.import <- function(pathway, genes, omit.species = TRUE) {
	if (!file.exists(pathway)) {
		stop("Error: KEGG pathway information data file has not been specified!")
	}
	if (!file.exists(genes)) {
		stop("Error: KEGG pathway information data file has not been specified!")
	}
	
	f.intermediate <- tempfile()
	if (omit.species) {
		pretty <- 2
	} else {
		pretty <- 0
	}
	
	run.stats <- run.perl(cmd = "keggToList.pl", perl = NULL, pathway, genes, f.intermediate, pretty)
	if (run.stats > 0) {
		stop("Errors encountered during data import!")
	}
	if (file.exists(f.intermediate)) {
		keggTable <- read.table(file = f.intermediate, header = FALSE, sep = "\t")
		keggList <- vector(mode = "list", length = dim(keggTable)[1])
		names(keggList) <- keggTable$V1
		for (i in 1:dim(keggTable)[1]) {
			keggList[[i]]$Name <- keggTable$V3[i]
			keggList[[i]]$Genes <- unlist(strsplit(keggTable$V2[i], "::"))
		}
		return(keggList)
	}
	stop("Error: data processing was complete but could not be imported into workspace!")
}