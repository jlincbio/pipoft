prep.kegg <- function(prefix = "pipoft", threads = 1, organism = "hsa", temp = tempdir()) {	
	run.stats <- run.perl(cmd = "keggBuild.pl", perl = NULL, 
		prefix, organism, threads, temp) # success = 0

	res <- list(genes = NA, pathways = NA, summary = NA)
	if (run.stats == 0) {
		res$genes <- paste(prefix, '.geneKegg', sep = '')
		res$pathways <- paste(prefix, '.pathKegg', sep = '')
		res$summary <- paste(prefix, '.txt', sep = '')
	} else {
		stop("unable to build KEGG dataset (code ", run.stats, ")")
	}
	return(res)
}