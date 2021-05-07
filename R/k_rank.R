kegg.rank <- function(expr, targets, kegg = NULL, sites = NULL, prefix = "keggRank", path = getwd(), alpha = 0.05, perl = NULL) { 
	perl.path <- set.perl(perl = perl)
	perl.test <- system.file("k_rank.pl", package = "pipoft")
	
	args <- c(perl.test)
	if (file.valid(expr)) {
		args <- c(args, "--input", expr)
	}
	
	if (file.valid(kegg)) {
		args <- c(args, "--kegg", kegg)
	} else {
		stop("A KEGG gene database (.geneKegg) must be specified!")
	}
	
	if (!is.null(targets)) {
		targets <- paste(targets, collapse = ",")
		args <- c(args, "--target", targets)
	} else {
		stop("Gene(s) must be specified as primary target!")
	}
	
	if (!is.null(sites)) {
		if (file.valid(sites)) {
			args <- c(args, "--sites", sites)
		}
	}
	
	args <- c(args, "--alpha", alpha, "--path", path, "--prefix", prefix, "--rreturn")
	for (i in 1:length(args)) {
		if (grepl(" ", args[i])) {
			args[i] <- paste("\'", args[i], "\'", sep = "")
		}
	}
	
	cmd <- paste(args, collapse = " ")
	if (perl.path != "") {
		outputs <- system2(command = perl, args = args, stdout = TRUE)
	}
	
	if (file.valid(outputs)) {
		res <- read.table(file = outputs, header = FALSE, sep = "\t", comment.char = "#")
		colnames(res) <- c("Symbol", "Class", "relative.rank", "p.median", "p.target", "Significance")
		rownames(res) <- as.character(res$Symbol)
		return(res)
	}
	stop("Unknown error: cannot read from output\n\t", outputs)
}