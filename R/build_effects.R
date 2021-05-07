build.effects <- function(sider.id, sider.effects, kegg.filter = NULL,
	drugbank = NULL, drugbank.id = NULL, 
	drugbank.filter = (if (is.null(drugbank)) NULL else "inhibit"),
	cutoff = -1, path = getwd(), prefix = "pipoft", perl = NULL) {
	
	perl.path <- set.perl(perl = perl)
	perl.test <- system.file("build_mat.pl", package = "pipoft")
	
	args <- c(perl.test)
	if (file.valid(sider.id)) {
		args <- c(args, "--drugnames", sider.id)
	}
	
	if (file.valid(sider.effects)) {
		args <- c(args, "--effects", sider.effects)
	}
	
	if (file.valid(kegg.filter)) {
		# check if drugbank.id is specified
		if (!is.null(drugbank.id)) {
			if (file.valid(drugbank.id) && file.valid(kegg.filter)) {
				args <- c(args, "--genelist", kegg.filter, "--uniprot", drugbank.id)
			}
		} else {
			stop("drugbank.id is required when kegg.filter is specified.")
		}
	}
	if (file.valid(drugbank)) {
		args <- c(args, "--drugbank", drugbank)
	}
	
	if (!is.null(drugbank.filter)) {
		filters <- paste(drugbank.filter, collapse = ",")
		args <- c(args, "--filters", filters)
	}
	
	if (cutoff > 0) {
		args <- c(args, "--cutoff", cutoff)
	}
	
	args <- c(args, "--path", path, "--prefix", prefix, "--rreturn")
	for (i in 1:length(args)) {
		if (grepl(" ", args[i])) {
			args[i] <- paste("\'", args[i], "\'", sep = "")
		}
	}
	
	cmd <- paste(args, collapse = " ")
	
	if (perl.path != "") {
		outputs <- system2(command = perl, args = args, stdout = TRUE)
	}
	outputs <- unlist(strsplit(outputs, split = ","))
	
	mat <- as.matrix(read.csv(file = outputs[1], header = FALSE))
	drugList <- scan(file = outputs[2], what = character(), sep = "\n")
	seList <- scan(file = outputs[3], what = character(), sep = "\n")
	rownames(mat) <- drugList
	colnames(mat) <- seList
	return(mat)
}