file.valid <- function(file) {
	if (is.null(file)) {
		return(FALSE)
	}
	if (file.exists(file)) {
		return(TRUE)
	}
	stop("File not found: ", file)
}

tableToTemp <- function(x) {
	file.out <- tempfile(pattern = "pipoft_")
	write.table(x, file = file.out, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
	return(file.out)
}

run.perl <- function(cmd, perl = NULL, ...) {
	if (is.null(perl)) {
		perl.path <- set.perl(perl = NULL)
	}
	cmd.path <- system.file(cmd, package = "pipoft")
	run.stats <- system2(command = perl.path, args = c(cmd.path, ...))
	return(run.stats)
}

set.perl <- function(perl, silent = TRUE) {
	if (missing(perl) || is.null(perl)) {
		perl <- "perl"
	}
	
	perl <- Sys.which(perl)
	if (perl == "" || perl == "perl") {
		stop("Perl executable cannot be located; specify the location with \'perl = \'")
	}
	
	if (.Platform$OS == "windows") {
		if (length(grep("rtools", tolower(perl))) > 0) {
			perl.ftype <- shell("ftype perl", intern = TRUE)
			if (length(grep("^perl=", perl.ftype)) > 0) {
				perl <- sub('^perl="([^"]*)".*', "\\1", perl.ftype)
			}
		}
	}
	
	if (!silent) {
		cat("Path to Perl binary: ", perl, "\n")
	}
	return(perl)
}

prep.setup <- function(perl = NULL, cpan = NULL, silent = TRUE) {
	perl.path <- set.perl(perl = perl)
	perl.test <- system.file("p_tests.pl", package = "pipoft")
	perl.args <- c(perl.test)
	if (perl.path != "") {
		# no perl; although this is probably caught with stop() anyways
		exit.stats <- system2(command = perl.path, args = c(perl.args))
	}
	
	perl.stats <- FALSE
	if (exit.stats > 0) {
		# perl packages missing
		stop("Missing Perl modules found; check configuration and try again.\n")
	} else {
		# create a log file to suppress warning messages at startup
		perl.logger <- paste(dirname(perl.test), '.perlpath', sep = '/')
		cat(perl.path, file = perl.logger)
		perl.stats <- TRUE
	}

	if (!silent) {
		if (perl.stats) {
			cat("No missing Perl modules were detected.\n")
		}
		return(!as.logical(exit.stats))
	}
}