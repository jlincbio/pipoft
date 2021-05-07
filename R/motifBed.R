prep.motif <- function(fasta, motif, dest = NULL, perl = "perl") {
    cmd <- system.file("fastaMotif.pl", package = "pipoft")
	if (is.null(dest)) {
		bed.output <- paste(tempfile(), "bed", sep = ".")
	} else {
		bed.output <- dest
	}
	
	run.stats <- system2(command = perl, args = c(cmd, fasta, motif, bed.output))
	if (run.stats == 0) {
		res <- read.table(bed.output, sep = "\t", header = FALSE)
		if (is.null(dest)) {
			# output not specified; remove temp BED
			if (!file.remove(bed.output)) {
				message("Error removing temp file:\n", bed.output)
			}
		}
		return(res)
	}
}