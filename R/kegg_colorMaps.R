kegg.colors <- function(pathway, colors) {
	run.stats <- run.perl(cmd = "keggColorMaps.pl", perl = NULL, 
		pathway, colors)
	if (run.stats == 0) {
		kegg.output <- paste(pathway, "png", sep = ".")
		if (file.exists(kegg.output)) {
			message("Output: ", getwd(), "/", kegg.output)
			message("Attempting to display pathway...")
			kegg.img <- png:::readPNG(kegg.output)
			grid:::grid.raster(kegg.img)
			return(kegg.output)
		}
	} else {
		if (run.stats == 2) {
			message("Warning: file cannot be saved locally. Use URL above to open.")
			return(NULL)
		}	
	}
	return(NA)
}