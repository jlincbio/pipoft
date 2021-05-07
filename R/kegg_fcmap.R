kegg.fcmap <- function(fc, 
	col.neg = "red", col.pos = "green", col.zero = "white",
	col.na = "purple", append.border = TRUE, grain = 0.5, col.black = "black") {
		z.max <- max(abs(floor(min(fc, na.rm = TRUE))), abs(ceiling(max(fc, na.rm = TRUE))))
		z.range <- seq(from = -z.max, to = z.max, by = grain)
		y.cutoff <- colorRampPalette(c(col.neg, col.zero, col.pos))(length(z.range))
		y.bins <- .bincode(fc, breaks = z.range) + 1
		y.fgcol <- y.cutoff[y.bins]
		y.na <- is.na(y.fgcol)
		y.fgcol[is.na(y.fgcol)] <- gplots:::col2hex(col.zero)
		y.fg <- y.fgcol
		if (append.border) {
			y.fgborder <- ifelse(y.na, yes = gplots:::col2hex(col.na), no = gplots:::col2hex(col.black))
			y.fg <- paste(y.fgcol, y.fgborder, sep = ",")
		}
		return(y.fg)
}