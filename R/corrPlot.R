# internal functions, not exported
plotBlank <- function(x = NULL, cex.base = 0) {
	plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "")
	if (!is.null(text)) {
		text(5, 5, x, cex = 1.4 + cex.base)
	}
}
plotCorrs <- function(x, y, col, cex.base = 0) {
	m <- match(x$RefSeq, y$RefSeq)
	x1 <- x$logFC
	y1 <- y$logFC[m]
	rcor <- cor(x1, y1, method = "spearman")
	plot(x1, y1, xlab = "logFC", ylab = "logFC", col = col, pch = 19, cex.axis = 0.4 + cex.base)
	title(main = sprintf("%.3f", rcor), col.main = col, cex.main = 1 + cex.base, font.main = 1)
}
plotGenes <- function(x, y, col, cex.base = 0) {
	m <- match(x$Symbol, y$Symbol)
	x1 <- x$mFC
	y1 <- y$mFC[m]
	rcor <- cor(x1, y1, method = "spearman")
	plot(x1, y1, xlab = "", ylab = "", col = col, pch = 19, cex.axis = 0.4 + cex.base)
	title(main = sprintf("%.3f", rcor), col.main = col, cex.main = 1 + cex.base, font.main = 1)
}
nameTranslate <- function(x) {
	return(deparse(substitute(x)))
}

se.corrplot <- function(tables, samples = NULL, fill.both = FALSE) {
	if (!is.list(tables)) {
		stop("data (\'tables\') must be specified as a list; use tables = list(...) instead")
	}
	
	if (length(tables) > 5) {
		stop("correlation plot for n > 5 samples is not implemented.")
	}
	fc <- lapply(tables, function(x) data.table::as.data.table(x)[,list(mFC = mean(logFC)), c("Symbol")])
	tn <- length(tables)
	mat <- matrix(data = 1:tn^2, tn, tn)
	mt <- matrix(NA, tn, tn)
	t.col <- rainbow(tn)
	layout(mat)
	
	cex.base <- round(1/log(tn + 1), digits = 1)
	
	if (is.null(samples)) {
		# if no names specified, auto generate
		samples <- as.character(1:tn)
	}
	
	for (i in 1:tn) {
		for (j in 1:tn) {
			if (i == j) {
				# diagonal, print label
				mt[i,j] <- paste(i,"x",j)	
				plotBlank(samples[i], cex.base = cex.base)
			} else {
				if (!fill.both) {
					# do not fill both; print whitespace 
					if (!is.na(mt[j,i])) {
						# plotted - blank now
						plotBlank()
						next
					}
				}
			    plotGenes(fc[[i]], fc[[j]], t.col[j], cex.base = cex.base)
			}
			mt[i,j] <- t.col[j]
		}
	}
}