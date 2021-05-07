se.heatmap <- function(se, samples = NULL, hm.colors = c("lightgray", "darkorange"), ...) {
	se.all <- sort(unique(unlist(se)))
	mat <- matrix(0, ncol = length(se), nrow = length(se.all))
	rownames(mat) <- se.all
	if (!is.null(samples)) {
		colnames(mat) <- samples
	} else {
		colnames(mat) <- 1:length(se)
	}

	for (i in 1:length(se)) {
		mat[,i] <- as.numeric(se.all %in% se[[i]])
	}

	gplots::heatmap.2(mat, scale = "none", col = hm.colors, trace = "none", key = FALSE, ...)
}