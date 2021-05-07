se.predict <- function(training, effects, expr, output = NULL, ...) {	
	drug.data <- data.frame(training, effects)
	index.exp <- 1:dim(training)[2]
	index.eff <- (dim(training)[2] + 1):(dim(drug.data)[2])
	
	drug.rft <- vector(mode = "list", length = length(index.eff))
	drug.rff <- vector(mode = "list", length = length(index.eff))
	for (i in 1:length(index.eff)) {
		cat("Side Effect:", colnames(effects)[i], "\n")
		drug.rff[[i]] <- as.formula(paste(
			paste(colnames(drug.data)[index.eff[i]], collapse = " + "), " ~ ",
			paste(colnames(drug.data)[index.exp], collapse = " + ")))
		drug.rft[[i]] <- ranger::ranger(formula = drug.rff[[i]], data = drug.data, importance = "impurity", ...)
	}

	# prediction errors determined by ranger models
	drug.oob <- sapply(drug.rft, function(x) x$prediction.error)
	names(drug.oob) <- colnames(effects)
	mat.importance <- sapply(drug.rft, function(x) ranger::importance(x, expr.test))
	sum.importance <- apply(mat.importance, 1, sum)
	
	expr <- expr[,which(colnames(expr) %in% colnames(training))]
	expr.order <- match(colnames(training), colnames(expr))
	expr.test <- expr[,expr.order]
	colnames(expr.test) <- colnames(training)
	expr.test[,which(is.na(expr.order))] <- 0

	expr.predict <- matrix(NA, nrow = length(drug.rft), ncol = dim(expr.test)[1])
	for (i in 1:length(drug.rft)) {
		expr.predict[i,] <- ranger::predict(drug.rft[[i]], expr.test)$predictions
	}
	
	rfResults <- list()
	rfResults$trees <- drug.rft
	rfResults$formula <- drug.rff
	rfResults$err <- drug.oob
	rfResults$importance <- mat.importance
	rfResults$test <- expr.test
	rfResults$predict <- expr.predict
	
	if (!is.null(output)) {
		save(rfResults, file = output)
	}
	return(rfResults)
}