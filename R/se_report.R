se.report <- function(results, file = "Rplots.pdf", ...) {
	m <- diag(dim(results$test)[1])
	pdf(file = file, ...)
	for (i in 1:dim(results$test)[1]) {
		for (j in 1:dim(results$test)[1]) {
			if (m[i,j] < 1) {
				# plot has not been generated
				x.name <- rownames(results$test)[i]
				y.name <- rownames(results$test)[j]

				plot(x = results$test[i,], y = results$test[j,], col = "darkgreen",
					xlab = x.name, ylab = y.name, pch = 16,
					main = sprintf("Correlation of %s and %s Gene Expressions", x.name, y.name),
					sub = "Points: Candidate Genes in Side Effects Prediction Pool")
				plot(x = results$predict[,i], y = results$predict[,j], col = "darkred",
					xlab = x.name, ylab = y.name, pch = 16,
					main = "Correlation of Predicted Side Effects",
					sub = sprintf("Points: Predicted S/E from %s\\/%s\ Expr. Profiles", x.name, y.name)) 
				m[i,j] <- m[j,i] <- 1 # skip existing plots
			}
		}
	}
	
	err.indiv <- sapply(1:dim(results$predict)[2], function(x) which(results$predict[,x] > 0.5))
	err.common <- Reduce(intersect, err.indiv)
	barplot(results$err[err.common], 
		las = 2, cex.names = 0.6, col = "violet",
		main = "Common Side Effects and Respective Error",
		ylab = "Prediction Error")
	dev.off()
}