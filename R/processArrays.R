prep.exparrays <- function(targets, annotations = NULL, use.ensembl = TRUE, plots = TRUE, norm.method = NULL, offset = NULL) {
	# generate run barcode
	suppressWarnings(op.date <- format(Sys.Date(), format = "%Y%m%d"))
	suppressWarnings(op.runc <- intToUtf8(sample(c(65:90, 97:122, 48:57), 6)))

	if (is.null(annotations)) {
		# add annotations to package
		annotations <- system.file("keggBuild.pl", package = "pipoft")
	}

	temp.table <- read.table(file = targets, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	fc.con <- unique(names(table(temp.table$Condition)))
	fc.def <- paste(fc.con, collapse = "-")
	fc.alt <- paste(rev(fc.con), collapse = "-")

	cat("Changes in gene expression should be defined as \'", fc.def, "\' [y/n]: ", sep = "")
	user.ctrl <- readLines(con = "stdin", n = 1)
	if (user.ctrl == ""){
		user.ctrl <- "y"
	}
	user.ctrl <- match.arg(tolower(user.ctrl), c("y", "n"))

	if (user.ctrl == "y") {
		fc.design <- fc.def
	} else {
		fc.design <- fc.alt
	}
	
	cat("Expression change defined now as:", fc.design, "\n")
	temp.target <- targets

	cat("Preparing to read raw microarray data...\n")
	targets <- limma::readTargets(file = temp.target)
	x <- limma::read.maimages(targets, source = "agilent", green.only = TRUE)

	f <- factor(targets$Condition, levels = unique(targets$Condition))
	design <- stats::model.matrix(~0 + f)
	colnames(design) <- levels(f)
	contrast <- limma::makeContrasts(contrasts = fc.design, levels = design)

	# initializing optimization routines
	if (is.null(norm.method)) {
		testMethod <- c("scale", "quantile")
	} else {
		testMethod <- match.arg(arg = tolower(norm.method), choices = c("scale", "quantile"))
	}

	if (is.null(offset)) {
		offset.max <- ceiling(max(x$Eb))
		offset.step <- round(ifelse((offset.max < 100), yes = (offset.max/12), no = (offset.max/25))) 
		testOffset <- sort(c(1, seq(from = 0, to = offset.max, by = offset.step)))
	
		if (testOffset[length(testOffset)] < offset.max) {
		testOffset <- c(testOffset, offset.max)
		}
		cat("Selecting a range of 0 to", offset.max, "in increments of", offset.step, "to test offsets...\n")
	} else {
		testOffset <- offset 
	}

	if (!is.null(norm.method) & !is.null(offset)) {
		cat("Proceeding with", testMethod, "normalization at an offset of", testOffset, "...\n")
		testNormalization <- FALSE
	} else {
		cat("Determining optimal normalization and background correction conditions...\n")
		testNormalization <- TRUE
	}

	optimize.array <- vector(mode = "list", length = length(testMethod))
		 model.fit <- vector(mode = "list", length = length(testMethod))
	  model.eBayes <- vector(mode = "list", length = length(testMethod))
	 model.dfprior <- vector(mode = "list", length = length(testMethod))

	output.options <- TRUE
	if (testNormalization) {
		# no parameters specified - check input options
		if (!plots) {
			output.options <- FALSE
		}
	}

	for (i in 1:length(testMethod)) {
		optimize.array[[i]] <- vector(mode = "list", length = length(testOffset))
	         model.fit[[i]] <- vector(mode = "list", length = length(testOffset))
	      model.eBayes[[i]] <- vector(mode = "list", length = length(testOffset))
		 model.dfprior[[i]] <- vector(mode = "numeric", length = length(testOffset))
		
		for (j in 1:length(testOffset)) {
			cat("Checking method \'", testMethod[i], "\' (offset = ", testOffset[j], ")...", sep = "")
			optimize.array[[i]][[j]] <- limma::backgroundCorrect(x, method="normexp", offset = testOffset[j], verbose = FALSE)
			optimize.array[[i]][[j]] <- limma::normalizeBetweenArrays(optimize.array[[i]][[j]], method = testMethod[i])
			optimize.array[[i]][[j]] <- limma::avereps(optimize.array[[i]][[j]], ID = optimize.array[[i]][[j]]$genes$ProbeName)
			model.fit[[i]][[j]] <- limma::lmFit(object = optimize.array[[i]][[j]], design = design)
			model.eBayes[[i]][[j]] <- limma::contrasts.fit(model.fit[[i]][[j]], contrast)
			model.eBayes[[i]][[j]] <- limma::eBayes(model.eBayes[[i]][[j]])
			model.dfprior[[i]][j] <- model.eBayes[[i]][[j]]$df.prior
		
			if (output.options) {
				output <- topTable(model.eBayes[[i]][[j]], adjust="BH", coef= fc.design,
				   				   genelist = model.eBayes[[i]][[j]]$genes, number = Inf)
				if (!testNormalization) {
					# save output
					output.final <- output
				}
				outputFile.txt <- paste(
					fc.design, "_", op.date, "-", testMethod[i], testOffset[j], ".txt", sep = "")
				outputFile.RData <- paste(
					fc.design, "_", op.date, "-", testMethod[i], testOffset[j], ".RData", sep = "")
		
				dfPrior <- round(model.dfprior[[i]][[j]], digits = 3)			
				subText <- paste(testMethod[i], " normalization (offset ", testOffset[j], "); dfPrior = ", dfPrior, sep = "")
		
				write.table(output, file = outputFile.txt, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
				save(output, file = outputFile.RData)
				pdf(file = gsub(".txt", ".pdf", outputFile.txt), width = 8, height = 6)
				limma::plotSA(model.fit[[i]][[j]], sub = subText)
				invisible(dev.off())
				remove(output)
			}
			cat("\n")
		}
		if (output.options) {
			pdf(file = paste("dfprior_", testMethod[i], "-", op.date, ".pdf", sep = ""), width = 8, height = 6)
			plot(testOffset, model.dfprior[[i]], xlab = "Offset", ylab = "df.prior", main = testMethod[i], type = "b", pch = 19)
			dev.off()
		}
	}

	max.dfprior.value <- sapply(model.dfprior, max)
	max.dfprior.pos <- sapply(1:length(testMethod), function(x) (which(model.dfprior[[x]] == max.dfprior.value[x])))

	select.method <- which(max.dfprior.value == max(max.dfprior.value))
	select.offset <- max.dfprior.pos[select.method]

	if (testNormalization) {
		cat("Comparing prior variances suggests that over the range of [", min(testOffset), "-", max(testOffset), "],\n", sep = "")
		cat(testMethod[select.method], "normalization appears to be more optimal.\n")
		cat("The best offset for this normalization method is:", testOffset[select.offset],"\n")
		cat("Preparing output based on these parameters...\n")
		output.final <- limma::topTable(
			model.eBayes[[select.method]][[select.offset]], adjust = "BH", coef = fc.design,
			genelist = model.eBayes[[select.method]][[select.offset]]$genes, number = Inf)
		if (plots) {
			# regenerate expression table
			outputFile.txt <- paste(fc.design, "_", op.date, "-", 
				testMethod[select.method], testOffset[select.offset], "_summary.txt", sep = "")
			write.table(output.final, file = outputFile.txt, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
		}
	}

	cat("Annotating microarray results...\n")
	g.annot <- read.table(file = annotations, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	fname.output <- paste(fc.design, "_", op.date, "-", testMethod[select.method], testOffset[select.offset], "_sumAnnotated.txt", sep = "")
	fname.sixcol <- paste(fc.design, "_", op.date, "-", testMethod[select.method], testOffset[select.offset], "_psiAnnotated.txt", sep = "")
	if (length(unique(names(table(g.annot[,5])) %in% c("1", "2"))) == 1) {
		# proceed only if annotation file can be understood
		if (!use.ensembl) {
			# prune non-refSeq entries
			g.annot <- g.annot[(g.annot[,5] == 2),]
			output.final <- output.final[output.final$SystematicName %in% g.annot[,4],]
		}
	
		i.annot <- match(output.final$SystematicName, g.annot[,4])
		h.annot <- data.frame(output.final$SystematicName, output.final$logFC, g.annot[i.annot,c(1,2,3,7,4)])
		colnames(h.annot) <- c("Accession", "logFC", "Chrom", "Start", "End", "Symbol", "Annotation")
		h.annot$Symbol <- ifelse(test = is.na(h.annot$Symbol), yes = h.annot$Accession, no = h.annot$Symbol)
		h.annot <- h.annot[,c(3:6, 1,2)]
		colnames(h.annot) <- c("Chrom", "Start", "End", "Symbol", "RefSeq", "logFC")
		h.annot <- h.annot[!is.na(h.annot[,1]),]
		h.annot <- h.annot[order(h.annot$logFC, decreasing = TRUE),]
		write.table(x = h.annot, file = fname.sixcol, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
	
		output.final.colnames <- colnames(output.final)
		output.final <- data.frame(output.final, g.annot[i.annot,7])
		colnames(output.final) <- c(output.final.colnames, "Symbol")
		output.final <- output.final[order(output.final$logFC, decreasing = TRUE),]
		write.table(x = output.final, file = fname.output, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
			
		cat(sprintf("BED6 output: %s\n", basename(fname.sixcol)))
		cat(sprintf("Full output: %s\n", basename(fname.output)))	
	} else {
		warning("Error parsing annotation file; leaving results as is.\n")
	}
}