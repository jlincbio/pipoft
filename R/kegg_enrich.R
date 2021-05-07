kegg.enrich <- function(genes, kegg, background = NULL, p.cutoff = NULL,
	pretty = FALSE, use.chisq = FALSE, trim = TRUE, save.genes = FALSE, simulate.p.value = FALSE, ...) {
	genes <- as.vector(genes)
	keggSpace <- unique(unlist(lapply(kegg, function(x) x$Genes)))
		if (is.null(background)) {
		background <- keggSpace
	} else {
		background <- as.vector(background)
	}
	if (trim) {
		background <- intersect(background, keggSpace)
	}
	message(c(
		sprintf("Number of genes in KEGG data : %d\n", length(keggSpace)),
		sprintf("Number of genes in background: %d\n", length(background)),
		sprintf("Number of genes in input list: %d\n", length(genes))))
	
	p.preset <- 0.01
	if (!is.null(p.cutoff)) {
		p.current <- as.numeric(p.cutoff)
		if (is.na(p.critical)) {
			message("Warning: invalid p-value specified; significance level is reverted to p = 0.01")
			p.current <- p.preset
		}
		p.critical <- p.current 
		p.preset <- p.current
	} else {
		p.critical <- p.preset / length(background)
	}
		
	z <-  data.frame(Accession = names(kegg), log_p = 0, Significance = "N/S",
			log_odds = 0, Enrichment = "Similar", fgIn = 0, fgOut = 0, bgIn = 0,
			bgOut = 0, Pathway = sapply(kegg, function(x) x$Name))
			
	t <- vector(mode = "list", length = length(z$Accession))
	g <- vector(mode = "list", length = length(z$Accession))
	names(g) <- z$Accession
		
	for (i in 1:length(z$Accession)) {
		a0 <- kegg[[i]]$Genes[kegg[[i]]$Genes %in% genes]
		a1 <- kegg[[i]]$Genes[kegg[[i]]$Genes %in% background]
		b0 <- genes[!(genes %in% a0)]
		b1 <- background[!(background %in% a1)]
		
		z$fgIn[i]  <- length(a0)
		z$bgIn[i]  <- length(a1)
		z$fgOut[i] <- length(b0)
		z$bgOut[i] <- length(b1)
		g[[i]] <- a0
		
		mat <- rbind(c(z$fgIn[i], z$bgIn[i]), c(z$fgOut[i], z$bgOut[i]))
		if (use.chisq) {
			t[[i]] <- chisq.test(x = mat, simulate.p.value = simulate.p.value)
		} else {
			t[[i]] <- fisher.test(x = mat, simulate.p.value = simulate.p.value)
		}
	}
	
	p <- sapply(t, function(x) x$p.value)
	z$log_p <- log(p)
	z$Significance <- ifelse((p < p.critical), yes = "Pass", 
		no = ifelse((p < p.preset), yes = "Marginal", no = "N/S"))
	
	z$log_odds <- log(((z$fgIn + 0.01)/(z$bgIn + 0.01))/((z$fgOut + 0.01)/(z$bgOut + 0.01)))
	z$Enrichment <- ifelse(z$log_odds > 0.1, yes = "Over", 
		no = ifelse(z$log_odds < -0.1, yes = "Under", no = "Similiar"))
	z <- z[order(z$log_p),] 
	if (pretty) {
		z <- z[z$Significance != "N/S",]
	}
	if (save.genes) {
		return(list(Summary = z, Genes = g[z$Accession]))
	}
	return(z)
}