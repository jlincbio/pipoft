library(biomaRt)
library(Rsubread)
library(limma)
library(edgeR)
library(gplots)
library(pipoft)
library(data.table)

if (file.exists("seqCounts_P29mtP29-20200622.RData")) {
	load("seqCounts_P29mtP29-20200622.RData")
} else {
	bam <- list.files(path = getwd(), full.names = TRUE, pattern = "bam")
	bam <- bam[grep(x = bam, pattern = "bai$", invert = TRUE)]
	cell.types <- c("P29mtP29","P29mtB82M")
	targets <- data.frame(list(Type = cell.types, BAM = bam))
	fc <- featureCounts(targets$BAM, annot.inbuilt = "mm10", nthreads = 16)
	x <- DGEList(counts = fc$counts, genes = fc$annotation[,c("GeneID","Length")])
	# uncorrected x; should be fine? corr = 1 comp. to voom
	fx <- log(cpm(x) + 0.5, base = 2)
	save.image(file = "seqCounts_P29mtP29-20200622.RData")
}

pMean <- function(x, lookup) {
	y <- which(lookup$Symbol %in% x)
	return(mean(lookup$logFC[y], na.rm = TRUE)/length(which(!is.na(lookup$logFC[y]))))
}
load("mmu_geneList-kegg.RData")

# entrezgene apparently outdated
mm10.ensembl <- useMart("ensembl", "mmusculus_gene_ensembl")
fx.annot <- getBM(c("external_gene_name", "entrezgene_id"), "entrezgene_id", rownames(fx), mm10.ensembl)

fx.rmsd <- sqrt((fx[,2] - fx[,1])^2)
fx.symbols <- fx.annot$external_gene_name[match(rownames(fx), fx.annot$entrezgene_id)]
fx.logFC <- as.vector(scale(fx[,2] - fx[,1]))
	#fx.table$P29mtB82M - fx.table$P29mtP29))
fx.table <- na.omit(data.frame(fx.symbols,fx, fx.logFC, fx.rmsd))
colnames(fx.table) <- c("Symbol", cell.types, "logFC", "RMSD")

fx.lm <- lm(fx.table$P29mtB82M ~ fx.table$P29mtP29)
fx.ht <- as.data.table(fx.table[abs(fx.table$logFC) >= 1.5,])
fx.ht <- data.frame(fx.ht[,list(xP29mtP29 = mean(P29mtP29), xP29mtB82M = mean(P29mtB82M)), c("Symbol")])
rownames(fx.ht) <- fx.ht$Symbol
fx.ht <- as.matrix(fx.ht[,2:3])
colnames(fx.ht) <- gsub("P29mt", "", cell.types)

cat(rownames(fx.ht), sep = "\n", file = "diffGenes_p29_fc1.5-20200622.txt")
if (!file.exists("diffGenes_p29_fc1.5_summary-20200622.txt")) {
	system("keggEnrich --input diffGenes_p29_fc1.5-20200622.txt --pvalue 0.01 --database mmu_keggDb-20200620 --pretty") # no sig
	file.rename("output_062220_summary.txt", "diffGenes_p29_fc1.5_summary-20200622.txt")
	file.rename("output_062220_geneList.txt", "diffGenes_p29_fc1.5_geneList-20200622.txt")
}

write.table(
	fx.table[order(fx.table$logFC, decreasing = TRUE),],
	file = "testReport_expressions-20200622.txt", quote = FALSE,
	sep = "\t", col.names = TRUE, row.names = FALSE)
pdf("testReport_revised-20200622.pdf", paper = "a4")
x.breaks <- seq(from = floor(min(fx.table$RMSD)), to = ceiling(max(fx.table$RMSD)), by = 0.5)
x.palette <- rev(colorRampPalette(c("darkblue", "blue", "lightblue"))(length(x.breaks) - 1))
x.col <- x.palette[.bincode(fx.table$RMSD, breaks = x.breaks, include.lowest = TRUE, right = FALSE)]

for (j in 1:2) {
	plot(fx.table$P29mtP29, fx.table$P29mtB82M, xlab = cell.types[1], ylab = cell.types[2],
		pch = 19, main = "Expression Correlation", col = x.col,
		sub = expression("Values: log"[2]*italic(" Merged CPM")))

	abline(fx.lm, col = "red") # R2 = .8986
	legend("topleft", 
		legend = c(
			expression(italic("R")^2*" = 0.899"),
			sprintf("RMSD <= %.1f", x.breaks[2]),
			sprintf("RMSD <= %.1f", x.breaks[ceiling((2 + length(x.breaks))/2)]),
			sprintf("RMSD <= %.1f", x.breaks[length(x.breaks)])
			),
		col = c("red", 
			x.palette[2],
			x.palette[ceiling((2 + length(x.breaks))/2)],
			x.palette[length(x.palette)]),
		pch = c(95, 19, 19, 19))
	if (j > 1) {
		# additional labeling
		lab.idx <- which(fx.table$Symbol == "Slc16a3")
		points(fx.table$P29mtP29[lab.idx], fx.table$P29mtB82M[lab.idx], col = "red", pch = 19)
		text(fx.table$P29mtP29[lab.idx], fx.table$P29mtB82M[lab.idx], col = "red", labels = "Slc16a3", pos = 4)
	}
}

hist(fx.table$logFC,
	main = "Fold-Change Distribution",
	xlab = expression("Values: log"[2]*italic(" FC")),
	sub = "FC = P29mtB82M/P29mtP29")
hist(fx.table$logFC[abs(fx.table$logFC) >= 1.5], col = "darkgreen",
	main = expression("Fold-Change Distribution (Cutoff: log"[2]*italic(" FC")*" >= ±1.5)"),
	xlab = expression("Values: log"[2]*italic(" FC")),
	sub = "FC = P29mtB82M/P29mtP29")

y.breaks <- seq(from = floor(min(fx.ht)), to = ceiling(max(fx.ht)), by = 1)
y.palette <- unique(c(colorRampPalette(c("black", "darkgreen"))(3), 
	colorRampPalette(c("darkgreen", "green"))(n = length(y.breaks) - 3)))
heatmap.2(as.matrix(fx.ht), scale = "none", key = TRUE, trace = "none", 
	col = y.palette, breaks = y.breaks, cexCol = 1.2,
	keysize = 1.5, key.title = "", key.xlab = expression("log"[2]*italic(" CPM")),
	density.info = "none", symbreaks = FALSE, symkey = FALSE,
	main = expression("Genes with log"[2]*italic(" FC")*" >= ±1.5"))

keggList.fc <- sapply(keggList, 
	function(x) mean(fx.table$logFC[which(fx.table$Symbol %in% x)], na.rm = TRUE))
#keggList.fcn <- sapply(keggList, pMean, lookup = fx.table)	
keggList.fc.top <- sort(keggList.fc, decreasing = TRUE)[1:30]
barplot(keggList.fc, main = "Mean Fold Changes by KEGG Pathways", axisnames = FALSE)
barplot(keggList.fc.top, las = 2, cex.names= .7, col = "darkgreen", main = "Top 30 Most Upregulated KEGG Pathways")
dev.off()

for (i in 1:length(keggList.fc.top)) {
	cur.id <- names(keggList.fc.top)[i]
	if (!file.exists(sprintf("%s.png", cur.id))) {
		system(sprintf("keggPsi %s", cur.id))
		cur.match <- match(keggList[[match(cur.id, names(keggList))]], fx.table$Symbol)
		cur.fc <- fx.table$logFC[cur.match]
		names(cur.fc) <- fx.table$Symbol[cur.match]
		cur.cl <- kegg.fcmap(cur.fc)
		cur.map <- read.table(sprintf("keggKey_%s.txt", cur.id), header = FALSE, sep = "\t")
		cur.tab <- data.frame(cur.map$V1[match(names(cur.fc), cur.map$V2)], cur.cl)
		write.table(cur.tab, file = sprintf("colorMap_%s.txt", cur.id), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
		kegg.colors(pathway = cur.id, colors = sprintf("colorMap_%s.txt", cur.id))
		readLines("stdin", n = 1)
		remove(cur.id, cur.fc, cur.cl, cur.map, cur.tab)
	}
}
save(fx.table, fx.annot, keggList.fc, file = "repCounts_P29mtP29-20200622.RData")

# 6/25
keggList.fc.low <- sort(keggList.fc)[1:30]
pdf("downReg_pathways-lowest30.pdf", paper = "a4")
barplot(keggList.fc.low, col = "darkblue", las = 2, cex.names=.7, main = "Top 30 Most Downregulated KEGG Pathways")
dev.off()
for (i in 1:length(keggList.fc.low)) {
	cur.id <- names(keggList.fc.low)[i]
	if (!file.exists(sprintf("%s.png", cur.id))) {
		system(sprintf("keggPsi %s", cur.id))
		cur.match <- match(keggList[[match(cur.id, names(keggList))]], fx.table$Symbol)
		cur.fc <- fx.table$logFC[cur.match]
		names(cur.fc) <- fx.table$Symbol[cur.match]
		cur.cl <- kegg.fcmap(cur.fc)
		cur.map <- read.table(sprintf("keggKey_%s.txt", cur.id), header = FALSE, sep = "\t")
		cur.tab <- data.frame(cur.map$V1[match(names(cur.fc), cur.map$V2)], cur.cl)
		write.table(cur.tab, file = sprintf("colorMap_%s.txt", cur.id), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
		kegg.colors(pathway = cur.id, colors = sprintf("colorMap_%s.txt", cur.id))
		readLines("stdin", n = 1)
		remove(cur.id, cur.fc, cur.cl, cur.map, cur.tab)
	}
}

# BAM files were processed into counts per million (CPM) by Subread [R1] based on the mm10 reference mouse genome annotations. Fold-change values were defined as the difference in log2 CPM for a given transcript between P29mtB82M and P29mtP29 such that a log FC > 0 indicates a higher CPM in P29mtB82M compared to the same transcript in P29mtP29. For the subsequent pathway analysis, results from individual transcripts were then combined by gene symbols using BioMart [R2] by querying Entrez GeneID and external gene names. The root-square deviation (RSD) for each transcript is defined as the square root of the square of the difference between P29mtB82M and P29mtP29, and is used as a metric to determine the relative extent of difference in magnitude of gene expression. Gene set pathway enrichment analysis were made against all pathways under the mouse ("mmu") category in KEGG [R3] and statistical significance was evaluated by Fisher's Exact Test under the significance level p < 0.01. KEGG pathway maps were made via API calls to the KEGG web interface (last accessed June 2020) via the R package "pipoft" [R4] Mean expression changes per KEGG pathway were calculated from the mean logFC for matching symbols within a particular KEGG pathway.

#[R1] Y. Liao, G.K. Smyth and W. Shi. "The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote." Nucleic Acids Research, 41(10):e108, 2013
#[R2] S. Durinck, Y. Moreau, A. Kasprzyk, S. Davis, B. De Moor, A. Brazma, W. Huber. "BioMart and Bioconductor: A Powerful Link Between Biological Databases and Microarray Data Analysis." Bioinformatics 21(16): 3439-3440, 2005.
#[R3] M. Kanehisa, S. Goto. "KEGG: Kyoto Encyclopedia of Genes and Genomes." Nucleic Acids Res. 28(1): 27-30, 2000.
#[R4] J. Lin, S. Krishnamurthy, H. Yoda, Y. Shinozaki, T. Watanabe, N. Koshikawa, A. Takatori, P. Horton, H. Nagase. "Estimating genome-wide off-target effects for pyrrole-imidazole polyamide binding by a pathway-based expression profiling approach." PLoS ONE 14(4): e0215247, 2019.