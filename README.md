# *pipoft*
R package for estimating genes affected by off-target binding of pyrrole-imidazole polyamides and side effect estimation.

~~Mar 30, 2019, Version 0.0.0.1110~~

May  7, 2021, Version 0.0.0.1200

## System Requirements
  1. A Perl interpreter (version 5.16 or above)
  2. A system running macOS or Linux; WSL may also be compatible
  3. R (version 3.3 or above recommended)
  
### Perl Dependencies
  * `Parallel::ForkManager` (but multithreading is not required)
  * `Statistics::TTest`
  * `File::Type`
  * `PerlIO::gzip`
  * `XML::Simple`
  * `XML::Twig`

Missing Perl packages may be installed with CPAN, including certain core modules missing in a minimal Perl installation. On Linux distributions, `PerlIO::gzip` may require additional zlib libraries & headers.

### R Dependencies
  * `limma`
  * `ranger`
  * `data.table`
  * `optparse`
  * `gplots`
  * `png`
  * `grid`

## Installation
__Note: from this point on all code blocks are R codes__
```
library("devtools") # install.packages("devtools")
install_github("jlincbio/pipoft")
```

Alternatively, from the release tarball:

`R CMD INSTALL pipoft-0.0.0.1200.gz`

Installation will fail if there are missing Perl packages; in this case, please run CPAN to install the missing dependencies prior to re-installation.

## Manuals

### Caveats
Initial release is not object-oriented but implementation is planned at a later time. This package provides an interactive front-end for the Perl scripts used in the publication; nonetheless, the package may be revamped in native R in the future.

### Configurations
When `prep.setup()` is run without arguments, the function will search for a Perl interpreter and output a list of missing Perl modules; in the current version, however, upon package installation this is automatically invoked so there is no need to run this unless you wish to change the Perl interpreter.

### Search Genomic Binding Sites for a PI Polyamide by Motif 
`prep.motif` will do this for you; use `dest` to specify an filename if you want to save the results to a 4-column BED file.   
```
myBed <- prep.motif(hg19.fa, motif = "TGWWGGCGW", dest = "~/tgwwggcgw.bed")
```

### Pathway Analysis using KEGG
__Updated May 2021__
`pipoft` also provides a basic interface to characterize gene set enrichment by pathway overrepresentation, using a list containing the names of differentially expressed genes as the input. For example, RNA-seq data can be processed using a mixture of `limma` and `edgeR` to generate a candidate gene list for downstream analysis.

`kegg.enrich` will perform a statistical overrepresentation test when a vector of gene symbols is entered as input. Statistical significance is then evaluated by either Fisher's Exact Test (default) or chi-square. By default the significance level is set to be *Np* = 0.01, where *N* is the number of genes in the background list.
```
# optional: use `kegg.import()` to load the generated KEGG pathway data into R workspace
keggData <- kegg.import("hsa_pathway", "hsa_geneList")
resEnrichment <- kegg.enrich(myGenes, kegg = keggData, background = backgroundGenes) # if "background" is not specified, defaults to all genes in KEGG data
```

### Side Effect Prediction and Report
Use `se.predict()` and `se.report()` to generate a tree and make expression-linked plots as PDF.  
```
# use all default settings to grow the forest
myTrees <- se.predict(myTrainingSet, myDrugBank, myArrays, "myTrees.RData")

# other ranger-specific parameters can also be passed directly
myTrees <- se.predict(myTrainingSet, myDrugBank, myArrays, num.trees = 1000, num.random.splits = 2)

se.report(myTrees, "arrayTrees.pdf")
se.report(myTrees, "arrayTrees-A4.pdf", paper = "a4r")
```

### Building drug-side effect training set
You will need to first create a tab-delimited version that includes compound names, URL, GEO platform links, and whether the data need to be log-transformed. Upon the creation of this file, pipoft will load the expression data to generate a matrix object ready for side effect prediction. Prior to that, you will need to download the files onto the hard drive and run a separate function call (this will be interactive) to perform preliminary data processing. We used previously published gene lists to limit the search scope, and considered for the purpose of discussion only drugs with inhibitory actions. In the end 705 side effects were obtained as of August 2018 from parsing available data from SIDER (version 4.1) and DrugBank (April 2017).

## Afterthoughts
This is a wrapper package that tries to integrate the Perl code used for the manuscript together with the R code and therefore depends externally on Perl and some Perl modules. Due to licensing issues for DrugBank and KEGG, precompiled reference datasets cannot be distributed. Please acquire licenses first if you do not qualify for free-use.

## Citations
Lin J, Krishnamurthy S, Yoda H, Shinozaki Y, Watanabe T, Koshikawa N, Takatori A, Horton P, Nagase H. Estimating genome-wide off-target effects for pyrrole-imidazole polyamide binding by a pathway-based expression profiling approach. *PLoS ONE* **14**(4): e0215247, 2019. doi: 10.1371/journal.pone.0215247

In R:
`citation("pipoft")`
