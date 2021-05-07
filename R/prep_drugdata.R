prep.drugdata <- function(drugs, platforms, effects = NULL) {
    load(system.file("keggCompounds.RData", package = "pipoft"))
    table.platforms <- read.table(file = platforms, header = TRUE, sep = "\t")
    table.drugs <- read.table(file = drugs, header = TRUE, sep = "\t")

    drugExp <- matrix(NA, nrow = dim(table.drugs)[1], ncol = length(pipoft.keggList))
    rownames(drugExp) <- table.drugs$Compound
    colnames(drugExp) <- pipoft.keggList

    myPlatforms <- vector(mode = "list", length = dim(table.platforms)[1])
    names(myPlatforms) <- table.platforms$Platform
    for (i in 1:length(myPlatforms)) {
        myPlatforms[[i]] <- read.table(
            file = table.platforms$File[i], header = FALSE, sep = "\t")
    }

    myDrugs <- vector(mode = "list", length = dim(table.drugs)[1])
    names(myDrugs) <- table.drugs$Compound
     for (i in 1:length(myDrugs)) {
        myDrugs[[i]] <- read.table(
            file = table.drugs$File[i], header = FALSE, sep = "\t")
			#file = table.platforms$File[i], header = FALSE, sep = "\t")

        # building index for platform
        pindex <- match(table.drugs$Platform[i], names(myPlatforms))
        p.curr <- myPlatforms[[pindex]]
		
		myDrugs[[i]]$Symbol <- p.curr[match(myDrugs[[i]][,1], p.curr[,1]),2]
		myDrugs[[i]] <- na.omit(myDrugs[[i]])

        if (table.drugs$LogTransformed == TRUE) {
            myDrugs[[i]]$FC <- myDrugs[[i]]$Treatment - myDrugs[[i]]$Control
        } else {
            # need to log transform the data to calculate FC
            myDrugs[[i]]$FC <- 
                log(myDrugs[[i]]$Treatment/myDrugs[[i]]$Control, base = 2)
        }
        
        myDrugs[[i]]$nFC <-
            (myDrugs[[i]]$FC - min(myDrugs[[i]]$FC))/
            (max(myDrugs[[i]]$FC) - min(myDrugs[[i]]$FC))

        # re-sorting gene symbols
        cur.match <- myDrugs[[i]]$FC[match(pipoft.keggList, myDrugs[[i]]$Symbol)]
        cur.match[is.na(cur.match)] <- 0
        cur.match[cur.match == -Inf] <- min(cur.match[is.finite(cur.match)])
        cur.match[cur.match ==  Inf] <- max(cur.match[is.finite(cur.match)])

        drugExp[match(names(myDrugs)[i], rownames(drugExp)),] <- cur.match
        remove(cur.match, pindex, p.curr)
    }
    
    na.kegg <- apply(drugExp, 2, function(x) sum(is.na(x)))
    na.kegg <- which(na.kegg == dim(drugExp)[1])
    
    if (length(na.kegg) > 0) {
        drugExp <- drugExp[,-na.kegg]
    }
    
    which(unname(apply(drugExp, 2, function(x) sum(is.na(x)))) != 0)
    which(unname(apply(drugExp, 2, function(x) sum(is.infinite(x)))) != 0)

    drugData <- list()
    if (!is.null(effects)) {
        load(file = effects)
        drugData$records <- data.frame(drugExp, drugEffect)
        drugData$indexExp <- 1:dim(drugExp)[2]
        drugData$indexEff <- (dim(drugExp)[2] + 1):(dim(drug.data)[2])
    } else {
        # drug effects data missing, leave blank
        drugData$records <- data.frame(drugExp)
        drugData$indexExp <- 1:dim(drugExp)[2]
        drugData$indexEff <- NULL
    }
    return(drugData)
}

    