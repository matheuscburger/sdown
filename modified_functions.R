# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: down
# Date......: 25/05/2012
# Author(s).: matheus
###############################################################################


###TESTE###
#q <- qFilt.2
#design = NULL
#contrasts
#sort = TRUE
#stringent = FALSE
#NA.threshold = 0.6                                                                                                                                                                        
#ndups = 1
#spacing = NULL
#dupcor
###########

my.limmaCtData <- function (q, design = NULL, contrasts, sort = TRUE, stringent = FALSE, NA.threshold = NA,                                                                                                                                                                        
		ndups = 1, spacing = NULL, dupcor, ...)                                                                                                                                                                                                  
{
	if( (!is.na(NA.threshold)) && stringent==TRUE){ ## Opcoes mutuamente exclusivas: NA.tresshold OU stringent
		warning("Using stringent=FALSE")
		stringent = FALSE
	}
	data <- exprs(q)
	featPos <- featurePos(q)
	if (missing(dupcor)) {
		if (ndups > 1) {
			dup.cor <- duplicateCorrelation(data, ndups = ndups,
					spacing = spacing, design = design)
			temp <- unwrapdups(featPos, ndups = ndups, spacing = spacing)
			featPos <- apply(temp, 1, paste, collapse = ";")
		}
		else {
			dup.cor <- NULL
		}
	}
	fit <- lmFit(data, design = design, ndups = ndups, spacing = spacing,
			correlation = dup.cor$consensus, ...)
	if (!missing(contrasts))
		fit <- contrasts.fit(fit, contrasts = contrasts)
	fit2 <- eBayes(fit)
	out <- list()
	if (!missing(contrasts)) {
		coefs <- colnames(contrasts)
		cont <- design %*% contrasts
	}
	else {
		coefs <- colnames(design)
		cont <- design
	}
	for (coef in coefs) {
		res <- topTable(fit2, coef = coef, number = nrow(fit2), sort.by = "none", ...)
		both.means <- both.cats <- array(0, c(nrow(res), 2),
				list(res$ID, c("Test", "Reference")))
		for (i in c(-1, 1)) {
			sample <- ifelse(i == 1, "Test", "Reference")
			index <- cont[, coef] == i
			mean <- rowMeans(unwrapdups(data[, index], ndups = ndups,
							spacing = spacing))
			new.cat <- rep("OK", length(mean))
			old.cat <- unwrapdups(featureCategory(q)[, index],
					ndups = ndups, spacing = spacing)
			count.cat <- apply(old.cat, 1, function(x) sum(x %in%
										c("Undetermined", "Unreliable")))
			cutoff <- ifelse(stringent, 1, floor(sum(index)*NA.threshold))
			new.cat[count.cat >= cutoff] <- "Undetermined"
			both.means[, sample] <- mean
			both.cats[, sample] <- new.cat
		}
		res.out <- cbind(res$ID, featPos, res[, c("t", "P.Value",
								"adj.P.Val", "logFC")], 2^(-res$logFC), both.means,
				both.cats)
		colnames(res.out) <- c("genes", "feature.pos", "t.test",
				"p.value", "adj.p.value", "ddCt", "FC", "meanTarget",
				"meanCalibrator", "categoryTarget", "categoryCalibrator")
		if (sort)
			res.out <- res.out[order(res.out$adj.p.value), ]
		out[[coef]] <- res.out
	}
	res <- decideTests(fit2, ...)
	rownames(res) <- topTable(fit2, sort = "none", n = nrow(fit2))$ID
	out[["Summary"]] <- res
	out
}



## TESTES ##
#q <- qFilt.2
#groups = levels(files$Treatment)[files$Treatment]
#calibrator =  "Control"
#alternative = "two.sided"
#paired = FALSE
#replicates = TRUE
#sort = TRUE
#stringent = TRUE
#p.adjust = "fdr"
############

my.ttestCtData <- function (q, groups = NULL, calibrator, alternative = "two.sided", paired = FALSE, replicates = TRUE, sort = TRUE, stringent = TRUE, p.adjust = "BH", ...) 
{
	data <- exprs(q)
	if (replicates) {
		split.by <- rownames(data)
	}else{
		split.by <- featurePos(q)
	}
	data2 <- split(as.data.frame(data), split.by)
	feats <- split(featureNames(q), split.by)
	featPos <- split(featurePos(q), split.by)
	if (length(groups) != ncol(data)) 
		stop("Dimensions of data and groups doesn't match\n")
	groups <- factor(groups, levels = unique(groups))
	if (length(levels(groups)) != 2) 
		stop("Two factor levels required for 'groups'\n")
	if (missing(calibrator)) 
		calibrator <- groups[1]
	g1 <- groups == calibrator
	g2 <- groups != calibrator
	
	
	###TESTES###
	#x <- data2[[1]]
	############
	t.tests <- lapply(data2, function(x) {
				x <- as.matrix(x)
				if (all(x == x[1, 1], na.rm=T)) {
					list(p.value = 1, statistic = NA, estimate = c(x[1, 1], x[1, 1]))
				}else {
					res <- t.test(x[, g1], x[, g2], alternative = alternative, paired = paired, ...)
					res[["estimate"]] <- c(mean(x[, g1]), mean(x[, g2]))
					res
				}
			})
	means <- t(sapply(t.tests, "[[", "estimate"))
	colnames(means) <- c("meanCalibrator", "meanTarget")
	p.value <- sapply(t.tests, "[[", "p.value")
	t.value <- sapply(t.tests, "[[", "statistic")
	genes <- sapply(feats, "[[", 1)
	featurePos <- sapply(featPos, paste, collapse = ";")
	adj.p.value <- p.adjust(p.value, method = p.adjust)
	cal <- grep(calibrator, colnames(means))
	FC <- means[, "meanTarget"] - means[, "meanCalibrator"]
	FC2 <- 2^(-FC)
	out <- data.frame(genes, featurePos, t.value, p.value, adj.p.value, 
			FC, FC2, means, row.names = 1:length(genes))
	for (l in unique(groups)) {
		new.cat <- rep("OK", length(data2))
		old.cat <- split(featureCategory(q[, groups == l]), split.by)
		count.cat <- sapply(old.cat, function(x) sum(unlist(x) %in% 
									c("Undetermined", "Unreliable")))
		cutoff <- ifelse(stringent, 1, ceiling(sum(groups == 
										l)/2))
		new.cat[count.cat >= cutoff] <- "Undetermined"
		out[, paste("category", ifelse(l == calibrator, "Calibrator", 
								"Target"), sep = "")] <- new.cat
	}
	names(out) <- c("genes", "feature.pos", "t.test", "p.value", 
			"adj.p.value", "ddCt", "FC", colnames(means), grep("category", 
					colnames(out), value = TRUE))
	if (sort) 
		out <- out[order(out$p.value), ]
	out
}


###TESTES###
#q <- qFilt.2
#groups = files$Treatment 
#calibrator = "Control" 
#p.adjust="fdr" 
#stringent=TRUE
#replicates = TRUE
#alternative = "two.sided"
#paired = FALSE
############

my.mannwhitneyCtData <- function (q, groups = NULL, calibrator, alternative = "two.sided", paired = FALSE, replicates = TRUE, sort = TRUE, stringent = TRUE, p.adjust = "BH", ...) 
{
	data <- exprs(q)
	if (replicates) {
		split.by <- rownames(data)
	}else {
		split.by <- featurePos(q)
	}
	data2 <- split(as.data.frame(data), split.by)
	feats <- split(featureNames(q), split.by)
	featPos <- split(featurePos(q), split.by)
	if (length(groups) != ncol(data)) 
		stop("Dimensions of data and groups doesn't match\n")
	groups <- factor(groups, levels = unique(groups))
	if (length(levels(groups)) != 2) 
		stop("Two factor levels required for 'groups'\n")
	if (missing(calibrator)) 
		calibrator <- groups[1]
	g1 <- groups == calibrator
	g2 <- groups != calibrator
	mw.tests <- lapply(data2, function(x) {
				x <- as.matrix(x)
				if(length(na.exclude(x[,g1])) < 2 || length(na.exclude(x[,g2])) < 2){
					res <- list()
					res[["p.value"]] <- NA
					res[["estimate"]] <- c(mean(x[, g1]), mean(x[, g2]))
				}else{
					res <- wilcox.test(x[, g1], x[, g2], alternative = alternative,	paired = paired, exact = FALSE)
					res[["estimate"]] <- c(mean(x[, g1]), mean(x[, g2]))
				}
				res
			})
	means <- t(sapply(mw.tests, "[[", "estimate"))
	colnames(means) <- c("meanCalibrator", "meanTarget")
	p.value <- unlist(sapply(mw.tests, "[[", "p.value"))
	mw.value <- sapply(mw.tests, function(x) paste(names(x[["statistic"]]), 
						"=", x[["statistic"]]))
	genes <- sapply(feats, "[[", 1)
	featurePos <- sapply(featPos, paste, collapse = ";")
	adj.p.value <- p.adjust(p.value, method = p.adjust)
	cal <- grep(calibrator, colnames(means))
	FC <- means[, "meanTarget"] - means[, "meanCalibrator"]
	FC2 <- 2^(-FC)
	out <- data.frame(genes, featurePos, mw.value, p.value, adj.p.value, 
			FC, FC2, means, row.names = 1:length(genes))
	for (l in unique(groups)) {
		new.cat <- rep("OK", length(data2))
		old.cat <- split(featureCategory(q[, groups == l]), split.by)
		count.cat <- sapply(old.cat, function(x) sum(unlist(x) %in% 
									c("Undetermined", "Unreliable")))
		cutoff <- ifelse(stringent, 1, ceiling(sum(groups == 
										l)/2))
		new.cat[count.cat >= cutoff] <- "Undetermined"
		out[, paste("category", ifelse(l == calibrator, "Calibrator", 
								"Target"), sep = "")] <- new.cat
	}
	names(out) <- c("genes", "feature.pos", "MW.value", "p.value", 
			"adj.p.value", "ddCt", "FC", colnames(means), grep("category", 
					colnames(out), value = TRUE))
	if (sort) 
		out <- out[order(out$p.value), ]
	out
}

		
