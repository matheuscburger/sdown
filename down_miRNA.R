# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: down
# Date......: 10/10/2011
# Author(s).: daniel
###############################################################################


library('HTqPCR')
library('genefilter')
library(AgiMicroRna)

dir_base = "/work/projects/bit_data"
projeto = "down"
analises = "htqpcr/miRNA"

mywd <- paste(dir_base,projeto,analises,sep="/") 

setwd(mywd)
getwd()

#files <- read.delim(file.path(mywd, "targets.txt"))
#raw <- readCtData(files = files$File, path = mywd, n.features=768, feature=2, type=3, Ct=4, header=TRUE)
#sampleNames(raw) <- gsub('raw/AB/', '', sampleNames(raw))


filesA <- read.delim(file.path(mywd, "targetsA384.txt"))
rawA <- readCtData(files = filesA$File, path = mywd, n.features=384, feature=2, type=3, Ct=4, header=TRUE)
sampleNames(rawA) <- gsub('raw/A/', '', sampleNames(rawA))
filesA$Names <- gsub('raw/A/|.txt', '', filesA$File)
featureClass(rawA)	<- filesA$Treatment
show(rawA)

filesB <- read.delim(file.path(mywd, "targetsB384.txt"))
rawB <- readCtData(files = filesB$File, path = mywd, n.features=384, feature=2, type=3, Ct=4, header=TRUE)
sampleNames(rawB) <- gsub('raw/B/', '', sampleNames(rawB))
filesB$Names <- gsub('raw/B/|.txt', '', filesB$File)
featureClass(rawB)	<- filesB$Treatment
show(rawB)

#g <- featureNames(raw)[1:10]
#plotCtOverview(rawA, genes = g, xlim = c(0, 50), groups = filesA$Treatment,
#				conf.int = TRUE, ylim = c(0, 55))

#plotCtOverview(rawA, genes = g, xlim = c(0, 50), groups = filesA$Treatment,
#				calibrator = "Control")

		
#plotCtCard(raw,card=3, col.range = c(10, 35), well.size = 2.6)

#
#cv.value <- function (x, ...) 
#{
#	sdx <- sd(x, ...)
#	if (is.na(sdx) || (mean(x, ...) == 0)) 
#		return(NA)
#	val <- sdx/abs(mean(x, ...))
#	return(val)
#}
#
#cv.fun <- function(x, cv.threshold = 0.1, ...) {
#	cv <- cv.value(x,...)
#	if (is.na(cv)) {
#		return('Flagged');
#	}
#	else if (cv < cv.threshold) {
#		return('Passed');
#	}
#	else {
#		return('Failed')
#	}
#}

#
#
#erawA <- exprs(rawA)
##hist(as.numeric(apply(erawA[,c('CC10A-1','CC10A-2')], 1, cv.value, na.rm=TRUE)),na.rm=TRUE)
##summary(as.numeric(apply(erawA[,c('CC10A-1','CC10A-2')], 1, cv.value, na.rm=TRUE)))
#
#flagA.df <- matrix(nrow=dim(erawA)[1],ncol=dim(erawA)[2])
#colnames(flagA.df) <- colnames(erawA)
#flagA.df <- as.data.frame(flagA.df)
#
#for (r in levels(as.factor(filesA$Rep))) {
#	#print(r)
#	print(paste(r,filesA$Name[which(filesA$Rep==r)],collapse=';',sep=":"))
#	#print(summary(as.numeric(apply(erawA[,filesA$Name[which(filesA$Rep==r)]], 1, cv.value))))
#	#print(dim(erawA[,filesA$Name[which(filesA$Rep==r)]]))
#	cv.result <- as.character(unlist(apply(erawA[,filesA$Name[which(filesA$Rep==r)]], 1, cv.fun, na.rm=TRUE)))
#	
#	for (f in filesA$Name[which(filesA$Rep==r)]) {
#		flagA.df[[f]] <- cv.result
#	}
#}
#flag(rawA) <- flagA.df
#
#erawB <- exprs(rawB)
##hist(as.numeric(apply(erawA[,c('CC10A-1','CC10A-2')], 1, cv.value, na.rm=TRUE)),na.rm=TRUE)
##summary(as.numeric(apply(erawA[,c('CC10A-1','CC10A-2')], 1, cv.value, na.rm=TRUE)))
#
#flagB.df <- matrix(nrow=dim(erawB)[1],ncol=dim(erawB)[2])
##rownames(flagB.df) <- rownames(erawB)
#colnames(flagB.df) <- colnames(erawB)
#flagB.df <- as.data.frame(flagB.df)
#
#for (r in levels(as.factor(filesB$Rep))) {
#	#print(r)
#	print(paste(r,filesB$Name[which(filesB$Rep==r)],collapse=';',sep=":"))
#	#print(summary(as.numeric(apply(erawB[,filesB$Name[which(filesB$Rep==r)]], 1, cv.value))))
#	#print(dim(erawB[,filesB$Name[which(filesB$Rep==r)]]))
#	cv.result <- as.character(unlist(apply(erawB[,filesB$Name[which(filesB$Rep==r)]], 1, cv.fun, na.rm=TRUE)))
#	
#	for (f in filesB$Name[which(filesB$Rep==r)]) {
#		flagB.df[[f]] <- cv.result
#	}
#}
#flag(rawB) <- flagB.df
#
#show(rawB)
#
#colnames(flag(rawA))
#colnames(exprs(rawA))
#colnames(flag(rawB))
#colnames(exprs(rawB))
#
#
#
#my.setCategory <- function (q, Ct.max = 35, Ct.min = 10, replicates = TRUE, quantile = 0.9, 
#		groups, flag = TRUE, flag.out = "Failed", verbose = TRUE, 
#		plot = FALSE, ...) 
#{
#	out <- q
#	data <- exprs(q)
#	featureCategory(out)[data > Ct.max] <- "Undetermined"
#	featureCategory(out)[data < Ct.min] <- "Unreliable"
#	if (verbose) {
#		feats <- featureCategory(out)
#		cats <- sort(unique(unlist(feats)))
#		count <- array(0, c(length(cats), n.samples(out)), list(cats, 
#						sampleNames(out)))
#		for (i in 1:ncol(feats)) {
#			tab <- table(feats[, i])
#			count[names(tab), i] <- tab
#		}
#		cat("Categories after Ct.max and Ct.min filtering:\n")
#		print(count)
#	}
#	if (flag) {
#		flags <- flag(q)
#		featureCategory(out)[flags == flag.out] <- "Unreliable"
#	}
#	if (replicates) {
#		split.by <- rownames(data)
#	}
#	else {
#		split.by <- paste(rownames(data), "_well", 1:nrow(data), 
#				sep = "")
#	}
#	data2 <- split(as.data.frame(data), split.by)
#	if (!is.null(quantile)) {
#		if (missing(groups)) {
#			warning("Sample groups must be supplied to filter by standard deviation.")
#			invisible(out)
#		}
#		l.groups <- unique(groups)
#		SD <- AV <- N <- array(0, dim = c(length(data2), length(l.groups)), 
#				dimnames = list(names(data2), l.groups))
#		for (g in l.groups) {
#			SD[, g] <- sapply(data2, function(d) sd(unlist(d[, 
#												groups == g])))
#			AV[, g] <- sapply(data2, function(d) mean(unlist(d[, 
#												groups == g])))
#			N[, g] <- sapply(data2, function(d) length(unlist(d[, 
#												groups == g])))
#		}
#		if (plot) {
#			for (g in seq_along(l.groups)) {
#				dev.new()
#				par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
#				hist(SD[, g], n = 25, main = l.groups[g], xlab = "Standard deviation")
#				plot(AV[, g], SD[, g], pch = 20, xlab = "Mean", 
#						ylab = "Standard deviation", main = l.groups[g], 
#						...)
#			}
#		}
#		for (d in seq_along(data2)) {
#			gene <- names(data2)[d]
#			for (g in l.groups) {
#				sd <- SD[gene, g]
#				n <- N[gene, g]
#				av <- AV[gene, g]
#				bounds <- 0.5 + quantile/2
#				conf <- vector("numeric", 2)
#				conf[1] <- qnorm(bounds, mean = av, sd = sd)
#				conf[2] <- qnorm(1 - bounds, mean = av, sd = sd)
#				Ct <- data[split.by == gene, groups == g]
#				index <- Ct > conf[1] | Ct < conf[2]
#				featureCategory(out)[split.by == gene, groups == 
#								g][index] <- "Unreliable"
#				featureCategory(out)[data > Ct.max] <- "Undetermined"
#			}
#		}
#		if (verbose & !is.null(groups)) {
#			feats <- featureCategory(out)
#			cats <- sort(unique(unlist(feats)))
#			count <- array(0, c(length(cats), n.samples(out)), 
#					list(cats, sampleNames(out)))
#			for (i in 1:ncol(feats)) {
#				tab <- table(feats[, i])
#				count[names(tab), i] <- tab
#			}
#			cat("Categories after standard deviation filtering:\n")
#			print(count)
#		}
#	}
#	if (nrow(getCtHistory(out)) == 0) 
#		out@history <- data.frame(history = "Manually created qPCRset object.", 
#				stringsAsFactors = FALSE)
#	out@history <- rbind(out@history, capture.output(match.call(setCategory)))
#	invisible(out)
#}

#############
##### "Seta" caracteristicas - Unreliable, Undetermined, Ok


## Grupo A
# plotDensityMicroRna(exprs(rawA), "rawA") # Distribuição BiModal
# quantile entre replicas
par(ask=T)
sink(paste("results/setCategoryA_groups_Rep_Ctmax_32_quantile_0.75.txt",sep=""))
rawA.cat <- setCategory(	rawA, 
						Ct.max = 32, 
						groups = filesA$Rep, 
						quantile = 0.75, 
						replicates=TRUE, 
						verbose=TRUE,
						flag=FALSE,
						plot=TRUE)

sink()


## Grupo B
# plotDensityMicroRna(exprs(rawB), "rawB") # Distribuição BiModal 
# quantile entre replicas
sink(paste("results/setCategoryB_groups_Rep_Ctmax_32_quantile_0.75.txt",sep=""))				
rawB.cat <- setCategory(	rawB, 
						Ct.max = 32, 
						groups = filesB$Rep, 
						quantile = 0.75,
						replicates=TRUE, 
						verbose=TRUE,
						flag=FALSE,						
						plot=TRUE)
sink()

## Grupo A
# quantile entre condicoes
sink(paste("results/setCategoryA_groups_Rep_Treatment_Ctmax_32_quantile_0.90.txt",sep=""))
rawA.cat.2 <- setCategory(	rawA.cat, 
						Ct.max = 32, 
						groups = filesA$Treatment,
						quantile = 0.90, 
						replicates=TRUE, 
						verbose=TRUE,
						flag=FALSE,
						plot=TRUE)
sink()

## Grupo B
# quantile entre condicoes
sink(paste("results/setCategoryB_groups_Rep_Treatment_Ctmax_32_quantile_0.90.txt",sep=""))
rawB.cat.2 <- setCategory(	rawB.cat, 
						Ct.max = 32, 
						groups = filesB$Treatment, 
						quantile = 0.90,
						replicates=TRUE, 
						verbose=TRUE,
						flag=FALSE,						
						plot=TRUE)
sink()				
				
########
#### Plot Categories
				
system(paste("mkdir -p results/plots/ 2>&1",sep=""), intern=TRUE)
	
png(filename=paste("results/plots/plotCtCategoryA_groups_Rep.png",sep=""), 
						bg="white", res=300, width=3000, height=3000)				
plotCtCategory(rawA.cat)
dev.off()

png(filename=paste("results/plots/plotCtCategoryB_groups_Rep.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)				
plotCtCategory(rawB.cat)
dev.off()


png(filename=paste("results/plots/plotCtCategoryA_groups_Rep_Treatment.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)				
plotCtCategory(rawA.cat.2)
dev.off()

png(filename=paste("results/plots/plotCtCategoryB_groups_Rep_Treatment.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)				
plotCtCategory(rawB.cat.2)
dev.off()

#########
### Filtro - Undetermined

rawA.filt <- filterCategory(rawA.cat.2, na.categories = c("Undetermined"))
#rawA.filt

rawB.filt <- filterCategory(rawB.cat.2, na.categories = c("Undetermined"))
#rawB.filt


###rawA.filt <- rawA.cat.2
###rawB.filt <- rawB.cat.2

dim(exprs(rawA.filt))
dim(exprs(rawB.filt))

png(filename=paste("results/plots/Clustering_CtCategoryA.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)
plotCtCategory(rawA.filt, by.feature = TRUE, cexRow = 0.1)
dev.off()

png(filename=paste("results/plots/Clustering_CtCategoryB.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)
plotCtCategory(rawB.filt, by.feature = TRUE, cexRow = 0.1)
dev.off()

#########
#### Normalização

sink(paste("results/normalizeCtDataA_deltaCt.txt",sep=""))
dA.norm <- normalizeCtData(	rawA.filt, 
							norm = "deltaCt", 
							deltaCt.genes = c("MammU6-4395470", "RNU44-4373384", "RNU48-4373383"),
							#deltaCt.genes = c("MammU6-4395470")
							)
sink()

sink(paste("results/normalizeCtDataB_deltaCt.txt",sep=""))										
dB.norm <- normalizeCtData(	rawB.filt, 
							norm = "deltaCt", 
							deltaCt.genes = c("U6 snRNA-001973", "RNU44-001094", "RNU48-001006")
							#deltaCt.genes = c("U6 snRNA-001973")
							)
sink()


#x11();plot(exprs(rawA.filt), 
#	exprs(dA.norm), 
#	pch = 20, main = "DeltaCt Normalisation - A", 
#	col = rep(brewer.pal(dim(exprs(rawA.filt))[2], "Spectral"), each = dim(exprs(rawA.filt))[1])); abline(0,1)


#plot(exprs(rawB.filt), 
#		exprs(dB.norm), 
#		pch = 20, main = "DeltaCt Normalisation - B", 
#		col = rep(brewer.pal(dim(exprs(rawB.filt))[2], "Spectral"), each = dim(exprs(rawB.filt))[1]))	


#########
#### Filtro
sink(paste("results/filterCtDataA_ncat_",as.integer(0.75*length(filesA$Name)),".txt",sep=""))
qFiltA <- filterCtData(	dA.norm, 
						remove.type = "Endogenous Control", 
						remove.name = c("MammU6-4395470", "RNU44-4373384", "RNU48-4373383"),
						#remove.category = c("Undetermined","Unreliable"),
						#remove.category = c("Unreliable"),
						#remove.category = c("Undetermined"),
						n.category = 0
						#n.category = as.integer(0.75*length(filesA$Name))				
						)
sink()

sink(paste("results/filterCtDataB_ncat_",as.integer(0.75*length(filesB$Name)),".txt",sep=""))						
qFiltB <- filterCtData(	dB.norm, 
		remove.type = "Endogenous Control", 
		remove.name = c("MammU6-4395470", "RNU44-4373384", "RNU48-4373383"),
		#remove.category = c("Undetermined","Unreliable"),
		#remove.category = c("Unreliable"),
		#remove.category = c("Undetermined"),
		n.category = 0
		#n.category = as.integer(0.75*length(filesB$Name))				
)
sink()


summary(exprs(qFiltA)) ; boxplot(exprs(qFiltA)) 							
summary(exprs(qFiltB)) ; boxplot(exprs(qFiltB))						
		
eqFiltA <- exprs(qFiltA)
eqFiltB <- exprs(qFiltB)


####### "Seta" 40 em NA (Unreliables)  "Saturação"
eqFiltA[eqFiltA %in% NA] <- 40
eqFiltB[eqFiltB %in% NA] <- 40

summary(eqFiltA)
summary(eqFiltB)



png(filename=paste("results/plots/IQRA.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)

iqrA.values <- apply(eqFiltA, 1, IQR, na.rm=TRUE)
hist(iqrA.values, n = 20, main = "", xlab = "IQR across samples - A")
abline(v = median(iqrA.values, na.rm=TRUE), col = 2)
dev.off()

png(filename=paste("results/plots/IQRB.png",sep=""), 
		bg="white", res=300, width=3000, height=3000)
iqrB.values <- apply(eqFiltB, 1, IQR, na.rm=TRUE)
hist(iqrB.values, n = 20, main = "", xlab = "IQR across samples - B")
abline(v = median(iqrB.values, na.rm=TRUE), col = 2)
dev.off()

sink(paste("results/filterCtDataA_rmIQR_",median(iqrB.values, na.rm=TRUE),".txt",sep=""))						
qFiltA.2 <- filterCtData(qFiltA, remove.IQR = median(iqrA.values, na.rm=TRUE))
sink()

sink(paste("results/filterCtDataB_rmIQR_",median(iqrB.values, na.rm=TRUE),".txt",sep=""))						
qFiltB.2 <- filterCtData(qFiltB, remove.IQR = median(iqrB.values, na.rm=TRUE))
sink()

eqFiltA.2 <- exprs(qFiltA.2)
eqFiltB.2 <- exprs(qFiltB.2)

dim(eqFiltA.2)

eqFiltA.2[eqFiltA.2 %in% NA] <- 40
eqFiltB.2[eqFiltB.2 %in% NA] <- 40

eqFiltA.2 <- eqFiltA 
eqFiltB.2 <- eqFiltB 
dim(eqFiltA.2)
dim(eqFiltB.2)

# heatmap
library('gplots')
library('cluster')

hclust.ave <- function (d,...) {
	return(hclust(d=d,method="average",...))
}


hclust.diana <- function (x,...) {
	return(as.hclust(diana(x=x,diss=TRUE,...)))
}


hv1<-heatmap.2(
		as.matrix(eqFiltA.2),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		hclustfun=hclust.diana,
		key=TRUE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=T,
		Colv=T,
		cexRow=0.45,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "Group A"
)


hv1<-heatmap.2(
		as.matrix(eqFiltB.2),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		hclustfun=hclust.diana,
		key=TRUE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		breaks=100,
		Rowv=T,
		Colv=T,
		cexRow=0.45,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "Group B"
)

eqFiltA.2 <- exprs(qFiltA.2)
eqFiltB.2 <- exprs(qFiltB.2)

########## SAVE
source("/work/projects/manalysis/myArray/myArrayFunctions.R")	
saveArray('saves/')

write.table(eqFiltA.2, file="saves/tables/groupA.txt", col.names=T, row.names=T, sep="\t")
write.table(eqFiltB.2, file="saves/tables/groupB.txt", col.names=T, row.names=T, sep="\t")
################

summary(eqFiltA.2)

colnames(eqFiltA.2)
colnames(eqFiltB.2)

eqFilt <- matrix(ncol=dim(eqFiltA.2)[2],nrow=dim(eqFiltA.2)[1]+dim(eqFiltB.2)[1])
colnames(eqFilt) <- gsub('A-([12])$','-\\1',colnames(eqFiltA.2))
rownames(eqFilt) <- c(rownames(eqFiltA.2),rownames(eqFiltB.2))

for (r in levels(as.factor(filesA$Rep))) {
	print(r)
	
	for (i in c(1,2)) {
		eqFilt[,paste(r,i,sep="-")] <- c(as.numeric(eqFiltA.2[,paste(r,'A-',i,sep="")]),as.numeric(eqFiltB.2[,paste(r,'B-',i,sep="")]))
	}
}

dim(eqFilt)

dim(na.exclude(eqFilt))

summary(eqFilt)

hv1<-heatmap.2(
		as.matrix((eqFilt)),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		hclustfun=hclust.diana,
		key=TRUE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=T,
		Colv=T,
		cexRow=0.45,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "All Samples (A and B)"
)


hv1<-heatmap.2(
		as.matrix(eqFilt),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		hclustfun=hclust.ave,
		key=TRUE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=T,
		Colv=T,
		cexRow=0.45,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "All Samples (A and B)"
)


#################
plotCtDensity(qFiltA.2)
plotCtDensity(qFiltB.2)

#plotCtCor(qFiltA.2, main = "Ct correlation - A")
#plotCtCor(qFiltB.2, main = "Ct correlation - B")

#summary(exprs(qFilt.2))

plotCtScatter(dA.norm, cards = c(1, 2), col = "type", diag = TRUE)
plotCtScatter(dB.norm, cards = c(1, 2), col = "type", diag = TRUE)

plotCtPairs(dA.norm, col = "type", diag = TRUE)
plotCtPairs(dB.norm, col = "type", diag = TRUE)

qDE.ttest.A <- ttestCtData(qFiltA.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr" )
#dim(subset(qDE.ttest.A, adj.p.value <= 0.08))
#head(subset(qDE.ttest.A, adj.p.value <= 0.08))

qDE.mwtest.A <- mannwhitneyCtData(qFiltA.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr")
#head(qDE.mwtest.A)
#dim(subset(qDE.mwtest.A, adj.p.value <= 0.15))
#head(subset(qDE.mwtest.A, adj.p.value <= 0.15))

designA <- model.matrix(~0 + filesA$Treatment)
colnames(designA) <- as.character(unique(filesA$Treatment))
contrastsA <- makeContrasts(Down - Control, levels = designA)
colnames(contrastsA) <- c("D-C")
#print(contrastsA)

qFiltA.3 <- qFiltA.2[order(featureNames(qFiltA.2)),]
qDE.limma.A <- limmaCtData(qFiltA.3, design = designA, contrasts = contrastsA, ndups = 1, spacing = 1)
#dim(subset(qDE.limma.A[['D-C']], adj.p.value <= 0.09))
#head(subset(qDE.limma.A[['D-C']], adj.p.value <= 0.09))
head(qDE.limma.A[['D-C']])

qDE.ttest.B <- ttestCtData(qFiltB.2, groups = filesB$Treatment, calibrator = "Control", p.adjust="fdr")
qDE.mwtest.B <- mannwhitneyCtData(qFiltB.2, groups = filesB$Treatment, calibrator = "Control",p.adjust="fdr")


designB <- model.matrix(~0 + filesB$Treatment)
colnames(designB) <- as.character(unique(filesB$Treatment))
contrastsB <- makeContrasts(Down - Control, levels = designB)
colnames(contrastsB) <- c("D-C")

qFiltB.3 <- qFiltB.2[order(featureNames(qFiltB.2)),]
qDE.limma.B <- limmaCtData(qFiltB.3, design = designB, contrasts = contrastsB, ndups = 1, spacing = 1)



write.table(qDE.ttest.A, 
		file="/tmp/NEW_ttestA.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		

write.table(qDE.ttest.B, 
		file="/tmp/NEW_ttestB.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		


write.table(qDE.mwtest.A, 
		file="/tmp/NEW_mwtestA.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		

write.table(qDE.mwtest.B, 
		file="/tmp/NEW_mwtestB.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		


write.table(qDE.limma.A[['D-C']], 
		file="/tmp/NEW_limmaA.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		

write.table(qDE.limma.B[['D-C']], 
		file="/tmp/NEW_limmaB.txt",
		sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)		




colnames(qDE.limma.A[['D-C']])

selected.A <- as.character((subset(qDE.ttest.A, p.value <= 0.08 & adj.p.value <= 0.95 & categoryTarget == 'OK' & categoryCalibrator == 'OK'))[,'genes'])
selected.B <- as.character((subset(qDE.ttest.B, p.value <= 0.08 & adj.p.value <= 0.95 & categoryTarget == 'OK' & categoryCalibrator == 'OK'))[,'genes'])

length(selected.A)
length(selected.B)

dim(eqFilt[union(selected.A,selected.B),])

hv1<-heatmap.2(
		as.matrix(
				eqFilt[union(selected.A,selected.B),]
				),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		breaks=10,
		hclustfun=hclust.ave,
		key=TRUE,
		symkey=FALSE,
		symbreak=FALSE,
		density.info="none",
		trace="none",
		Rowv=T,
		Colv=T,
		cexRow=0.75,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "All Groups (A and B)"
)

min(eqFilt[union(selected.A,selected.B),])


#########
# Análises Joice
#########

joice.analysis <- read.delim(file="Resultados_miRNA.csv", header=TRUE, sep="\t")
dim(joice.analysis)
colnames(joice.analysis)

joice.analysis[,'miRNA']


hv1<-heatmap.2(
		as.matrix(
				eqFilt[which(rownames(eqFilt) %in% joice.analysis[,'miRNA']),]
		),
		na.rm=TRUE,
		scale="none",		
		col=greenred,
		breaks=100,
		hclustfun=hclust.ave,
		key=TRUE,
		symkey=FALSE,
		symbreak=FALSE,
		density.info="none",
		trace="none",
		Rowv=T,
		Colv=T,
		cexRow=0.75,
		cexCol=1,
		keysize=1,
		margins = c(10,10),
		dendrogram=c("both"),
		main = "All Groups (A and B)"
)


joice.miRNAs <- gsub('-[0-9]+$','',joice.analysis[,'miRNA'],perl=TRUE)

length(joice.miRNAs)


library(RmiR.Hs.miRNA)


dbGetQuery(RmiR.Hs.miRNA_dbconn(), "SELECT * FROM tarbase WHERE mature_miRNA='hsa-miR-21'")


tables <- dbListTables(RmiR.Hs.miRNA_dbconn())

dbReadTable(RmiR.Hs.miRNA_dbconn(), 'targetscan')[,1:2]

joice.df <- data.frame(miRNA=unique(joice.miRNAs))

head(joice.df)
rownames(joice.df) <- unique(as.character(joice.df$miRNA))

joice.list <- list()

for( t in c(1:length(tables)) ) {
	print(paste(t,tables[t]))
	targetsdb <- dbReadTable(RmiR.Hs.miRNA_dbconn(), tables[t])[,1:2]
	
	joice.list[[tables[t]]] <- data.frame()
	
	for (m in joice.miRNAs) {
		joice.list[[tables[t]]][m,'targets'] <- vector() 
		
		joice.df[m,tables[t]] <- targetsdb[targetsdb$mature_miRNA %in% m,'gene_id']
	}
	
	#joice.miRNAs <- cbind(joice.miRNAs, targetsdb[targetsdb$mature_miRNA %in% joice.miRNAs])
}


