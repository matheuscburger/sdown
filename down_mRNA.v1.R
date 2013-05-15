# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: down
# Date......: 27/08/2012
# Author(s).: matheus
###############################################################################
analysis.name = "down_mRNA.v1"
## Carregar bibliotecas
#library('dataframes2xls')
library('HTqPCR')
library('genefilter')
library("AgiMicroRna")
library("gplots")
library('ggplot2')
library('reshape')
#library("clusterProfiler")
#library("ReactomePA")
library('RmiR.Hs.miRNA')
library("biomaRt")
library('org.Hs.eg.db')
annpkg <- 'hgug4112a.db'
library(package=annpkg,character.only = TRUE)

source("/work/projects/manalysis/myArray/myArrayFunctions.R")	

## Definir Diretorio base
dir_base = "/work/projects/bit_data"
projeto = "down"
analises = "htqpcr/mRNA"

mywd <- paste(dir_base,projeto,analises,sep="/") 

setwd(mywd)
getwd()



## Ler arquivos
files <- read.delim(file.path(mywd, "targets.txt"))
raw <- readCtData(files = files$File, path = mywd, n.features=96, column.info = list(feature=2, type=3, Ct=4), header=TRUE)
sampleNames(raw) <- gsub('raw/array/', '', sampleNames(raw))
files$Names <- gsub('raw/array/|.txt', '', files$File)
featureClass(raw)	<- files$Treatment
show(raw)


## Plots densidade
png(paste("quality/samples","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(raw)
abline(v=37)
dev.off()

## Definir categorias - Unreliable, Undetermined, Ok
# quantile entre replicas

sink(paste("results/setCategory_groups_Rep_Ctmax_37_quantile_0.75.txt",sep=""))
raw.cat <- setCategory(	raw, 
		Ct.max = 37, 
		groups = files$Rep, 
		quantile = 0.75, 
		replicates=TRUE, 
		verbose=TRUE,
		flag=FALSE,
		plot=T)

sink()

# quantile entre condicoes
sink(paste("results/setCategory_groups_Rep_Treatment_Ctmax_37_quantile_0.90.txt",sep=""))
raw.cat.2 <- setCategory(	raw.cat, 
		Ct.max = 37, 
		groups = files$Treatment,
		quantile = 0.90, 
		replicates=TRUE, 
		verbose=TRUE,
		flag=FALSE,
		plot=TRUE)
sink()


## Plot Categories
system(paste("mkdir -p results/plots/ 2>&1",sep=""), intern=TRUE)
png(filename=paste("results/plots/plotCtCategory_groups_Rep.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(raw.cat)
dev.off()

png(filename=paste("results/plots/plotCtCategory_groups_Rep_Treatment.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(raw.cat.2)
dev.off()

## Filtro - Undetermined
raw.filt <- filterCategory(raw.cat.2, na.categories = c("Undetermined"))	## Inclusao de NAs

## Plots densidade
png(paste("quality/group","filterUndetermined","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(raw.filt)
dev.off()

## Cluster Categories
png(filename=paste("results/plots/Clustering_CtCategory.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)
plotCtCategory(raw.filt, by.feature = TRUE, cexRow = 0.1)
dev.off()

## Normalização
sink(paste("results/normalizeCtData_deltaCt.txt",sep=""))
d.norm <- normalizeCtData(	raw.filt, 
		norm = "deltaCt", 
		deltaCt.genes = c('18S-Hs99999901_s1', 'GAPDH-Hs99999905_m1', 'GUSB-Hs99999908_m1', 'HPRT1-Hs99999909_m1'),
)
sink()


## Plots densidade
png(paste("quality/group","filterUndetermined","norm","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(d.norm)
dev.off()

## Filtro
#sink(paste("results/filterCtDataA_ncat_",as.integer(0.75*length(filesA$Name)),".txt",sep=""))
sink(paste("results/filterCtData",".txt",sep=""))
qFilt <- filterCtData(	d.norm, 
		remove.type = "Endogenous Control", 
		remove.name = c('18S-Hs99999901_s1', 'GAPDH-Hs99999905_m1', 'GUSB-Hs99999908_m1', 'HPRT1-Hs99999909_m1'),
		#remove.category = c("Undetermined","Unreliable"),
		#remove.category = c("Unreliable"),
		#remove.category = c("Undetermined"),
		n.category = 0
#n.category = as.integer(0.75*length(filesA$Name))				
)
sink()

## Plots densidade
png(paste("quality/group","filterUndetermined","norm","filterEndoControl","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(qFilt)
dev.off()

## Obtem valores de expressao
eqFilt <- exprs(qFilt)


## Designa 40 para NAs "Saturação"
#eqFiltA[eqFiltA %in% NA] <- 40
#eqFiltB[eqFiltB %in% NA] <- 40

## Plots densidade
png(paste("quality/groupA","filterUndetermined","norm","saturate","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(dA.norm)
dev.off()

png(paste("quality/groupB","filterUndetermined","norm","saturate","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(dB.norm)
dev.off()

## Avaliação do Interquartile Range
png(filename=paste("results/plots/IQR.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)
iqr.values <- apply(eqFilt, 1, IQR, na.rm=TRUE)
hist(iqr.values, n = 20, main = "", xlab = "IQR across samples")
abline(v = median(iqr.values, na.rm=TRUE), col = 2)
dev.off()


# Nenhum corte sobre o IQR
eqFilt.2 <- eqFilt 
 

## Heatmap
library('gplots')
library('cluster')
library(marray)

rbg = maPalette(low = "green", high = "red", mid = "black", k = 5000)

hclust.ave <- function (d,...) {
	return(hclust(d=d,method="average",...))
}


hclust.diana <- function (x,...) {
	return(as.hclust(diana(x=x,diss=TRUE,...)))
}

#dist.euclidean <- function(x, ...){
#	dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#}

na.reduce <- function(x, max.NA){
	ifelse(sum(is.na(x))>= max.NA,FALSE, TRUE)
}

genes <- apply(eqFilt.2, 1, na.reduce, max.NA=ncol(eqFilt.2)/2)

length(which(genes))

x11();
hv1<-heatmap.2(
		as.matrix(eqFilt.2[which(genes),]),
		na.rm=TRUE,
		scale="col",		
		col=rbg,
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
		main = "Hierarchical Clustering",
		na.color="gray"
)


## Save
source("/work/projects/manalysis/myArray/myArrayFunctions.R")	
#saveArray('saves/')
loadArray('saves/')

write.table(eqFilt.2, file="saves/tables/data.txt", col.names=T, row.names=T, sep="\t")


## Teste estatistico

## Teste MY T

qFilt.2 <- qFilt
exprs(qFilt.2) <- eqFilt.2
qDE.ttest <- my.ttestCtData(qFilt.2, groups = files$Treatment, calibrator = "Control", p.adjust="fdr", na.rm=TRUE, stringent=FALSE)

#dim(subset(qDE.ttest, adj.p.value <= 0.05))
#head(subset(qDE.ttest, adj.p.value <= 0.05))
write.table(qDE.ttest, paste("results","tables","my.ttest.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)

### Teste MY Mann-Whitney
#
qDE.mwtest <- my.mannwhitneyCtData(qFilt.2, groups = files$Treatment, calibrator = "Control", p.adjust="fdr", stringent=FALSE)
write.table(qDE.mwtest, paste("results","tables","my.mwtest.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
#tail(qDE.mwtest)
#summary(qDE.mwtest$FC)

## Limma

design <- model.matrix(~0 + files$Treatment)
colnames(design) <- as.character(unique(files$Treatment))
contrasts <- makeContrasts(Down - Control, levels = design)
colnames(contrasts) <- c("D-C")

qDE.limma <- limmaCtData(qFilt.2, design = design, contrasts = contrasts, ndups = 1, spacing = NULL, stringent=FALSE)
head(qDE.limma[['D-C']])
write.table(qDE.limma[['D-C']], paste("results","tables","limma.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)


## Limma modified
qDE.my.limma <- my.limmaCtData(qFilt.2, design = design, contrasts = contrasts, ndups = 1, spacing = NULL, NA.threshold = 0.6)
write.table(qDE.my.limma[['D-C']], paste("results","tables","my.limma.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
my.DE <- subset(qDE.my.limma[['D-C']], categoryCalibrator=="OK" & categoryTarget=="OK" )
write.table(my.DE, paste("results","tables","my.limma.A.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.filtered.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)


## Save
source("/work/projects/manalysis/myArray/myArrayFunctions.R")	
#saveArray('saves/')
loadArray('saves/')

#dim(qDE.my.limma[['D-C']]) - dim(my.DE)



de.probes <- levels(my.DE$genes)[my.DE$genes[which(my.DE$adj.p.value <= 0.05)]]
de.genes <- unlist(lapply(strsplit(de.probes, "-"), function(x){x[[1]]}))


t.CC <- apply(aux[,grep("CC", colnames(aux))], 1, t.test)
ci.CC<- data.frame(ci.inf=unlist(lapply(t.CC, function(x){x$conf.int[1]})), ci.sup=unlist(lapply(t.CC, function(x){x$conf.int[2]})))
ci.CC$name <- row.names(ci.CC)

t.CD <- apply(aux[,grep("CD", colnames(aux))], 1, t.test)
ci.CD<- data.frame(ci.inf=unlist(lapply(t.CD, function(x){x$conf.int[1]})), ci.sup=unlist(lapply(t.CD, function(x){x$conf.int[2]})))
ci.CD$name <- row.names(ci.CD)

mean.CC <- as.data.frame(apply(aux[,grep("CC", colnames(aux))], 1, mean, na.rm=T))
mean.CC$name <- rownames(mean.CC)
mean.CD <- as.data.frame(apply(aux[,grep("CD", colnames(aux))], 1, mean, na.r=T))
mean.CD$name <- rownames(mean.CD)

cc <- merge(mean.CC, ci.CC, by='name')
colnames(cc) <- c("name", "mean", "ci.inf", "ci.sup")
cc$cond <- "control"
cd <- merge(mean.CD, ci.CD, by='name')
colnames(cd) <- c("name", "mean", "ci.inf", "ci.sup")
cd$cond <- "down"

plot.data <- rbind(cc,cd)




aux <- eqFilt.2[de.probes,]
rownames(aux) <- de.genes 
aux.m <- melt(aux)
colnames(aux.m) <- c("gene", "sample", "value")
aux.m$cond <- apply(aux.m, 1, function(x){
			res <- "Down"			
			if(length(grep("CD",x["sample"])) == 0 )
				res <- "Control"
			res
		})

#subset(aux.m, cond=="Down")
write.table(aux.m, "/tmp/de.genes.limma.txt", col.names=T, quote=F, sep="\t")
aux.m.sum <- summarySE(aux.m, measurevar="value", groupvars=c("gene","cond"), na.rm=T)
aux.m.sum$order <- order(aux.m.sum$se)

# Standard error of the mean
ggplot(plot.data, aes(x=name, y=mean, colour=cond, order=mean)) + 
		geom_errorbar(aes(ymin=ci.inf, ymax=ci.sup), width=.2) +
		geom_point(size=2.5) +
		opts(axis.text.x=theme_text(angle=-45))
ggsave(file = paste("results/plots/de.genes.ic.png",sep=''), width=10, height=10)


		
### INTERSECT (miR, mR)

miR.intersect <- read.table("/data/projects/bit_data/down/htqpcr/intersect/genes.txt", header=F, colClasses="character")[,1]

genes <-  levels(qDE.ttest$genes)[qDE.ttest$genes]
genes <- unlist(lapply(strsplit(genes, "-"), function(x){x[1]}))
qDE.ttest[which(genes %in% miR.intersect), ]	# resultado do testes para genes em intersecção

genes <-  levels(qDE.my.limma[['D-C']]$genes)[qDE.my.limma[['D-C']]$genes]
genes <- unlist(lapply(strsplit(genes, "-"), function(x){x[1]}))
qDE.my.limma[['D-C']][which(genes %in% miR.intersect), ]	# resultado do testes para genes em intersecção

int.probes <- levels(qDE.my.limma[['D-C']][which(genes %in% miR.intersect), "genes"])[qDE.my.limma[['D-C']][which(genes %in% miR.intersect), "genes"]]

aux <- eqFilt.2[int.probes,]
rownames(aux) <- unlist(lapply(strsplit(int.probes, "-"), function(x){x[1]}))
aux.m <- melt(aux)
colnames(aux.m) <- c("gene", "sample", "value")
aux.m$cond <- apply(aux.m, 1, function(x){
			res <- "Down"			
			if(length(grep("CD",x["sample"])) == 0 )
				res <- "Control"
			res
		})

#subset(aux.m, cond=="Down")

aux.m.sum <- summarySE(aux.m, measurevar="value", groupvars=c("gene","cond"), na.rm=T)
aux.m.sum$order <- order(aux.m.sum$se)

# Standard error of the mean
g <- ggplot(aux.m.sum, aes(x=gene, y=value, colour=cond, order=value)) + 
		geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2) +
		geom_point(size=2.5 )
ggsave("/data/projects/bit_data/down/htqpcr/intersect/finaldot.png", width=10, height=10) 

