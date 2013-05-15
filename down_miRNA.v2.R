# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: down
# Date......: 22/05/2012
# Author(s).: matheus
###############################################################################

## Carregar bibliotecas
library('dataframes2xls')
library('HTqPCR')
library('genefilter')
library("AgiMicroRna")
library("miRNApath")
library("gplots")
library('reshape')
#library("clusterProfiler")
#library("ReactomePA")
library('RmiR.Hs.miRNA')

library('org.Hs.eg.db')
annpkg <- 'hgug4112a.db'

library(package=annpkg,character.only = TRUE)

source("/work/projects/manalysis/myArray/myArrayFunctions.R")	

## Definir Diretorio base
dir_base = "/work/projects/bit_data"
projeto = "down"
analises = "htqpcr/miRNA"

mywd <- paste(dir_base,projeto,analises,sep="/") 

setwd(mywd)
getwd()

## Ler arquivos
# Grupo A
filesA <- read.delim(file.path(mywd, "targetsA384.txt"))
rawA <- readCtData(files = filesA$File, path = mywd, n.features=384, feature=2, type=3, Ct=4, header=TRUE)
sampleNames(rawA) <- gsub('raw/A/', '', sampleNames(rawA))
filesA$Names <- gsub('raw/A/|.txt', '', filesA$File)
featureClass(rawA)	<- filesA$Treatment
show(rawA)

# Grupo B
filesB <- read.delim(file.path(mywd, "targetsB384.txt"))
rawB <- readCtData(files = filesB$File, path = mywd, n.features=384, feature=2, type=3, Ct=4, header=TRUE)
sampleNames(rawB) <- gsub('raw/B/', '', sampleNames(rawB))
filesB$Names <- gsub('raw/B/|.txt', '', filesB$File)
featureClass(rawB)	<- filesB$Treatment
show(rawB)


## Plots densidade

png(paste("quality/groupA","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(rawA)
dev.off()

png(paste("quality/groupB","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(rawB)
dev.off()

## Definir categorias - Unreliable, Undetermined, Ok
## Grupo A
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
		plot=F)

sink()

## Grupo B
# quantile entre replicas
sink(paste("results/setCategoryB_groups_Rep_Ctmax_32_quantile_0.75.txt",sep=""))				
rawB.cat <- setCategory(	rawB, 
		Ct.max = 32, 
		groups = filesB$Rep, 
		quantile = 0.75,
		replicates=TRUE, 
		verbose=TRUE,
		flag=FALSE,						
		plot=F)
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
		plot=F)
sink()

par(ask=F, mfrow=c(1,1))
## Plot Categories

system(paste("mkdir -p results/plots/ 2>&1",sep=""), intern=TRUE)

png(filename=paste("results/plots/plotCtCategoryA_groups_Rep.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(rawA.cat)
dev.off()

png(filename=paste("results/plots/plotCtCategoryB_groups_Rep.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(rawB.cat)
dev.off()


png(filename=paste("results/plots/plotCtCategoryA_groups_Rep_Treatment.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(rawA.cat.2)
dev.off()

png(filename=paste("results/plots/plotCtCategoryB_groups_Rep_Treatment.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)				
plotCtCategory(rawB.cat.2)
dev.off()

## Filtro - Undetermined

rawA.filt <- filterCategory(rawA.cat.2, na.categories = c("Undetermined"))	## Inclusao de NAs
names(which(apply(exprs(rawA.filt),1,function(x){ is.element(NA,x) })==TRUE))
names(which(apply(featureCategory(rawA.cat.2),1,function(x){ is.element("Undetermined",x) })==TRUE))
write.table(names(which(apply(featureCategory(rawA.cat.2),1,function(x){ is.element("Undetermined",x) })==TRUE)),
		"/tmp/verify/mymicrosA.txt", row.names=F, col.names=F, quote=F)

rawB.filt <- filterCategory(rawB.cat.2, na.categories = c("Undetermined"))	## Inclusao de NAs 
write.table(names(which(apply(exprs(rawB.filt),1,function(x){ is.element(NA,x) })==TRUE)), "/tmp/microsUndeterminedB", row.names=F, quote=F)
## Plots densidade

png(paste("quality/groupA","filterUndetermined","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(rawA.filt)
dev.off()

png(paste("quality/groupB","filterUndetermined","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(rawB.filt)
dev.off()

## Cluster Categories
png(filename=paste("results/plots/Clustering_CtCategoryA.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)
plotCtCategory(rawA.filt, by.feature = TRUE, cexRow = 0.1)
dev.off()

png(filename=paste("results/plots/Clustering_CtCategoryB.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)
plotCtCategory(rawB.filt, by.feature = TRUE, cexRow = 0.1)
dev.off()

## Normalização

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

## Plots densidade

png(paste("quality/groupA","filterUndetermined","norm","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(dA.norm)
dev.off()

png(paste("quality/groupB","filterUndetermined","norm","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(dB.norm)
dev.off()

## Filtro

#sink(paste("results/filterCtDataA_ncat_",as.integer(0.75*length(filesA$Name)),".txt",sep=""))
sink(paste("results/filterCtDataA",".txt",sep=""))
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

#sink(paste("results/filterCtDataB_ncat_",as.integer(0.75*length(filesB$Name)),".txt",sep=""))						
sink(paste("results/filterCtDataB",".txt",sep=""))
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

## Plots densidade
png(paste("quality/groupA","filterUndetermined","norm","filterEndoControl","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(qFiltA)
dev.off()

png(paste("quality/groupB","filterUndetermined","norm","filterEndoControl","png", sep="."), res =300, height=2000, width=2000)
plotCtDensity(qFiltB)
dev.off()

## Obtem valores de expressao

eqFiltA <- exprs(qFiltA)
eqFiltB <- exprs(qFiltB)

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



######
rownames(eqFiltB) %in% rownames(eqFiltA)

write(rownames(eqFiltB),"/tmp/groupB")
write(rownames(eqFiltA),"/tmp/groupA")

## Avaliação do Interquartile Range
png(filename=paste("results/plots/IQRA.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)

iqrA.values <- apply(eqFiltA, 1, IQR, na.rm=TRUE)
hist(iqrA.values, n = 20, main = "", xlab = "IQR across samples - A")
abline(v = median(iqrA.values, na.rm=TRUE), col = 2)
dev.off()

png(filename=paste("results/plots/IQRB.png",sep=""), 
		bg="white", res=300, width=2000, height=2000)
iqrB.values <- apply(eqFiltB, 1, IQR, na.rm=TRUE)
hist(iqrB.values, n = 20, main = "", xlab = "IQR across samples - B")
abline(v = median(iqrB.values, na.rm=TRUE), col = 2)
dev.off()

# Nenhum corte sobre o IQR
eqFiltA.2 <- eqFiltA 
eqFiltB.2 <- eqFiltB 

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

genesA <- apply(eqFiltA.2, 1, na.reduce, max.NA=ncol(eqFiltA.2)/2)

length(which(genesA))

x11();
hv1<-heatmap.2(
		as.matrix(eqFiltA.2[which(genesA),]),
		na.rm=TRUE,
		scale="none",		
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
		main = "Group A",
		na.color="gray"
)


genesB <- apply(eqFiltB.2, 1, na.reduce, max.NA=ncol(eqFiltB.2)/2)

length(which(genesB))

hv1<-heatmap.2(
		as.matrix(eqFiltB.2[which(genesB),]),
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

## Save
source("/work/projects/manalysis/myArray/myArrayFunctions.R")	
saveArray('saves/')

write.table(eqFiltA.2, file="saves/tables/groupA.txt", col.names=T, row.names=T, sep="\t")
write.table(eqFiltB.2, file="saves/tables/groupB.txt", col.names=T, row.names=T, sep="\t")


# Junção Grupos A e B

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

## Teste estatistico

## Teste T
#
#qFiltA.2 <- qFiltA
#exprs(qFiltA.2) <- eqFiltA.2
#qDE.ttest.A <- ttestCtData(qFiltA.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr", na.rm=TRUE, stringent=TRUE)
##dim(subset(qDE.ttest.A, adj.p.value <= 0.08))
##head(subset(qDE.ttest.A, adj.p.value <= 0.08))
#write.table(qDE.ttest.A, paste("results","tables","ttest.A.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
#
#qFiltB.2 <- qFiltB
#exprs(qFiltB.2) <- eqFiltB.2
#qDE.ttest.B <- ttestCtData(qFiltB.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr" , stringent=FALSE)
#write.table(qDE.ttest.B, paste("results","tables","ttest.B.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
#
#
### Teste Mann-Whitney
#
#qDE.mwtest.A <- mannwhitneyCtData(qFiltA.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr", stringent=TRUE)
#write.table(qDE.mwtest.A, paste("results","tables","mwtest.A.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
#head(qDE.mwtest.A)
#
#qDE.mwtest.B <- mannwhitneyCtData(qFiltB.2, groups = filesA$Treatment, calibrator = "Control", p.adjust="fdr", stringent=TRUE)
#write.table(qDE.mwtest.B, paste("results","tables","mwtest.B.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
#head(qDE.mwtest.A)

## Limma

designA <- model.matrix(~0 + filesA$Treatment)
colnames(designA) <- as.character(unique(filesA$Treatment))
contrastsA <- makeContrasts(Down - Control, levels = designA)
colnames(contrastsA) <- c("D-C")

qDE.limma.A <- limmaCtData(qFiltA.2, design = designA, contrasts = contrastsA, ndups = 1, spacing = NULL, stringent=FALSE)
head(qDE.limma.A[['D-C']])
write.table(qDE.limma.A[['D-C']], paste("results","tables","limma.A.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)

designB <- model.matrix(~0 + filesB$Treatment)
colnames(designB) <- as.character(unique(filesB$Treatment))
contrastsB <- makeContrasts(Down - Control, levels = designB)
colnames(contrastsB) <- c("D-C")

qDE.limma.B <- limmaCtData(qFiltB.2, design = designB, contrasts = contrastsB, ndups = 1, spacing = NULL, stringent=FALSE)
write.table(qDE.limma.B[['D-C']], paste("results","tables","limma.B.cal_control.padj_fdr.stringent_F.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)

# x11();func.list$volcano.plus( qDE.limma.B[["D-C"]], "FC", "p.value", paste("Volcano Plot (limma) \n "), 0.01, 2, -2,
#		'adj.p.value' ,ncolors=5)

## Limma modified
	
qDE.my.limma.A <- my.limmaCtData(qFiltA.2, design = designA, contrasts = contrastsA, ndups = 1, spacing = NULL, NA.threshold = 0.6)
write.table(qDE.my.limma.A[['D-C']], paste("results","tables","my.limma.A.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
my.DE.A <- subset(qDE.my.limma.A[['D-C']], categoryCalibrator=="OK" & categoryTarget=="OK" )
write.table(my.DE.A, paste("results","tables","my.limma.A.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.filtered.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)

qDE.my.limma.B <- my.limmaCtData(qFiltB.2, design = designB, contrasts = contrastsB, ndups = 1, spacing = NULL, NA.threshold = 0.6)
write.table(qDE.my.limma.B[['D-C']], paste("results","tables","my.limma.B.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)
my.DE.B <- subset(qDE.my.limma.B[['D-C']], categoryCalibrator=="OK" & categoryTarget=="OK" )
write.table(my.DE.B, paste("results","tables","my.limma.B.cal_control.padj_fdr.stringent_F.NAThreshold_0.6.filtered.txt",sep="/"), col.names=T, row.names=F, sep="\t", quote=F)


## Avaliação

hist(my.DE.A$p.value, breaks=c(seq(0,1,0.01)))
abline(v = 0.05, col = 2)
hist(my.DE.A$adj.p.value, breaks=c(seq(0,1,0.01)))
abline(v = 0.05, col = 2)

hist(my.DE.B$p.value, breaks=c(seq(0,1,0.01)))
abline(v = 0.05, col = 2)
hist(my.DE.B$adj.p.value, breaks=c(seq(0,1,0.01)))
abline(v = 0.05, col = 2)

pv.cut <- 0.05
fc.cut <- 1.5

miR.DE.A <- subset(my.DE.A, p.value<pv.cut & ((FC >= fc.cut) | (FC <= 1/fc.cut)) )
head(miR.DE.A)
write.table(miR.DE.A$genes, file=paste("results","tables",paste("sig.miR.A.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=F, quote=F)
write.table(miR.DE.A, file=paste("results","tables",paste("sig.miR.values.A.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=T, quote=F, sep="\t")

miR.DE.B <- subset(my.DE.B, p.value<pv.cut & ((FC >= fc.cut) | (FC <= 1/fc.cut)) )
write.table(miR.DE.B$genes, file=paste("results","tables",paste("sig.miR.B.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=F, quote=F)
write.table(miR.DE.B, file=paste("results","tables",paste("sig.miR.values.B.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=T, quote=F, sep="\t")

dim(miR.DE.A)
dim(miR.DE.B)

#aux <- miR.DE.B[,c("genes","FC","p.value")]
#aux[which(aux$FC<1),'FC'] <- -1/aux[which(aux$FC<1),'FC']
#aux$FC <- round(aux$FC,2)
#aux$p.value <- round(aux$p.value,4)
#write.table(aux, "/tmp/myauxB.tmp", row.names=F, col.names=T, quote=F, sep="|",dec=",")
#
#aux <- miR.DE.A[,c("genes","FC","p.value")]
#aux[which(aux$FC<1),'FC'] <- -1/aux[which(aux$FC<1),'FC']
#aux$FC <- round(aux$FC,2)
#aux$p.value <- round(aux$p.value,4)
#write.table(aux, "/tmp/myauxA.tmp", row.names=F, col.names=T, quote=F, sep="|",dec=",")


## Heatmaps

x11()
aux.miRsA <- levels(miR.DE.A$genes)[miR.DE.A$genes]
hv1<-heatmap.2(
		as.matrix(exprs(qFiltA.2)[aux.miRsA,]),
		na.rm=TRUE,
		scale="none",		
		col=rbg,
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
		main = "Group A",
		na.color="gray"
)
rm(aux.miRsA)

x11()
aux.miRsB <- levels(miR.DE.B$genes)[miR.DE.B$genes]
hv1<-heatmap.2(
		as.matrix(exprs(qFiltB.2)[aux.miRsB,]),
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
		main = "Group B",
		na.color="grey"
)
rm(aux.miRsB)


source("/work/projects/manalysis/myArray/myArrayFunctions.R")	
saveArray('saves/')
#loadArray('saves/')

#########
loadArray('saves/')


## miRNApath
#dataA <- (my.DE.A)
#dataB <- (my.DE.B)
#
#head(dataA)
#
#colnames(dataA)
#
### Converter ids (coluna genes) / remover ultima coluna
#
#aux.ids <- levels(dataA$genes)[dataA$genes]
#head(aux.ids)
#aux.ids <- unlist(strsplit(aux.ids, "\\-\\d+$",perl=TRUE))
#sustenido <- grep("#$", aux.ids)
#aux.ids[sustenido] <- gsub(".$","",aux.ids[sustenido])
#dataA$Id <- aux.ids
#
#rm(aux.ids)
#
#aux.ids <- levels(dataB$genes)[dataB$genes]
#head(aux.ids)
#aux.ids <- unlist(strsplit(aux.ids, "\\-\\d+$",perl=TRUE))
#sustenido <- grep("#$", aux.ids)
#aux.ids[sustenido] <- gsub(".$","",aux.ids[sustenido])
#dataB$Id <- aux.ids
#
### Escrever arquivo
#
#dataA$Assay <- dataA$Id
#
#dataA$Group <- "DOWN.vs.CONTROL"
#
#aux.fc <- dataA$FC;  aux.fc[which(aux.fc < 1)] <- -1/aux.fc[which(aux.fc < 1)];  dataA$FC <- aux.fc
#
#file.dataA <- paste("/tmp",paste("mirnaA",format(Sys.time(),"%H%M%S"),"txt",sep="."),sep="/")
#write.table(dataA,file=file.dataA, quote=FALSE, row.names=FALSE, col.names=TRUE, na="", sep="\t")
#
#dataB$Assay <- dataB$Id
#
#dataB$Group <- "DOWN.vs.CONTROL"
#
#aux.fc <- dataB$FC;  aux.fc[which(aux.fc < 1)] <- -1/aux.fc[which(aux.fc < 1)];  dataB$FC <- aux.fc
#
#file.dataB <- paste("/tmp",paste("mirnaB",format(Sys.time(),"%H%M%S"),"txt",sep="."),sep="/")
#write.table(dataB,file=file.dataB, quote=FALSE, row.names=FALSE, col.names=TRUE, na="", sep="\t")
#
### Criar objeto
#
#mirnaobjA <- loadmirnapath( mirnafile=file.dataA, pvaluecol="p.value", groupcol="Group", mirnacol="Id", assayidcol= "Assay", foldchangecol = "FC")
#
#mirnaobjA
#mirnaobjA@columns["pvaluecol"] <- "p.value";
#mirnaobjA@columns["foldchangecol"] <- "FC";
#
#mirnaobjA <- filtermirnapath( mirnaobjA, pvalue = pv.cut, expression = NA, foldchange = fc.cut );
#
#mirnaobjA
#head(mirnaobjA@mirnaTable)
#tail(mirnaobjA@mirnaTable)
#
#
#mirnaobjB <- loadmirnapath( mirnafile=file.dataB, pvaluecol="p.value", groupcol="Group", mirnacol="Id", assayidcol= "Assay", foldchangecol = "FC")
#
#mirnaobjB
#mirnaobjB@columns["pvaluecol"] <- "p.value";
#mirnaobjB@columns["foldchangecol"] <- "FC";
#
#mirnaobjB <- filtermirnapath( mirnaobjB, pvalue = pv.cut, expression = NA, foldchange = fc.cut );
#
#mirnaobjB
#head(mirnaobjB@mirnaTable)
#tail(mirnaobjB@mirnaTable)
#
#
#### Obter Alvos
### Utilizando script em perl para gerar uma tabela contendo as interações miR-mR
###
#
#
#miRs <- c(dataA$Id, dataB$Id)
#length(miRs)	# 374
#write.table(as.character(miRs), file="results/tables/miRs.txt", row.names=F, col.names=F, quote=F)
#
#aux.ids <-as.character(c(levels(miR.DE.A$genes)[miR.DE.A$genes],levels(miR.DE.B$genes)[miR.DE.B$genes]))
#head(aux.ids)
#aux.ids <- unlist(strsplit(aux.ids, "\\-\\d+$",perl=TRUE))
#sustenido <- grep("#$", aux.ids)
#aux.ids[sustenido] <- gsub(".$","",aux.ids[sustenido])
#write.table(as.character(aux.ids), file="results/tables/miRs.selected.txt", row.names=F, col.names=F, quote=F)
#
#
##library("RmiR.Hs.miRNA")
##library("RmiR")
##
##
##con <- RmiR.Hs.miRNA_dbconn()
##
##tables <- dbListTables(con)
##
##ls()
##
##
##
##my.micros <- unlist(strsplit(c(levels(miR.DE.A$genes),levels(miR.DE.B$genes)), "\\-\\d+$",perl=TRUE))
##my.genes <- c()
##for(t in tables){
##	assign(t,unique(dbReadTable(con, t)[,c("mature_miRNA","gene_id")]))
##	my.genes <- unique(c(my.genes,get(t)[which(get(t)[,'mature_miRNA'] %in% my.micros), "gene_id"]))
##
##}
##miRmR <- matrix(0, nrow=length(my.micros), ncol=length(my.genes), dimnames = list(my.micros, as.character(my.genes) ))
##for(t in tables){
##	pos.in.table <- which(get(t)[,'mature_miRNA'] %in%  my.micros )
##	in.matrix <- get(t)[pos.in.table, c("mature_miRNA","gene_id")]
##	
###	miRmR[in.matrix[,1], in.matrix[,2]] <- 
###			miRmR[in.matrix[,1], in.matrix[,2]] + 1
##	
##	for(i in rownames(in.matrix) ){
##		miRmR[in.matrix[i,1], in.matrix[i,2]] <- miRmR[in.matrix[i,1], in.matrix[i,2]] + 1  
##		
##	}
##}
##miRmR[relation.mir,]
##relation.mir[1]
##get(t)[3,]
##miRmR[3, as.character(get(t)[3,'gene_id'])]
##
##
##summary(miRmR)
##head(miRmR)
##
##max(miRmR)
##
##
##which(my.micros %in% miranda$mature_miRNA) 
##head(miranda)
##
##
##head(my.micros)
#
#miRmR <- read.table("/work/matheus/getmiRmR/out/miRmR.txt")
#dim(miRmR)
#
#miRmR[1:10,1:10]
#
#max(miRmR)
#min(miRmR)
#
#aux.sum <- apply(miRmR, 2, sum)
#min(aux.sum)
#
#
#entrez <- colnames(miRmR); substr(entrez, 1,1) <- ""
#class(entrez)
#
#
#
### Data structure with Symbol names for each entrez id
#x <- org.Hs.egSYMBOL
## Get the gene names that are mapped to an entrez gene identifier
#mapped_genes <- mappedkeys(x)
## Convert to a list
#SYMBOL <- as.list(x[mapped_genes])
#my.genes <- unlist(SYMBOL[as.character(entrez)])
#
#table(my.genes)
#
#### Subset miRmR
#### MicroRNAs no set A presentes no banco de dados de interação miRmR
##relations.micro <- which(rownames(miRmR) %in% dataA$Id)
##length(relations.micro)
##length(dataA$Id)
##dataA$Id[which(!dataA$Id %in% rownames(miRmR))]         #"has-miR-155" => has != hsa
##
##miRmR.sub <- miRmR[relations.micro,]
##
##
##
##miRmR.sub.atl2 <- miRmR.sub; 
##miRmR.sub.atl2[miRmR.sub.atl2<2] <- 0 # at least 2 databases
##
##
##write.table(as.character(dataA$Id), file="results/tables/miRs.setA.txt", row.names=F, col.names=F, quote=F)
##
#
#
### Format Table
#
#dim(miRmR)
#
### Load the miRNA to gene associations
#mirnaobj <- loadmirnatogene( mirnafile="mirnaGene.txt", mirnaobj=mirnaobj, mirnacol="miRNA Name",
#	genecol="Entrez Gene ID",columns=c(assayidcol="ASSAYID") );
#mirnaobj;
#



#############
### New Analysis

## Obter microRNAs exclusivos

control <- grep("CC",colnames(featureCategory(qFiltA.2)))
down <- grep("CD",colnames(featureCategory(qFiltA.2)))
head(featureCategory(qFiltA.2))

## Verifica se todas as amostras ->> down <<- são categorizadas como Unreliable ou Undetermined
## então provavelmente esse microRNA eh exclusivo do ->> controle <<-
aux.control.exclusive.A <- names(which(apply(featureCategory(qFiltA.2), 1,function(x){
							is.na(table(as.character(x[down]))["OK"])
						})))
length( aux.control.exclusive.A )

aux.down.exclusive.A <- names(which(apply(featureCategory(qFiltA.2), 1,function(x){
					is.na(table(as.character(x[control]))["OK"]) ## Não existe OK ??
				})))
length( aux.down.exclusive.A )

control.exclusive.A <- (setdiff(aux.control.exclusive.A, aux.down.exclusive.A))
## Verificação
featureCategory(qFiltA.2)[control.exclusive.A, control]
apply(featureCategory(qFiltA.2)[control.exclusive.A,control],1,function(x){table(x)["OK"]}) ## Não deve possuir NA
apply(featureCategory(qFiltA.2)[control.exclusive.A,down],1,function(x){table(x)["OK"]}) ## NA

exprs(qFiltA.2)[control.exclusive.A,control]
exprs(qFiltA.2)[control.exclusive.A,down]

mean(exprs(rawA)[control.exclusive.A,control])	## Media dados BRUTOS
sd(as.numeric(exprs(rawA)[control.exclusive.A,control]))
mean(exprs(rawA)[control.exclusive.A,down])	## Media dados BRUTOS
sd(as.numeric(exprs(rawA)[control.exclusive.A,down]))

down.exclusive.A <- (setdiff(aux.down.exclusive.A, aux.control.exclusive.A))
## Verificação
featureCategory(qFiltA.2)[down.exclusive.A, down]
apply(featureCategory(qFiltA.2)[down.exclusive.A,down],1,function(x){table(x)["OK"]}) ## Não deve possuir NA
apply(featureCategory(qFiltA.2)[down.exclusive.A,control],1,function(x){table(x)["OK"]}) ## NA

exprs(qFiltA.2)[control.exclusive.A,control]
exprs(qFiltA.2)[control.exclusive.A,down]

mean(exprs(rawA)[down.exclusive.A,control])	## Media dados BRUTOS
sd(as.numeric(exprs(rawA)[down.exclusive.A,control]))
mean(exprs(rawA)[down.exclusive.A,down])	## Media dados BRUTOS
sd(as.numeric(exprs(rawA)[down.exclusive.A,down]))

#######
colnames(qDE.my.limma.A[["D-C"]])
control.exclusive.A <- subset(qDE.my.limma.A[["D-C"]],categoryTarget=="Undetermined" & categoryCalibrator=="OK")
control.exclusive.A$LFC <- log2(control.exclusive.A$FC)
control.exclusive.A.sig <- subset(control.exclusive.A, p.value <= 0.05 & (LFC>1 | LFC< -1))
featureCategory(qFiltA.2)[levels(control.exclusive.A.sig$genes)[control.exclusive.A.sig$genes], down]
exprs(qFiltA.2)[levels(control.exclusive.A.sig$genes)[control.exclusive.A.sig$genes], down]
exprs(qFiltA.2)[levels(control.exclusive.A.sig$genes)[control.exclusive.A.sig$genes], control]

control.exclusive.B <- subset(qDE.my.limma.B[["D-C"]],categoryTarget=="Undetermined" & categoryCalibrator=="OK")
control.exclusive.B$LFC <- log2(control.exclusive.B$FC)
control.exclusive.B.sig <- subset(control.exclusive.B, p.value <= 0.05 & (LFC>1 | LFC< -1))

down.exclusive.A <- subset(qDE.my.limma.A[["D-C"]],categoryTarget=="OK" & categoryCalibrator=="Undetermined")
down.exclusive.A$LFC <- log2(down.exclusive.A$FC)
down.exclusive.A.sig <- subset(down.exclusive.A, p.value <= 0.05 & (LFC>1 | LFC< -1))

down.exclusive.B <- subset(qDE.my.limma.B[["D-C"]],categoryTarget=="OK" & categoryCalibrator=="Undetermined")
down.exclusive.B$LFC <- log2(down.exclusive.B$FC)
down.exclusive.A.sig <- subset(down.exclusive.A, p.value <= 0.05 & (LFC>1 | LFC< -1))


####################

## pegar miRs com NA em mais do que 20 % das amostras em um grupo
## pegar a expressão no outro grupo e verificar se esta acima da mediana do array
median.A <- median(as.numeric(exprs(qFiltA.2)), na.rm=T)
median.B<- median(as.numeric(exprs(qFiltB.2)), na.rm=T)


x <- down.exclusive.A
qPCRset <- qFiltA.2
expr.cut <- median.A
group.idx=control
group.names=sampleNames(qFiltA)[control]

#'qPCRset
#'group.names = nome das amostras do grupo com exclusividade de expressão
#'down.cut = valor de 0 a 1 que indica a parcela minima com valores categorizados como OK no grupo com exclusividade
#'up.cut = valor de 0 a 1 que indica a parcela minima com valores categorizados como Undetermined no outro grupo
verify.exclusive <- function(qPCRset, group.names=sampleNames(qFiltA)[control], expr.cut=10, up.cut=0.8, down.cut=0.6){
	if(as.character(class(qPCRset))!="qPCRset")
		stop("Please, think a little ...")
	## Verificar microRNAs com pelo menos down.cut% de expressão no grupo exclusivo
	min.ok <- ceiling(0.6*length(group.names))	# quantidade minima de ok no grupo exlusivo
	miR.ok <- names(which(apply(featureCategory(qPCRset)[group.names], 1, function(x){
					table(x)["OK"]>= min.ok
				})))
	other.grp <- setdiff(sampleNames(qPCRset),group.names)
	min.undet <- ceiling(0.8*length(other.grp))		## quantidade minima de valores indeterminados para grupo ser considerado sem expressao
	miR.pos <- names(which(apply(featureCategory(qPCRset)[miR.ok,other.grp],1,
			function(x){
				table(x)["Undetermined"]>= min.undet
			}
			)))
	
	print(group.names)
	print(other.grp)
	print(miR.pos)
	if(length(miR.pos)>0){
		if(length(miR.pos)==1){
			if(mean(exprs(qPCRset)[miR.pos, group.names], na.rm=TRUE) <= expr.cut )
				return(miR.pos)
			else
				return(NA)
		}
		return(names(which(apply(exprs(qPCRset)[miR.pos, group.names],1,function(x){
					mean(x, na.rm=T) <= expr.cut
				}))))
	}else{
		
		return(NA)
	}
	
}

median(as.numeric(exprs(qFiltA.2)), na.rm=T)
median(as.numeric(exprs(qFiltB.2)), na.rm=T)
my.par <- list(  list(qPCRset=qFiltA.2, group.names=sampleNames(qFiltA.2)[control], expr.cut=median(as.numeric(exprs(qFiltA.2)), na.rm=T)),
		list(qPCRset=qFiltA.2, group.names=sampleNames(qFiltA.2)[down], expr.cut=median(as.numeric(exprs(qFiltA.2)), na.rm=T)),
		list(qPCRset=qFiltB.2, group.names=sampleNames(qFiltB.2)[control], expr.cut=median(as.numeric(exprs(qFiltB.2)), na.rm=T)),
		list(qPCRset=qFiltB.2, group.names=sampleNames(qFiltB.2)[control], expr.cut=medialibrary('HTqPCR')


i <- my.par[[1]]
for(i in my.par){	
	print(i[[1]])
	print(verify.exclusive(i[['qPCRset']], i[['group.names']], i[['expr.cut']], 0.5, 0.6))
}
		

## Obter relação miRmR
# para todos os miRNAs presentes na análise (setA e setB)
## qFiltA.2 ## dados utilizados para o teste estatistico
## qFiltB.2 ## dados utilizados para o teste estatistico

## todos os microRNAs
all.miRs <- c(featureNames(qFiltA.2), featureNames(qFiltB.2))
all.miRs <- unlist(strsplit(all.miRs, "\\-\\d+$",perl=TRUE))
sustenido <- grep("#$", all.miRs)
all.miRs[sustenido] <- gsub(".$","",all.miRs[sustenido])
length(all.miRs)

all.miRs <- unique(all.miRs)

length(all.miRs)

write.table(all.miRs,"/data/matheus/getmiRmR/in/all.miRs.txt", col.names=F, quote=F, row.names=F)

## Rodar script em perl (process.pl) para obter tabela miRmR

miRmR <- read.table("/work/matheus/getmiRmR/out/miRmR.txt")

dim(miRmR)

length(all.miRs) - nrow(miRmR) # miRs not found
miRs.not <- setdiff(all.miRs, rownames(miRmR))

miRmR[1:10,1:10]


## todos os genes
entrez <- colnames(miRmR); substr(entrez, 1,1) <- ""
class(entrez)

universe <- entrez
length(universe)

miR.DE <- c(levels(miR.DE.A$genes)[(miR.DE.A$genes)],levels(miR.DE.B$genes)[(miR.DE.B$genes)])
substr(miR.DE[grep("has", miR.DE)],1,3) <- "hsa"
length(miR.DE)
miR.DE <- unlist(strsplit(miR.DE, "\\-\\d+$",perl=TRUE))
sustenido <- grep("#$", miR.DE)
miR.DE[sustenido] <- gsub(".$","",miR.DE[sustenido])

miRmR.sel <- subset(miRmR, rownames(miRmR) %in% miR.DE)
dim(miRmR.sel)

setdiff(miR.DE, rownames(miRmR.sel)) ## microRNAs DE não presentes no banco ou com nomes incompativeis

miRmR.sel[1:10, 1:10]

length(which(apply(miRmR.sel, 2, sum) == 0))

miRmR.sel <- miRmR.sel[, -which(apply(miRmR.sel, 2, sum) == 0)]

dim(miRmR.sel)

sig.genes <- colnames(miRmR.sel); substr(sig.genes, 1,1) <- ""
length(sig.genes)



miRmR.sel.bin <- miRmR.sel
position <- which(miRmR.sel.bin>0, arr.ind = TRUE)
for(i in 1:nrow(position)){
	miRmR.sel.bin[position[i,1], position[i,2]] <- 1	
}

# Obter miRs ups e downs
dim(miR.DE.A)
miR.up <- c(levels(miR.DE.A[which(miR.DE.A$FC>1), "genes"])[miR.DE.A[which(miR.DE.A$FC>1), "genes"]], 
		levels(miR.DE.B[which(miR.DE.B$FC>1), "genes"])[miR.DE.B[which(miR.DE.B$FC>1), "genes"]])
substr(miR.up[grep("has", miR.up)],1,3) <- "hsa"
miR.down <- c(levels(miR.DE.A[which(miR.DE.A$FC<1), "genes"])[miR.DE.A[which(miR.DE.A$FC<1), "genes"]], 
		levels(miR.DE.B[which(miR.DE.B$FC<1), "genes"])[miR.DE.B[which(miR.DE.B$FC<1), "genes"]])
substr(miR.down[grep("has", miR.down)],1,3) <- "hsa"

correctNames <- function(x){
	x <- unlist(strsplit(x, "\\-\\d+$",perl=TRUE))
	sustenido <- grep("#$", x)
	x[sustenido] <- gsub(".$","",x[sustenido])
	x
}

miR.up <- correctNames(miR.up)
miR.down <- correctNames(miR.down)

length(miR.down)
length(which((row.names(miRmR.sel.bin) %in% miR.down)))
miR.down[which((miR.down %in% row.names(miRmR.sel.bin)) == FALSE)] ## MicroRNAs não encontrados 

length(miR.up)
length(which((row.names(miRmR.sel.bin) %in% miR.up)))
miR.up[which((miR.up %in% row.names(miRmR.sel.bin)) == FALSE)] ## MicroRNAs não encontrados

miRmR.sel.bin.up <- miRmR.sel.bin[which(row.names(miRmR.sel.bin) %in% miR.up),]
miRmR.sel.bin.up <- miRmR.sel.bin.up[, -which(apply(miRmR.sel.bin.up, 2, sum) == 0)]
miRmR.sel.bin.down <- miRmR.sel.bin[which(row.names(miRmR.sel.bin) %in% miR.down),]
miRmR.sel.bin.down <- miRmR.sel.bin.down[, -which(apply(miRmR.sel.bin.down, 2, sum) == 0)]
dim(miRmR.sel.bin.up)
dim(miRmR.sel.bin.down)

miRmR.sel.bin.neg <- miRmR.sel.bin 
miRmR.sel.bin.neg[which((row.names(miRmR.sel.bin) %in% miR.down)),] <- -1*miRmR.sel.bin.neg[which((row.names(miRmR.sel.bin) %in% miR.down)),]


genes.sum <- apply(miRmR.sel.bin, 2, sum)
genes.sum <- genes.sum[order(genes.sum, decreasing=T)]

genes.sum.neg <- apply(miRmR.sel.bin.neg, 2, sum)
genes.sum.neg <- genes.sum.neg[order(genes.sum.neg, decreasing=T)]

dif <- merge(genes.sum,abs(genes.sum.neg), by="row.names")
head(dif)
dif$dif <- apply(dif, 1, function(x){as.numeric(x["x"])-as.numeric(x["y"])})
head(dif)
summary(dif)

genes.sum.up <- apply(miRmR.sel.bin.up, 2, sum)
genes.sum.up <- genes.sum.up[order(genes.sum.up, decreasing=T)]
head(genes.sum.up)

genes.sum.down <- apply(miRmR.sel.bin.down, 2, sum)
genes.sum.down <- genes.sum.down[order(genes.sum.down, decreasing=T)]


## Genes atacados simultaneamente por micros Up e Down
length(which(dif$dif != 0))
dim(dif)
dif[which(dif$dif != 0),"Row.names"]
## Remover esses genes da lista down, mantendo na up
genes.sum.down <- genes.sum.down[(setdiff(names(genes.sum.down),dif[which(dif$dif != 0),"Row.names"]))]
length(genes.sum.down)

##### Write Tables
write.table(dif[which(dif$dif != 0),"Row.names"],paste("results/tables/miRmR/","genes.miRUpDown.txt"), quote=F)
write.table(miRmR.sel.bin.up, paste("results/tables/miRmR/","miRmR.up.txt",sep=""), quote=F, sep="\t")
write.table(miRmR.sel.bin.down, paste("results/tables/miRmR/","miRmR.down.txt",sep=""), quote=F, sep="\t")

write.table(genes.sum.down, paste("results/tables/miRmR/","genes.sum.down.txt",sep=""), quote=F, sep="\t")
write.table(genes.sum.up, paste("results/tables/miRmR/","genes.sum.up.txt",sep=""), quote=F, sep="\t")

genes.aux <- names(genes.sum.down)
substr(genes.aux, 1, 1) <- ""
head(genes.aux)
genes.sum.down.symbol <- genes.sum.down
names(genes.sum.down.symbol) <- genes.aux
names(genes.sum.down.symbol) <- as.character(unlist(SYMBOL[genes.aux]))
##rbind(names(genes.sum.down),names(genes.sum.down.symbol))[,1:10]	# verify
write.table(genes.sum.down.symbol, paste("results/tables/miRmR/","genes.sum.down.symbol.txt",sep=""), quote=F, sep="\t")

genes.aux <- names(genes.sum.up)
substr(genes.aux, 1, 1) <- ""
head(genes.aux)
genes.sum.up.symbol <- genes.sum.up
names(genes.sum.up.symbol) <- genes.aux
names(genes.sum.up.symbol) <- as.character(unlist(SYMBOL[genes.aux]))
##rbind(names(genes.sum.down),names(genes.sum.down.symbol))[,1:10]	# verify
write.table(genes.sum.up.symbol, paste("results/tables/miRmR/","genes.sum.up.symbol.txt",sep=""), quote=F, sep="\t")


genes.sum.up.symbol.gt3 <- names(which(genes.sum.up.symbol>=3))
write.table(genes.sum.up.symbol.gt3, paste("results/tables/miRmR/","genes.sum.up.symbol.gt3.txt",sep=""), 
		quote=F, sep="\t", col.names=F, row.names=F)
genes.sum.down.symbol.gt3 <- names(which(genes.sum.down.symbol>=3))
write.table(genes.sum.down.symbol.gt3, paste("results/tables/miRmR/","genes.sum.down.symbol.gt3.txt",sep=""),
		quote=F, sep="\t", col.names=F, row.names=F)


## 






#library(ReactomePA)
#
#
### Get the gene names that are mapped to an entrez gene identifier
#mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
### Convert to a list
#SYMBOL <- as.list(org.Hs.egSYMBOL[mapped_genes])
#
#
######UP
#my.up <- names(which(genes.sum.up>3))
#substr(my.up, 1, 1) <- ""
#x.up <- enrichPathway(gene=my.up, 
#			#organism="human",
#			pvalueCutoff=0.05, 
#			#qvalueCutoff=0.05, 
#			readable=T)
#summary(x.up)
#plot(x.up, showCategory=5)
#
##### DOWN
#my.down <- names(which(genes.sum.down>3))
#substr(my.up, 1, 1) <- ""
#x.down <- enrichPathway(gene=my.down, organism="human",pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
#summary(x.down)
#plot(x.down, showCategory=5)
#
#
#
##library(GOstats)
##
##GOparams <- new("GOHyperGParams", geneIds=sig.genes, universeGeneIds=universe,
## annotation=annpkg, ontology="BP", pvalueCutoff=0.001, conditional=TRUE, testDirection="over")
##
##hgOver <- hyperGTest(GOparams)
##class(hgOver)
##summary(hgOver)
#
##saveArray('saves/')
##loadArray('saves/')
#
#
#
#library("clusterProfiler")
#GO.BP.up <- enrichGO(gene=my.up,
#		organism="human",
#		ont="BP",
#		pvalueCutoff = 0.01,
#		qvalueCutoff = 0.05,
#		readable=TRUE)
#head(summary(GO.BP.up))
#class(GO.BP.up)
#
#plot(GO.BP.up)
#plot(GO.BP.up, type="cnet")


#### Grafico IPA Canonical Pathway
getwd()
canpath.up <- read.table("IPAanalysis/pathway.pvalue.up", header=F,dec=",")
canpath.down <- read.table("IPAanalysis/pathway.pvalue.down", header=F, sep="\t", dec=",")

#canpath.up <- canpath.up[order(canpath.up$V2),]
#rownames(canpath.up) <- 1:nrow(canpath.up) 
#canpath.down <- canpath.down[order(canpath.down$V2),]
#rownames(canpath.dowm) <- 1:nrow(canpath.down)
#
##head(canpath.up)
##canpath.up <- transform(canpath.up, V1 = reorder(V1, order(V2)))
##
##canpath.up$V1 <- factor(canpath.up$V1, levels=canpath.up[order(canpath.up$V1, -canpath.up$V2), ]$V1) #reorder by grp/value
canpath.up$V1 <- levels(canpath.up$V1)[canpath.up$V1]

png(paste("IPAanalysis/CanPath.up","png", sep="."), res =300, height=2000, width=2000) 
p <- ggplot(canpath.up, aes(y=V2, x=V1, 
     ymin = 0, ymax = V2)) + geom_point() + coord_flip() 
p + geom_linerange() + scale_x_discrete(limits=canpath.up$V1[order(canpath.up$V2, decreasing = T)[1:50]])+
		xlab("Canonical Pathways") + ylab("-log(p-value)")
dev.off()

png(paste("IPAanalysis/CanPath.down","png", sep="."), res =300, height=2000, width=2000)
p <- ggplot(canpath.down, aes(V1, V2, 
				ymin = 0, ymax = V2)) + geom_point() + coord_flip() 
p + geom_linerange() + 
		scale_x_discrete(limits=canpath.down$V1[order(canpath.down$V2, decreasing = T)[1:50]])+
		xlab("Canonical Pathways") + ylab("-log(p-value)")

dev.off()


func.down <- read.table("IPAanalysis/MB_projDown_down_g3_functions.txt", header=T, sep="\t", dec=",")
png(paste("IPAanalysis/Functions.down","png", sep="."), res =300, height=2000, width=2000)
qplot(reorder(factor(Category),factor(Category),length),data=subset(func.down, p.Value<=0.01),geom="bar", xlab="Functions", ylab= "Category")
last_plot() + coord_flip()
dev.off()

func.up <- read.table("IPAanalysis/MB_projDown_up_g3_functions.txt", header=T, sep="\t", dec=",")
png(paste("IPAanalysis/Functions.up","png", sep="."), res =300, height=2000, width=2000)
qplot(reorder(factor(Category),factor(Category),length),data=subset(func.up, p.Value<=0.01),geom="bar", xlab="Functions", ylab= "Category")
last_plot() + coord_flip()
dev.off()





#####
## Data 04/07/12
## respondendo requisição joice 3 de julho de 2012 12:28
## duvidas em relação a correção do p-value
head(qDE.my.limma.A[['D-C']])
hist(qDE.my.limma.A[['D-C']][, 'p.value'])

#x11()


png("results/plots/pvalues.png", 1400,1000)
par(mfrow=c(4,4))
for(j in c('A','B')){
	for(i in p.adjust.methods){
		hist(p.adjust(get(paste("qDE.my.limma.",j, sep=""))[['D-C']][, 'p.value'], method=i), breaks=100,
				xlab="Adjusted p-value", main=paste("Set ",j," -",i))
		abline(v=0.05, col="red")
	}
}
dev.off()





#####
## Data 16/07/12
## respondendo requisição joice 16 de julho de 2012 
## alvos dos microRNAs com p-valor ajustado suficiente

loadArray('saves/')

## Tabela de associação miRmR
miRmR[1:10, 1:10]

## entrez ids
entrez <- colnames(miRmR); substr(entrez, 1,1) <- ""

## obter novo set de microRNAs diferencialmente utilizando o p-value ajustado
miR.DE.A.adj <- subset(my.DE.A, adj.p.value<pv.cut & ((FC >= fc.cut) | (FC <= 1/fc.cut)) )
write.table(miR.DE.A.adj$genes, file=paste("results","tables",paste("sig.miR.A.adj.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=F, quote=F)
write.table(miR.DE.A.adj, file=paste("results","tables",paste("sig.miR.values.A.adj.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=T, quote=F, sep="\t")

miR.DE.B.adj <- subset(my.DE.B, adj.p.value<pv.cut & ((FC >= fc.cut) | (FC <= 1/fc.cut)) )
write.table(miR.DE.B.adj$genes, file=paste("results","tables",paste("sig.miR.B.adj.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=F, quote=F)
write.table(miR.DE.B.adj, file=paste("results","tables",paste("sig.miR.values.B.adj.p_",pv.cut,".txt",sep=""),sep="/"), row.names=F, col.names=T, quote=F, sep="\t")


selected.miR <- c(levels(miR.DE.B.adj$genes)[miR.DE.B.adj$genes], levels(miR.DE.A.adj$genes)[miR.DE.A.adj$genes])
selected.miR <- unlist(strsplit(selected.miR, "\\-\\d+$",perl=TRUE))
sustenido <- grep("#$", selected.miR)
selected.miR[sustenido] <- gsub(".$","",selected.miR[sustenido])


new.rel <- miRmR[which(rownames(miRmR) %in% selected.miR), ]	## subset
col.sum <- apply(new.rel, 2, sum)
new.rel <- miRmR[which(rownames(miRmR) %in% selected.miR), which(col.sum!=0)]


##convert xgeneid em symbols
my.genes <- colnames(new.rel)
substr(my.genes, 1, 1) <- ''

##geneid discontinued 115669 - new geneid 319089
my.genes[which(my.genes == '115669')] <- '319089'

## Data structure with Symbol names for each entrez id
x <- org.Hs.egSYMBOL
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x)
# Convert to a list
SYMBOL <- as.list(x[mapped_genes])
my.genes <- unlist(SYMBOL[as.character(my.genes)])

length(my.genes)
length(colnames(new.rel))


colnames(new.rel) <- my.genes

new.rel[1:6, 1:10]

apply(new.rel, 1, sum)


class(new.rel)

atleast <- 3
length(unique(subset(melt(as.matrix(new.rel)), value >=  atleast )[,2]))	#nro de genes
length(unique(subset(melt(as.matrix(new.rel)), value >=  atleast )[,1]))	#nro de micros
write.table(subset(melt(as.matrix(new.rel)), value >=  atleast ), file=paste('results/tables/miRmR/','miRmR_bd_atleast',atleast,'.txt',sep=""), quote=F, sep="\t", row.names=F, col.names=F)
aux<-subset(melt(as.matrix(new.rel)), value >=  atleast ); write.xls(aux, file=paste('results/tables/miRmR/','miRmR_bd_atleast',atleast,'.xls',sep=""), row.names=F, col.names=F)


saveArray('saves/')


##############

###  EMAIL ERICA  3 de agosto de 2012 15:34
#Olá Wilson, tudo bem?
#		
#Liguei hj para vc, mas não consegui te encontrar.
#
#Estou “em desespero” em função do prazo que temos para fazer a análise funcional,
#que depende dos resultados da análise de bioinformática. A FAPESP só prorrogou 
#meu auxílio por 05 meses e, como não temos experiência em análise funcional e tb 
#temos que comprar todo o material importado, temos pouco tempo.
#Portanto, precisamos dos dados com urgência. Temos Skype e, como vc sugeriu, 
#podemos marcar um horário para conversarmos na segunda?
#		
#Bom final de semana!
#ÉRIKA.



################

### Sessão Skype com Joyce e Jorge 13/08/12
## Com base nos ~3500 genes definidos como alvos dos 6 micros DE
## buscaremos mais microRNAs que atuam sobre os genes a fim de
## aumentar a confiabilidade sobre os genes afetados, devemos também
## avaliar os micros utilizando o TARBASE

loadArray('saves/')

## Ler alvos

targets <- as.character(unique(read.table(paste("alvosAfetados/targets_symb.csv",sep=""), sep=",", header=T)[,"mR"]))

## Utilizar RmiR.mR

mirbase <- unique(dbReadTable(RmiR.Hs.miRNA_dbconn(), "mirbase")[,1:2])
dim(mirbase)
mirbase <- mirbase[grep('hsa',mirbase$mature_miRNA),]
dim(mirbase)

allRel <- which(mirbase$gene_id %in% targets)

## table(mirbase[allRel,"gene_id"]) # numero de micros interagindo com um mesmo gene



correctNames <- function(x){
	x <- unlist(strsplit(x, "\\-\\d+$",perl=TRUE))
	sustenido <- grep("#$", x)
	x[sustenido] <- gsub(".$","*",x[sustenido])
	x
}

eqFiltA.aux <- eqFiltA.2
rownames(eqFiltA.aux) <- correctNames(rownames(eqFiltA.aux))
eqFiltB.aux <- eqFiltB.2
rownames(eqFiltB.aux) <- correctNames(rownames(eqFiltB.aux))

mean.ctrl <- 17-c(apply(eqFiltA.aux[,control],1, mean, na.rm=T), apply(eqFiltB.2[,control],1, mean, na.rm=T)) ########### OBSERVE!!! Valores mais altos equivalem a maior expressão pois esta sendo utilizado 17-Valor
mean.down <- 17-c(apply(eqFiltA.aux[,down],1, mean, na.rm=T), apply(eqFiltB.2[,down],1, mean, na.rm=T))

summary(mean.ctrl)


#tarbase <- read.table("/data/matheus/TCC/targetPrediction/tarbase/v5.0/TarBase_V5.0.human.tab", header=T)
#tarbase$Ensembl <- levels(tarbase$Ensembl)[tarbase$Ensembl]
tarbase <- read.table("/data/projects/bit_data/down/htqpcr/miRNA/tarbase/v6.0/tarbase.rel", header=F, colClasses=c("character", "character"))
colnames(tarbase) <- c("entrez", "miRNA")
head(tarbase)
#dim(tarbase)


# Get Ensembl Gene ID
#library(biomaRt)
#ensembl <- useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl"))
#attr <- listAttributes(ensembl); attr[grep("entrez",attr[,"description"], ignore.case=T),]
#result <- getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters=c("entrezgene"), targets, ensembl, uniqueRows = T)
#result$entrezgene <- as.character(result$entrezgene) 


result <- data.frame(entrezgene=targets, stringsAsFactors=FALSE)

#write.table(targets, "tarbase/v6.0/entrezid.txt", col.names=F, row.names=F, quote=F )

for(line in 1:nrow(result)){
	aux <- result[line, "entrezgene"]
	#idx <- which(result$entrezgene==aux)
	idx <- line
	rel <- which(mirbase$gene_id == aux)	# obter todas as relacoes de um micro 
	result[idx, "miR"] <- length(rel)
	
	micros.aux.ctrl <- (na.omit(mean.ctrl[mirbase$mature_miRNA[rel]]))	# micros que regulam realmente algum gene
	result[idx, "miR.ctrl"] <- length(micros.aux.ctrl)
	result[idx, "mean.ctrl"] <- mean(micros.aux.ctrl)
	result[idx, "sum.ctrl"] <- sum(micros.aux.ctrl)
	
	micros.aux.down <- (na.omit(mean.down[mirbase$mature_miRNA[rel]]))	# micros que regulam realmente algum gene
	result[idx, "miR.down"] <- length(micros.aux.down)
	result[idx, "mean.down"] <- mean(micros.aux.down)
	result[idx, "sum.down"] <- sum(micros.aux.down)
	result[idx, "miR.intersect"] <- length(intersect(names(micros.aux.down),names(micros.aux.ctrl)))
	
	## Verifying tarbase
	#names(micros.aux.ctrl)
	#names(micros.aux.down)
	
	#inters <- which(tarbase$entrez %in% result[line, "entrezgene"])
	inters <- which(tarbase$entrez %in% result[line, "entrezgene"])
	if(length(inters) > 0){
		#print(tarbase[inters,])
		#grep(tarbase[inters,"miRNA"], names(micros.aux.ctrl))
		#grep(tarbase[inters,"miRNA"], names(micros.aux.down))
		tarbase.down <- names(micros.aux.down)[unlist(lapply(tarbase[inters,"miRNA"], function(x){ na.exclude(match(toupper(x),toupper(names(micros.aux.down))))}))]
		tarbase.ctrl <- names(micros.aux.ctrl)[unlist(lapply(tarbase[inters,"miRNA"], function(x){ na.exclude(match(toupper(x),toupper(names(micros.aux.ctrl))))}))]
		#tarbase.ctrl <- print(names(micros.aux.ctrl)[unlist(lapply(tarbase[inters,"miRNA"], function(x){grep(x, names(micros.aux.ctrl))}))])
		result[line, "tarbase.down"] <- length(tarbase.down)
		if(length(tarbase.down) > 0)
			result[line, "tarbase.down.names"] <- paste(tarbase.down, sep=",", collapse=",")
		else
			result[line, "tarbase.down.names"] <- NA
		result[line, "tarbase.ctrl"] <- length(tarbase.ctrl)
		if(length(tarbase.ctrl) > 0)
			result[line, "tarbase.ctrl.names"] <- paste(tarbase.ctrl, sep=",", collapse=",")
		else
			result[line, "tarbase.ctrl.names"] <- NA
		#break
	}

}


aux <- subset(result,tarbase.ctrl > 0 | tarbase.down > 0)
colnames(result)

# Realiza fold-change sobre soma dos valores de expressão dos miRNAs
result$fold <- result$sum.down / result$sum.ctrl
result[which(result$fold < 1), "fold.change"] <- -1/result[which(result$fold < 1), "fold"] 

## Obtem gene symbols
aux <- getBM(attributes=c("hgnc_symbol", "entrezgene"),filters=c("entrezgene"), result$entrezgene, ensembl, uniqueRows = T)
result <- merge(result, aux, by.x = "entrezgene")
result <- result[,c("entrezgene", "hgnc_symbol", "miR", "miR.ctrl", "mean.ctrl", "sum.ctrl", "miR.down", "mean.down", "sum.down", "miR.intersect", "tarbase.down", "tarbase.down.names", "tarbase.ctrl",
		"tarbase.ctrl.names", "fold", "fold.change")]        

write.table(result, "alvosAfetados/fold_tarbase/result.txt", col.names=T, quote=F, row.name=F, sep="\t")

# corte utilizando fold.change 
result.fc1.5 <- subset(result, abs(fold.change) > 1.5 )
write.table(result.fc1.5[order(result.fc1.5$fold.change), c("hgnc_symbol", "entrezgene", "sum.down", "sum.ctrl", "fold.change")], "alvosAfetados/fold_tarbase/result.fc1.5.txt", col.names=T, quote=F, row.name=F, sep="\t")

#saveArray('saves/')
loadArray('saves/')



##load('/data/projects/bit_data/down/htqpcr/mRNA/saves/R_03_09_2012_2.RData')
#de.genes
#
#inter <- data.frame()
#for(myg in de.genes){
#	res <- result[grep(myg, result$hgnc_symbol),]
#	if (ncol(res) > 0){
#		inter <- rbind(inter, res)
#	}
#}
#
#inter
#
#selected.miR


### Ler valores de expressao genica

mRNA.data <- read.table("/data/projects/bit_data/down/htqpcr/mRNA/saves/tables/data.txt", header=T)
gene.names <- unlist(lapply(strsplit(rownames(mRNA.data), "-"), function(x) x[1]))

x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx[["10"]]
(unlist(unique(xx[result$entrezgene])))

## Avaliar intersecção entre result.fc1.5 e genes

head(result)


dim(result[order(result$entrezgene), c("entrezgene", "hgnc_symbol")])
dim(unique(result[order(result$entrezgene), c("entrezgene", "hgnc_symbol")]))
length(unique(result$hgnc_symbol))


miRintersect <-  unlist(unique(xx[result$entrezgene]))[ (which((unlist(unique(xx[result$entrezgene]))) %in% gene.names)) ]

write.table(miRintersect, '/data/projects/bit_data/down/htqpcr/intersect/genes.txt', col.names=F, row.names=F, quote=F)



#write.table(unlist(unique(xx[result$entrezgene])) , "/tmp/symbol.fmiRNA", quote=F, row.names=F, col.names = F)
#write.table(unique(result$hgnc_symbol), "/tmp/symbol.fmiRNA", quote=F, row.names=F, col.names = F)
#write.table(unique(gene.names), "/tmp/symbol.fmRNA", quote=F, row.names=F, col.names = F)


#c("TNFRSF18", "VEGFA") %in% gene.names

#inter <- data.frame()
#for(myg in gene.names){
#	res <- result[grep(myg, result$hgnc_symbol),]
#	if (ncol(res) > 0){
#		inter <- rbind(inter, res)
#	}
#}
#
#inter







