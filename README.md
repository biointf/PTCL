# PTCL
Peripheral T-cell Non-Hodgkin Lymphoma  classification

---
title: "PTCL classification"
author:
- affiliation: University of Milan
  name: Luca Agnelli
- affiliation: Wellcome Trust Sanger Institute, UK and University of Milan, IT
  name: Francesco Maura
date: "08/09/2017"
output:
  html_notebook: default
  html_document: default
---
Built with R version:  
`r getRversion()`

# Libraries

Load necessary libraries

```{r setup, include=TRUE}
library(affy)
library(AnnotationDbi)
library(biomaRt)
library(circlize)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(glmnet)
library(gplots)
library(heatmap.plus)
library(Hmisc)
library(lattice)
library(limma)
library(MASS)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(rgl)
library(rtracklayer)
library(survival)
library(VennDiagram)

set1 = c(brewer.pal(9,"Set1"), brewer.pal(8, "Dark2"))

violinJitter <- function(x, magnitude=1){
  d <- density(x)
  data.frame(x=x, y=runif(length(x),-magnitude/2, magnitude/2) * approxfun(d$x, d$y)(x))
}

rotatedLabel <- function(x0 = seq_along(labels), y0 = rep(par("usr")[3], length(labels)), labels, pos = 1, cex=1, srt=45, ...) {
  w <- strwidth(labels, units="user", cex=cex)
  h <- strheight(labels, units="user",cex=cex)
  u <- par('usr')
  p <- par('plt')
  f <- par("fin")
  xpd <- par("xpd")
  par(xpd=NA)
  text(x=x0 + ifelse(pos==1, -1,1) * w/2*cos(srt/360*2*base::pi), y = y0 + ifelse(pos==1, -1,1) * w/2 *sin(srt/360*2*base::pi) * (u[4]-u[3])/(u[2]-u[1]) / (p[4]-p[3]) * (p[2]-p[1])* f[1]/f[2] , labels, las=2, cex=cex, pos=pos, adj=1, srt=srt,...)
  par(xpd=xpd)
}

avefc = function (y, log=TRUE, replace= FALSE) {
	 if (log) y = 2^y
   if (replace) y = y + (1-min(y))
   m = apply(y,1,mean)
	 y.n = y/m
	 y.n2 = y.n
	 y.n2 [y.n2 < 1] = 1/ (y.n2 [y.n2 < 1])
	 ave.fc = apply (y.n2, 1, mean)
	 return(ave.fc)
}

rle.custom = function (a, logged2 = TRUE, file = NULL, colorbox= NULL, labels=NULL , legend = NULL ) {
	a.m <- apply(a,1,median)
if (logged2) {
	for (i in 1:dim(a)[2]) {
   		 a [,i] <-  a [,i] - a.m
	}
	} else {
		for (i in 1:dim(a)[2]) {
   		 a [,i] <-  log (a [,i] / a.m )
	}
	}
	png(file,10240,3840)
	par(mar=c(10,4,6,2))
	boxplot (a, ylim= c(-5,5), outline=F, col=colorbox, xlab="pts", names=labels, las=2, cex.axis = 1.5, main="RLE", xlim = c(1,600), cex.main = 5 )
	legend("bottomright",legend = c(levels(as.factor(pts.info.data$final.molec))),   
	  fill = colorz, # 6:1 reorders so legend order matches graph
 	  title = "Legend",
 	  cex = 5)
	#abline(h=0.5, lty=2, col="grey30"); abline(h=-0.5, lty=2, col="grey30")
	dev.off()

	a.c = apply(a, 2, stats::quantile)
	return(a.c)
}

```

# Define Ensembl library for  gene convertion from array to HUGO

```{r}
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
```

# Gene expression data: upload or generate GEP normalized matrix

```{r}
### choice 1: import processed matrix
data.dir="~/Dropbox/on.going.papers/ptcl-nos.INT/paper/script/" # set custom data.dir

setwd(data.dir)
load (file.path(data.dir,"/Rmd.files/541_PTCL_batch_adjusted_geo.id.Rdata"))

geneExpr = adj.data
# import batch and re-order accordingly
load(file.path(data.dir,"/Rmd.files/PTCL.batch.Rdata"))
batch = batch [order(batch$nameNEW),]
batch.series = as.vector(batch$center)

### end of choice 1

### choice 2: generate your own affy object and custom data

### # download CEL files from GEO series GSE6338, GSE19067, GSE19069, GSE40160, GSE58445, GSE65823 and EBI series ETABM702, ETABM783
### # GSM368580.CEL, GSM368582.CEL, GSM368584.CEL, GSM368586.CEL, GSM368589.CEL, GSM368591.CEL, GSM368594.CEL, GSM472164.CEL, GSM1411278.CEL, GSM1411284.CEL, GSM1411285.CEL, GSM1411287.CEL, GSM1411355.CEL, GSM1411364.CEL, GSM1411368.CEL, GSM1411425.CEL, GSM1411427.CEL excluded from the analysis (see Methods for explaination")
### celfiles <- dir("~/Documents/DATI/PTCL.nos/GSE6338-GSE19067-GSE19069-GSE40160-GSE58445-GSE65823-ETABM702-ETABM783/", pattern = ".CEL")
### library(hgu133plus2hsentrezgcdf)
### library(affy)
### gset = justRMA(celfile.path = "/Users/emagene/Documents/DATI/PTCL.nos/GSE6338-GSE19067-GSE19069-GSE40160-GSE58445-GSE65823-ETABM702-ETABM783/", ### filenames = celfiles, sampleNames = gsub(".CEL","", celfiles), cdfname = "hgu133plus2hsentrezgcdf")
### geneExpr = exprs(gset)
### batch adjustment
### library(sva)  
### # import batch and re-order accordingly
### load("./Rmd.files/PTCL.batch.Rdata")
### batch = batch [order(rownames(batch)),]
### batch.series = as.vector(batch$center)
### geneExprNEW = geneExpr [ , order(colnames(geneExpr)) ]
### geneExprNEW = geneExprNEW[grep("AFFX",rownames(geneExprNEW), invert=TRUE),]
### # check order correspondence and, if correct, adjust data
### if (all(colnames(geneExprNEW) == rownames(batch))) {
###   adj.data = ComBat (geneExprNEW, batch.series, mod = NULL, par.prior = TRUE, prior.plots = TRUE)
### } else {
###   cat("Error: colnames and batch did not correspond")
### }
### geneExpr = adj.data
### colnames(geneExpr) = as.vector(batch$nameNEW)
### end of choice 2
```

# Clinical Data

```{r}
pts.info.data <- read.table("./Rmd.files/541_paz_info_MUT.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors = F)

# customize colors for categories
levels(as.factor(pts.info.data$final.molec))
# "AITL"     "ALCL.neg" "ALCL.pos" "ATLL"     "NKT"      "PTCL.nos" "T.CD30"   "T.CD4"    "T.CD8"    "T.DR"     "T.reg"    "TCR-HL"  
colorz = c("black", "yellow","dodgerblue2","brown2","forestgreen", "orange", "grey42", "grey52","grey62","grey72","grey82","grey92")
temp = split (  pts.info.data$sample.nameNEW, pts.info.data$final.molec )
colorx = colnames(geneExpr)
length(colorz)
length(temp)
for (i in 1:length(colorz)) colorx [ which(colorx %in% unlist(temp[i])) ] = colorz[i]
library(gplots)
colorx = col2hex(colorx)
```

# Pie Chart with Percentages

```{r}
slices <- table(pts.info.data$final.molec)
lbls <- names(table(pts.info.data$final.molec))
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, ": ", slices, " (", pct, "%)", sep="" ) # add percents to labels
par(mfrow=c(1,1))
par(mar=c(6,6,6,6), xpd=F)
pie(slices,labels = lbls, init.angle = 0, col=colorz,  main="", cex=0.8, radius=1)
```

# PCA

```{r}
# apply variational filter

afc2 = avefc (geneExpr, log=TRUE, replace=FALSE)
data541exprs.vf = geneExpr [afc2 >= 2, ]
dim(data541exprs.vf	)
# retry PCA on shorted gene list
data541m = t(as.matrix(data541exprs.vf))
pca<-prcomp(data541m,scale=T)
plot3d(pca$x,rgl.use=T,col=colorx,size=0.6,type="s")
# texts3d(pca$x,text=pts.info.data$final.molec,adj=-0.5,cex=1) # run to add labels to samples

```

# Heatmap

```{r}
mat = as.matrix(data541exprs.vf)
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))
types = pts.info.data$final.molec
color.annot = col2hex(colorz); names (color.annot)= names(temp)  
ha = HeatmapAnnotation(df = data.frame(type = types) , col = list(type = c( color.annot ) ) )
ha@anno_list[[1]]@color_mapping@colors = col2hex(colorz)
names(ha@anno_list[[1]]@color_mapping@colors) = names(temp)

ht = Heatmap(mat_scaled, name = "expression", km = 7, clustering_method_columns = "ward.D", col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")), top_annotation = ha, top_annotation_height = unit(4, "mm"), show_row_names = FALSE, show_column_names = FALSE)

column_order(ht)
```


# Check relative log expression after batch correction

```{r}
rle.medians = rle.custom(geneExpr, colorbox=colorx, file="./RLE.541.png", labels=pts.info.data$sample.nameNEW )
plot(rle.medians[3,], type="l", xlab="pts", ylab="RLE median", ylim= c(-1,1) )
abline(h=0.5, lty=2, col="grey50"); abline(h=-0.5, lty=2, col="grey50")
```


# Final gene expression matrix: define design file and filter geneExpr for patients included in design data frame

```{r}
design <- pts.info.data[,c(1:2,6:8,14:17)]
rownames(design)<- design[,1]
design<- design[,-c(1:2)]
design<-na.omit(design) ### select onyl patients with all mutations data available (n=53)
design$age<- as.numeric(as.character(design$age))
design$age<- design$age - median(design$age)
design[design == "WT"] <- 0
design[design == "MUT"] <- 1
design$final.molec[design$final.molec=="AITL"] <- 0
design$final.molec[design$final.molec=="PTCL.nos"] <- 1
design$gender[design$gender=="M"] <- 1
design$gender[design$gender=="F"] <- 0
design$offset <- rep(1, nrow(design))
design<-design[,c(8,1:7)]

all(pts.info.data$sample.nameNEW == colnames(geneExpr)) ## check correspondence
# geneExpr = geneExpr [ , order (pts.info.data$geo.id)] ### do only to set correspondence in case of custom procedure
# colnames(geneExpr) = pts.info.data$sample.nameNEW [ order (pts.info.data$geo.id)]

geneExpr2<- (geneExpr[, rownames(design)])
geneExpr2<- data.matrix(geneExpr2, rownames.force = NA)
design<- data.matrix(design, rownames.force = NA)
```

# Model fitting

```{r}
glm = lmFit(geneExpr2[,rownames(design)], design = design )
glm = eBayes(glm)
F.stat <- classifyTestsF(glm[,-1],fstat.only=TRUE)
glm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){
  glm$F.p.value <- pchisq(df1*glm$F,df1,lower.tail=FALSE)
}else
  glm$F.p.value <- pf(glm$F,df1,df2,lower.tail=FALSE)

set.seed(12345678)
rlm <- lmFit(geneExpr[,rownames(design)], apply(design, 2, sample))
rlm <- eBayes(rlm)
F.stat <- classifyTestsF(rlm[,-1],fstat.only=TRUE)
rlm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){
  rlm$F.p.value <- pchisq(df1*rlm$F,df1,lower.tail=FALSE)
}else
  rlm$F.p.value <- pf(rlm$F,df1,df2,lower.tail=FALSE)
F.stat <- classifyTestsF(glm[,2:5],fstat.only=TRUE)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
F.p.value <- pchisq(df1*F.stat,df1,lower.tail=FALSE)
R.stat <- classifyTestsF(rlm[,2:5],fstat.only=TRUE)
Rall = 1 - 1/(1 + glm$F * (ncol(design)-1)/(nrow(design)-ncol(design)))
Rgenetics = 1 - 1/(1 + F.stat * 4/(nrow(design)-ncol(design)))
Pgenetics = 1 - 1/(1 + R.stat * 4/(nrow(design)-ncol(design)))
names(Rgenetics) <- names(Pgenetics) <- names(Rall) <-  rownames(geneExpr)

```

# Differentially Expressed Genes

```{r}
par(bty="n", mgp = c(2,.33,0), mar=c(3,2.5,1,1)+.1, las=1, tcl=-.25, xpd=NA)
d <- density(Pgenetics,bw=1e-3)
f <- 40 #nrow(gexpr)/512
par(mfrow=c(1,1))
plot(d$x, d$y * f, col='grey', xlab=expression(paste("Explained variance per gene ", R^2)), main="", lwd=2, type="l", ylab="", xlim=c(0,1.2))
title(ylab="Density", line=1.5)
d <- density(Rgenetics, bw=1e-3)
r <- min(Rgenetics[p.adjust(F.p.value,"BH")<0.01]) ######## threshold to select 412 genes
x0 <- which(d$x>r)
polygon(d$x[c(x0[1],x0)], c(0,d$y[x0])* f, col=paste(set1[1],"44",sep=""), border=NA)
lines(d$x, d$y* f, col=set1[1], lwd=2)
text(d$x[x0[1]], d$y[x0[1]]*f, pos=4, paste(sum(Rgenetics > r), "genes q < 0.01"))
legend("topright", bty="n", col=c(set1[1], "grey"), lty=1, c("Observed","Random"), lwd=2)


glmPrediction <- glm$coefficients %*% t(design)
rlmPrediction <- rlm$coefficients %*% t(design)
```

# Print signficiant genes

```{r}
kk<-as.data.frame((p.adjust(F.p.value,"BH")<0.01))
kk$gene<- rownames(kk)
colnames(kk)[1]<-"code"
kk2<-kk[kk$code=="TRUE",]
### sort(kk2$gene) ##### if you want to print the entire list of differentially expressed genes

```

# Extract the list of differentially expressed genes by mutations

```{r}
# customize colors in colMutations
colMutations = col2hex(c("magenta", "purple","gray60","red","lightblue","green","orange"))
names(colMutations) <- colnames(design)[-1]

gene_code<- kk2$gene
tab=NULL
for(i in (1:length(kk2$gene)))
{
  gene_single<- gene_code[i]
  y <- glm$coefficients[gene_single,-1]+glm$coefficients[gene_single,1]
  w <- glm$p.value[gene_single,-1] < 0.01
  int<-c(gene_single, as.character(w))
  tab<- rbind(tab, int)
}
rownames(tab)<-seq(1:nrow(tab))
colnames(tab)<- c("gene",colnames(design)[-1])

# Write to disk a file with all significant genes
write.table(tab, "table_differentially_expressed_gene.txt",sep="\t", quote=F, row.names = F, col.names = T)
```

# Example of extraction

```{r}
gene_single= "10413_at" ### YAP1
par(bty="n", mgp = c(1.5,.33,0), mar=c(2.5,2.5,1,1)+.1, las=1, tcl=-.25, xpd=F)
  temp_name = unique(getBM( attributes = c("ensembl_transcript_id", "entrezgene", "external_gene_name"), filters = "entrezgene", values = gsub("_at","",gene_single), mart = ensembl )$external_gene_name)
  plot(glmPrediction[gene_single,], geneExpr[gene_single,rownames(design)], ylab=(paste("Observed ",temp_name, " expression")), xlab=(paste("Predicted ",temp_name, " expression")), pch=16, cex=.8)
  abline(0,1)
  u <- par("usr")
  par(xpd=NA)
  y <- glm$coefficients[gene_single,-1]+glm$coefficients[gene_single,1]
  u <- par("usr")
  x0 <- rep(u[3]+1,ncol(design)-1)
  y0 <- u[4] + 0.05*(u[4]-u[3]) - rank(-y)/length(y) * (u[4]-u[3])/1.2
  d <- density(y)
  lines(d$x, d$y/5+1+u[3], col="grey")
  lines(d$x, -d$y/5+1+u[3], col="grey")
points(x=y, y=x0+violinJitter(y, magnitude=0.25)$y, col=colMutations, pch=16)
  text(x=glm$coefficients[gene_single,1], y= 1.5 +u[3], "Model coefficients", cex=0.8)
legend("topleft",names(colMutations), col = colMutations, bty= "n", cex = 0.8, pch = 16)
```

# Plot significant effects per covariate (q<0.01)

```{r}
testResults <- decideTests(glm, method="hierarchical",adjust.method="BH", p.value=0.01)[,-1]
significantGenes <- sapply(1:ncol(testResults), function(j){
  c <- glm$coefficients[testResults[,j]!=0,j+1]
  table(cut(c, breaks=c(-5,seq(-1.5,1.5,l=7),5)))
})

colnames(significantGenes) <- colnames(testResults)
rownames(tab)<-c(1:nrow(tab))
tab2<- as.data.frame(tab)
tab2$gene<-as.character(as.character(tab2$gene))
tab2$final.molec<-as.character(as.character(tab2$final.molec))
tab2$TET2<-as.character(as.character(tab2$TET2))
tab2$RHOA<-as.character(as.character(tab2$RHOA))
tab2$IDH2<-as.character(as.character(tab2$IDH2))
tab2$DNMT3A<-as.character(as.character(tab2$DNMT3A))

par(mfrow=c(1,1))
par(bty="n", mgp = c(2.5,.33,0), mar=c(5,5.5,5,0)+.1, las=2, tcl=-.25)
b <- barplot(significantGenes, las=2, ylab = "Differentially expressed genes", col=brewer.pal(8,"RdYlBu"), legend.text=FALSE , border=0, xaxt="n")
rotatedLabel(x0=b-0.1, y0=rep(-2, ncol(significantGenes)), labels=colnames(significantGenes), cex=.7, srt=45, font=ifelse(grepl("[[:lower:]]", colnames(design))[-1], 1,3), col=colMutations)
rotatedLabel(b-0.1, colSums(significantGenes), colSums(significantGenes), pos=3, cex=, srt=45)
clip(0,30,0,1000)
x0 <- 7.5
image(x=x0+c(0,0.8), y=par("usr")[4]+seq(-1,1,l=9) -4, z=matrix(1:8, ncol=8), col=brewer.pal(8,"RdYlBu"), add=TRUE)
text(x=x0+1.1, y=par("usr")[4]+c(-1,0,1) -4, format(seq(-1,1,l=3),2), cex=0.66)
lines(x=rep(x0+0.9, 2), y=par("usr")[4]+c(-1,1) -4)
segments(x0+0.9,par("usr")[4] + 1-4,x0+0.95,par("usr")[4] + 1-4)
segments(x0+0.9,par("usr")[4] + 0-4,x0+0.95,par("usr")[4] + 0-4)
segments(x0+0.9,par("usr")[4] + -1-4,x0+0.95,par("usr")[4] + -1-4)
text(x0 + 0.45, par("usr")[4] + 1.5-4, "log2 FC", cex=.66)

```

# Print the list of differentially expressed genes using the Ensembl annotation

```{r}
select_hist<- pts.info.data[pts.info.data$final.molec == "AITL" |  pts.info.data$final.molec == "PTCL.nos",]
gene<- as.data.frame(testResults)
sig_genes<- gene[gene$final.molec!= 0 |gene$IDH2 != 0 | gene$TET2 != 0 | gene$DNMT3A != 0 | gene$RHOA != 0,]
list_genes<-sort(rownames(sig_genes)) ##### list of signficiant genes

geneannotation1 <- getBM( attributes = c("ensembl_transcript_id", "entrezgene", "external_gene_name"), filters = "entrezgene", values = gsub("_at","",list_genes), mart = ensembl)
sort(unique(geneannotation1$external_gene_name))

```

# Generate a heatmap with AITL, PTCL-NOS with the extracted differentially expressed genes.

```{r}
gep<- geneExpr[,select_hist$sample.nameNEW]
mat<- gep[list_genes,]
rownames(mat) = unique(getBM( attributes = c("ensembl_transcript_id", "entrezgene", "external_gene_name"), filters = "entrezgene", values = gsub("_at","",rownames(mat)), mart = ensembl)$external_gene_name)
mycol= c("red","white","yellow")
mylabel = select_hist[,c("sample.nameNEW","final.molec","IDH2","RHOA","TET2","DNMT3A")]
rownames(mylabel) = mylabel$sample.nameNEW
mylabel$sample.nameNEW = NULL
mylabel.nocol = mylabel
mylabel.col = mylabel
mylabel.col[is.na(mylabel.col)]<-0
#head(mylabel.col)
mylabel.col$final.molec[mylabel.col$final.molec == "AITL"] = "black"; mylabel.col$final.molec[mylabel.col$final.molec == "PTCL.nos"] = "orange"
for (a in 2:5) mylabel.col[,a] = factor(mylabel.col[,a], levels = levels(as.factor(mylabel.col[,a])), labels = mycol )

mat  <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = mylabel.nocol, annotation_colors = list(final.molec = c(AITL = "black", PTCL.nos = "orange"),
                                  IDH2 = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]),
                                  RHOA = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]),
                                  TET2 = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]),
                                  DNMT3A = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]) ) ,
         border_color= NA, color = colorRampPalette(rev(brewer.pal(n = 5 , name = "RdYlBu")))(100), scale = "row", clustering_method = "ward.D2",clustering_distance_cols = "euclidean" )
```

# Use ConsensusClusterPlus to extract most significant clusters

```{r}
gep<- geneExpr[,select_hist$sample.nameNEW]
mat<- gep[list_genes,]
title=tempdir()
d<- data.matrix(mat)
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=6,
                               pFeature=1,
                               title=title,
                               clusterAlg="hc",
                               innerLinkage="ward.D2",
                               finalLinkage="ward.D2",
                               distance="euclidean",
                               seed=123456789)
icl = calcICL(results,title=title,plot="png")
#icl # compare clusters & scores and choose appropriately
kk<- as.data.frame((results[[4]]$consensusClass)) ##### 4 significant cluster
kk$geo.id<- rownames(kk)
colnames(kk)[1]<- "cluster"
table(kk$cluster)

```

# Plot heatmap according to the 4-cluster stratification.

```{r}
## add ConsensusCluster label, samples orderedaccordingly
heat<- merge(t(mat), kk, by.x = 0, by.y="geo.id")
heat2<- merge(heat, pts.info.data, by.x = 1, by.y="sample.nameNEW")
heat2<- heat2[order(heat2$cluster),]
mycol= c("red","white","yellow")
mylabel = heat2[,c("Row.names","cluster","final.molec","TET2","RHOA","IDH2","DNMT3A")]
colnames(mylabel)<- c("sample.names","clusters","Histology","TET2","RHOA","IDH2","DNMT3A")
rownames(mylabel) = mylabel$sample.names
mylabel$sample.names  = NULL
mylabel.nocol = mylabel
mylabel.col = mylabel
mylabel.col[is.na(mylabel.col)]<-0
#head(mylabel.col)
mylabel.col$Histology[mylabel.col$Histology == "AITL"] = "black"; mylabel.col$Histology[mylabel.col$Histology == "PTCL.nos"] = "orange"
#mylabel.col$Histology[mylabel.col$Histology == "ALCL.neg"] = "yellow"
for (a in c(3:6)) mylabel.col[,a] = factor(mylabel.col[,a], levels = levels(as.factor(mylabel.col[,a])), labels = mycol )
mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
for (a in 1) mylabel.col[,a] = factor(mylabel.col[,a], levels = levels(as.factor(mylabel.col[,a])), labels = mycol_plus[1:4] )
mylabel.nocol$clusters<-as.numeric(as.character(mylabel.nocol$clusters))
mylabel.nocol$clusters<-as.character(paste("cluster",mylabel.nocol$clusters, sep=""))

# plot heatmap
par(mfrow=c(1,1))
mat3<- t(data.matrix(heat2[,2:17]))
colnames(mat3)<-heat2$Row.names
mat3= mat3[order(rownames(mat3)),]
temp_name = getBM( attributes = c("ensembl_transcript_id", "entrezgene", "external_gene_name"), filters = "entrezgene", values = gsub("_at","",rownames(mat3)), mart = ensembl)[,c(2:3)]
temp_name = temp_name[!duplicated(temp_name[,1]),]
rownames(mat3) = temp_name [,2]
mat3  <- mat3 - rowMeans(mat3)
pheatmap(mat3, annotation_col = mylabel.nocol, annotation_colors = list(clusters = c(cluster1= mycol_plus[1], cluster2 = mycol_plus[2], cluster3 = mycol_plus[3], cluster4 = mycol_plus[4]), Histology = c(AITL = "black", PTCL.nos = "orange"), IDH2 = c(MUT=mycol[1], "NA"=mycol[2],WT=mycol[3]), RHOA = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]), TET2 = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]), DNMT3A = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]) ) ,  border_color= NA, color = colorRampPalette(rev(brewer.pal(n = 5 , name = "RdYlBu")))(100), scale = "row", cluster_cols = FALSE, show_colnames= F)

```

# Analyze sample stratification based on the extracted differentially expressed genes betwee AILT and PTCL-nos and the ALCL ALK-negative 3-gene model.

```{r}
select_hist<- pts.info.data[pts.info.data$final.molec == "AITL" | pts.info.data$final.molec == "PTCL.nos" | pts.info.data$final.molec == "ALCL.neg",]
# Add three classifier genes for ALCL ALK-neg [Agnelli et al, Blood, 2012]
# Check on array
anaplastic_gene<- c("TNFRSF8","BATF3","TMOD1")
geneannotation2 <- getBM( attributes = c("entrezgene", "external_gene_name"), filters = "external_gene_name", values = anaplastic_gene, mart = ensembl )

anaplastic_gene_ARRAY<- paste0(geneannotation2$entrezgene, "_at")

# Append 16-gene model to 3-gene model
list_genes_all<- c(list_genes, anaplastic_gene_ARRAY)

# Redo consensus cluster analysis
gep<- geneExpr[,select_hist$sample.nameNEW]
mat2<- gep[list_genes_all,]
title=tempdir()
d<- data.matrix(mat2)
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=6,
                               pFeature=1,
                               title=title,
                               clusterAlg="hc",
                               innerLinkage="ward.D2",
                               finalLinkage="ward.D2",
                               distance="euclidean",
                               seed=123456789)
icl = calcICL(results,title=title,plot="png")
#icl # compare clusters & scores and choose appropriately
kk2<- as.data.frame((results[[5]]$consensusClass)) ##### 4 significant cluster
kk2$geo.id<- rownames(kk2)
colnames(kk2)[1]<- "cluster"
table(kk2$cluster)
```

# Plot heatmap  AITL, PTCL-NOS, ALCL-neg and the 19-gene model

```{r fig.width = 15, fig.height = 7}

heat_5<- merge(t(mat2), kk2, by.x = 0, by.y="geo.id")
heat2<- merge(heat_5, pts.info.data, by.x = 1, by.y="sample.nameNEW")
heat2<- heat2[order(heat2$cluster),]
mycol= c("red","white","yellow")
mylabel = heat2[,c("Row.names","cluster","final.molec","TET2","RHOA","IDH2","DNMT3A")]
colnames(mylabel)<- c("sample.names","clusters","Histology","TET2","RHOA","IDH2","DNMT3A")
rownames(mylabel) = mylabel$sample.names
mylabel$sample.names  = NULL
mylabel.nocol = mylabel
mylabel.col = mylabel
mylabel.col[is.na(mylabel.col)]<-0
#head(mylabel.col)
mylabel.col$Histology[mylabel.col$Histology == "AITL"] = "black"; mylabel.col$Histology[mylabel.col$Histology == "PTCL.nos"] = "orange"; mylabel.col$Histology[mylabel.col$Histology == "ALCL.neg"] = "yellow"
for (a in c(3:6)) mylabel.col[,a] = factor(mylabel.col[,a], levels = levels(as.factor(mylabel.col[,a])), labels = mycol )
mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
for (a in 1) mylabel.col[,a] = factor(mylabel.col[,a], levels = levels(as.factor(mylabel.col[,a])), labels = mycol_plus[1:5] )
mylabel.nocol$clusters<-as.numeric(as.character(mylabel.nocol$clusters))
mylabel.nocol$clusters<-as.character(paste("cluster",mylabel.nocol$clusters, sep=""))


par(mfrow=c(1,1))
mat3<- t(data.matrix(heat2[,2:20]))
colnames(mat3)<-heat2$Row.names
mat3= mat3[order(rownames(mat3)),]
temp_name = getBM( attributes = c("ensembl_transcript_id", "entrezgene", "external_gene_name"), filters = "entrezgene", values = gsub("_at","",rownames(mat3)), mart = ensembl)[,c(2:3)]
temp_name = temp_name[!duplicated(temp_name[,1]),]
rownames(mat3) = temp_name [,2]
mat3  <- mat3 - rowMeans(mat3)
par(mfrow=c(1,1))
pheatmap(mat3, annotation_col = mylabel.nocol, annotation_colors = list(clusters = c(cluster1= mycol_plus[1], cluster2 = mycol_plus[2], cluster3 = mycol_plus[3], cluster4 = mycol_plus[4], cluster5 = mycol_plus[5]), Histology = c(AITL = "black", PTCL.nos = "orange", ALCL.neg= "yellow"), IDH2 = c(MUT=mycol[1], "NA"=mycol[2],WT=mycol[3]), RHOA = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]), TET2 = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]), DNMT3A = c(MUT=mycol[1],"NA"=mycol[2],WT=mycol[3]) ) ,  border_color= NA, color = colorRampPalette(rev(brewer.pal(n = 5 , name = "RdYlBu")))(100), scale = "row", cluster_cols = FALSE, show_colnames= F, row_annotation =3)

table(heat2$final.molec, heat2$cluster)

```

# LOOCV on AILT, ALCL neg, PTCLnos based on 19-gene model

```{r}
y = t(mat3)
cl.orig = c()
for (u in 1:nrow(y)) cl.orig [u] = unlist(strsplit(rownames(y)[u],"\\."))[1]

perm.mother = rownames(y)
perm.son = combn (perm.mother, length(perm.mother)-1)

output <- cbind(perm.mother, NA)

for (i in 1:length(perm.mother)) {
  train <- y [ perm.son[,i], ]
  test <- y [ ! ( rownames(y) %in% perm.son[,i]) , ]
  cl <- cl.orig [which(rownames(y)%in%perm.son[,i])]
  z <- lda(train, cl)
  p <- predict(z,test)$class
  output  [ setdiff(1:340, which( rownames(y) %in% perm.son[,i]) ) , 2  ] = as.character(p)
#  output  [ output[,1] == rownames(test) , 3  ] = z$scaling [1,1]
#  output  [ output[,1] == rownames(test) , 4  ] = z$scaling [2,1]
#  output  [ output[,1] == rownames(test) , 5  ] = z$scaling [3,1]
}

colnames(output) = c("true","LOOCV.predicted")
output = as.data.frame(output)
output$true.class = cl.orig

table(output$true.class, output$LOOCV.predicted  )

```

# Cibersort algorithm to characterize the tumour and microenviroment composition of each cluster

```{r fig.width = 15, fig.height = 12}
load(file.path(data.dir,"./Rmd.files//cibersort.all.Rdata"))
ciber_all<-as.data.frame.matrix(t(cibersort.percentages))
ciber_all$sample.nameNEW <- rownames(ciber_all)
colnames(kk2)[2]<-"sample.nameNEW"
require(plyr)
final <-join(ciber_all, kk2, by = "sample.nameNEW",  type="left")
final2<-merge(pts.info.data[,c(1,6,14:17)], final, by="sample.nameNEW")
final3<- subset(final2, final.molec %in% c("AITL","ALCL.neg","ALCL.pos","ATLL","NKT","PTCL.nos"))
final3<- final3[order(final3$final.molec),]
library(RColorBrewer)
n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

par(mar=c(7,5,10,8), xpd=TRUE)
x<- barplot(t(final3[7:28]), names.arg = rep("", length(final3$final.molec)), cex.names = 0.7, col=col_vector, border=NA)
legend("topright",legend=colnames(final3)[7:28], col=col_vector, pch=c(15), inset=c(-0.2,0), pt.cex= 0.7,
cex = 0.7, bty = "n",  x.intersp = 0.7)
names_hist<- unique(final3$final.molec)
col_hist<- c("orange","yellow","dodgerblue2","brown2","forestgreen","black")
num<- as.numeric(table(final3$final.molec))
for(i in (1:length(num)))
{
  segments(x[sum(num[1:i])+1-num[i]], 1.05,x[sum(num[1:i])],1.05,lwd=4, col=col_hist[i])
  text(x[(sum(num[1:i])-num[i] +1+ sum(num[1:i]))/2], 1.1, names_hist[i], cex=1, srt=0)
}


```

#Boxplot comparing the contribution of each cibersort signature between all extracted clusters

```{r}
par(mfrow=c(2,3))
par(mar=c(3,3,3,3), xpd=F)
for(i in (7:28))
{
  k<- as.numeric(final2[,i])
  table_wilk<- pairwise.wilcox.test(k,final2$cluster,p.adjust.methods = "BH" )$p.value
  df_wilk <- data.frame(expand.grid(dimnames(table_wilk)),array(table_wilk))

  df_wilk2_sig<- df_wilk[df_wilk$array.table_wilk.<0.05,]
  df_wilk2_sig$Var1<-as.numeric(as.character(df_wilk2_sig$Var1))
  df_wilk2_sig$Var2<-as.numeric(as.character(df_wilk2_sig$Var2))
  df_wilk2<-na.omit(df_wilk2_sig)
  if(nrow(df_wilk2)>0)
  {
  thresholds = boxplot(k~final2$cluster, plot = FALSE)$stats
  boxplot(k~final2$cluster, ylim=c(0,max(thresholds)+0.8), main=colnames(final2)[i], cex.main=1, col=mycol_plus, las=2)$stats
  for(j in (1:nrow(df_wilk2)))
  {
    segments(df_wilk2$Var1[j], max(k)+j/16, df_wilk2$Var2[j],max(k)+j/16, col = "grey40")
    p<-df_wilk2$array.table_wilk.[j]
    if(p<0.00001){p2 = "<0.00001"}else{
    p2<-as.numeric(formatC(p,digits=6,format="f"))}
    pval <- paste("p =",p2,sep=" ")
    text((df_wilk2$Var1[j]+ df_wilk2$Var2[j])/2,  max(k) +0.02 +j/16, cex = 0.8, pval, col= "grey40")
    }
  }
}

```
