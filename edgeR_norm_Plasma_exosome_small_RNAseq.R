
library(edgeR)
library(locfit)
library(statmod)
library(ggfortify)
library(factoextra)
library(gplots)
library(RColorBrewer)
library(calibrate)

## Read the raw read counts
all_samples<-read.table("UMI_reads.txt", h=T,check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
head(all_samples)
condition <- c(rep("control",26), rep("case",26))
y<-DGEList(counts=all_samples[,2:53],genes=all_samples$miRNA, group=condition)
dim(y)
#[1] 2556   52
rownames(y$counts) <- rownames(y$genes) <- y$genes$genes

## Before QC Plots
#Set up treatment colours
colorCodes <- c(rep("blue",26), rep("red",26))
#Read counts per sample and counts per gene
pdf("miRNA_UMI_before.pdf", h=4, w=6)
barplot(y$samples$lib.size*1e-6, col=colorCodes,las=2, names=colnames(y), ylab="Library size (millions)")
legend("topleft", legend=c("T1D", "Control"), bty="n", col=c(2,4), pch=15, cex=0.8)
dev.off()

tiff("miRNA_UMI_before.tiff", width=12, height=6,units = "in", compression="lzw",res=300)
barplot(y$samples$lib.size*1e-6, col=colorCodes,las=2, names=colnames(y), ylab="Library size (millions)")
legend("topright", legend=c("T1D", "Control"), lty=1, bty="n", cex=0.8, lwd=4, col=c("red","blue"))
dev.off()

#The function plotMDS produces a plot in which distances between samples correspond to leading biological coefficient of variation (BCV) between those samples:
pdf("MDS_plot_before.pdf")
plotMDS(y,col=colorCodes, main="MDS Plot", cex=0.7)
legend("topleft", legend=c("T1D", "Control"), bty="n", cex=0.6,col=c(2,4), pch=15)
dev.off()

#In the plot, dimension 1 separates the treated from the normal samples, while dimension 2 roughly corresponds to patient number. This confirms the paired nature of the samples.


## Filter lowly expressed genes
selr <- rowSums(y$counts>10) >=20
summary(selr)
#   Mode   FALSE    TRUE 
#logical     267    2289 

#Recompute library sizes
y <- y[selr, , keep.lib.sizes=FALSE]

dim(y)
#[1] 2289   52

#TMM normalization
y <- calcNormFactors(y)
y$samples
#          group lib.size norm.factors
#C01_S1  control   673681   1.46771046
#C02_S2  control   505444   2.14217436
#C03_S3  control  2807119   0.38813040
#C04_S4  control   501996   1.54973658
#C05_S5  control   729594   1.30334705
#C06_S6  control   919238   1.36678594
#C07_S7  control  1230345   0.90742885
#C08_S8  control   419020   2.55960018
#C09_S9  control  1203216   0.98575722
#C10_S10 control   861173   1.15838006
#C11_S11 control  1003747   1.10488803
#C12_S12 control   509739   1.17205738
#C13_S13 control   471197   1.23551221
#C14_S14 control   728159   1.31165612
#C15_S15 control   556235   1.05656363
#C16_S16 control  1770207   0.55696559
#C17_S17 control   943382   1.14046956
#C18_S18 control   721926   1.19157782
#C19_S19 control   721427   1.07871663
#C20_S20 control   328415   1.19455560
#C21_S21 control   212901   1.47900310
#C22_S22 control   615161   0.88333000
#C23_S23 control  1200198   0.85050347
#C24_S24 control   699388   0.62348936
#C25_S25 control   340767   1.35570339
#C26_S26 control   464582   0.95133528
#G01_S1     case  7934276   0.26600779
#G02_S2     case  8130160   0.25589626
#G03_S3     case 16777284   0.05539319
#G04_S4     case  3211956   0.29489692
#G05_S5     case   461221   1.59491618
#G06_S6     case  1581813   0.58406006
#G07_S7     case  1263264   1.30960877
#G08_S8     case  2248209   1.37218621
#G09_S9     case  5859485   0.21856210
#G10_S10    case   951885   1.08681901
#G11_S11    case   732280   1.09101875
#G12_S12    case   541554   1.26052505
#G13_S13    case  2018830   0.78593244
#G14_S14    case   766227   1.86030647
#G15_S15    case   984015   1.37637708
#G16_S16    case   818312   2.14743607
#G17_S17    case   562615   1.85048276
#G18_S18    case  1132111   0.90543128
#G19_S19    case   565470   2.39650472
#G20_S20    case   661597   1.58344337
#G21_S21    case   514184   1.52342865
#G22_S22    case   375023   1.24223213
#G23_S23    case   740459   1.26336484
#G24_S24    case  1067236   0.75472352
#G25_S25    case   330092   1.33499231
#G26_S26    case   317689   1.43916182


logUMI=log2(y$counts+0.2)
logcpm <- cpm(y, prior.count=2, log=TRUE)
write.table(logcpm, file="logCPM_filtered_miRNAs.txt", sep="\t", quote=FALSE)
write.table(logUMI, file="logUMI_filtered_miRNAs.txt", sep="\t", quote=FALSE)

## Plots After filtering low expressed miRNAs
#Set up treatment colours
colorCodes <- c(rep("blue",26), rep("red",26))

#Read counts per sample and counts per gene
pdf("miRNA_UMI_after.pdf", h=4, w=6)
barplot(y$samples$lib.size*1e-6, col=colorCodes,las=2, names=colnames(y), ylab="Library size (millions)")
legend("topleft", legend=c("T1D", "Control"), bty="n", col=c(2,4), pch=15, cex=0.8)
dev.off()

tiff("miRNA_UMI_after.tiff", width=12, height=6,units = "in", compression="lzw",res=300)
barplot(y$samples$lib.size*1e-6, col=colorCodes,las=2, names=colnames(y), ylab="Library size (millions)")
legend("topright", legend=c("T1D", "Control"), lty=1, bty="n", cex=0.8, lwd=4, col=c("red","blue"))
dev.off()

#The function plotMDS produces a plot in which distances between samples
#correspond to leading biological coefficient of variation (BCV) between those samples:
pdf("MDS_plot_after")
plotMDS(y,col=colorCodes, main="MDS Plot", cex=0.7)
legend("topleft", legend=c("T1D", "Control"), bty="n", cex=0.6,col=c(2,4), pch=15)
dev.off()
#In the plot, dimension 1 separates the treated from the normal samples, while dimension 2 roughly corresponds to patient number. This confirms the paired nature of the samples.

#PCA plot
df <- t(y$counts)
pdf("PCA_plot.pdf")
autoplot(prcomp(df), colour = colorCodes, label = TRUE,shape = FALSE)

#2nd PCA plot using factoextra package
cond=as.factor(condition)
res.pca <- prcomp(df, scale = TRUE)
fviz_pca_ind(res.pca, habillage = cond, addEllipses = TRUE, ellipse.level = 0.68) + theme_minimal()
dev.off()


## DE analysis using GLM aproach
design=model.matrix(~0+group, data= y$samples)  ## group is from y$samples$group
# +0 in the model is an instruction not to include an intercept column and instead to include a column for each group

colnames(design)<-levels(y$samples$group)
design
#        case control
#C01_S1     0       1
#C02_S2     0       1
#C03_S3     0       1
#C04_S4     0       1
#C05_S5     0       1
#C06_S6     0       1
#C07_S7     0       1
#C08_S8     0       1
#C09_S9     0       1
#C10_S10    0       1
#C11_S11    0       1
#C12_S12    0       1
#C13_S13    0       1
#C14_S14    0       1
#C15_S15    0       1
#C16_S16    0       1
#C17_S17    0       1
#C18_S18    0       1
#C19_S19    0       1
#C20_S20    0       1
#C21_S21    0       1
#C22_S22    0       1
#C23_S23    0       1
#C24_S24    0       1
#C25_S25    0       1
#C26_S26    0       1
#G01_S1     1       0
#G02_S2     1       0
#G03_S3     1       0
#G04_S4     1       0
#G05_S5     1       0
#G06_S6     1       0
#G07_S7     1       0
#G08_S8     1       0
#G09_S9     1       0
#G10_S10    1       0
#G11_S11    1       0
#G12_S12    1       0
#G13_S13    1       0
#G14_S14    1       0
#G15_S15    1       0
#G16_S16    1       0
#G17_S17    1       0
#G18_S18    1       0
#G19_S19    1       0
#G20_S20    1       0
#G21_S21    1       0
#G22_S22    1       0
#G23_S23    1       0
#G24_S24    1       0
#G25_S25    1       0
#G26_S26    1       0

#attr(,"assign")
#[1] 1 1
#attr(,"contrasts")
#attr(,"contrasts")$group
#[1] "contr.treatment"


#Estimation of dispersion
disp <- estimateDisp(y, design, robust=TRUE)
disp$common.dispersion
#[1] 0.3713107

#plot dispersion values
pdf("BCV_plot_new.pdf")
plotBCV(disp, main="BCV Plot")
dev.off()

#Fit genewise glms
fit <- glmFit(disp, design)
lrt<- glmLRT(fit, contrast=c(1,-1))  ### -1 is for control +1 is for case 
topTags(lrt)

res = as.data.frame(topTags(lrt, n = Inf))

#Extract the significantly differentially expressed genes
resOrdered = res[order(res$PValue),]
resSig = subset(resOrdered, FDR<0.05)
#Print results to file
write.table(resOrdered, file='casevscontrol_DEResults.txt',sep='\t',quote=FALSE)
write.table(resSig, file='casevscontrol_DE_pVal0.05.txt',sep='\t',quote=FALSE)


## Use a cutoff of log2-fold-change of 1 for DE plots
is.de <- decideTestsDGE(lrt, p.value=0.05, lfc=1)
summary(is.de)

#summary(is.de)
#       1*case -1*control
#Down                   9
#NotSig              2113
#Up                   167

#Select genes with FDR < 0.05 and logFC <= 1 to highlight on plot

pdf("DE_genes_logFC-1.pdf")
plotSmear(lrt, de.tags=rownames(lrt) [is.de!=0], xlab="Average logCPM", ylab="logFC", cex=0.6)
abline(h = c(-1, 0, 1), col = c("dodgerblue", "red", "dodgerblue"), lty=2)
dev.off()

pdf("DE_genes_logFC-1_scale.pdf")
plotSmear(lrt, de.tags=rownames(lrt) [is.de!=0], xlab="Average logCPM", ylab="logFC", cex=0.6,ylim=c(-5,5))
abline(h = c(-1, 0, 1), col = c("dodgerblue", "red", "dodgerblue"), lty=2)
dev.off()


## Heatmaps and Volcano Plots
table1=read.table("logUMI_filtered_miRNAs.txt", h=T)
hmcol<- colorRampPalette(brewer.pal(10, "RdYlGn"))(256)

#Other optional palettes
#hmcol<- colorRampPalette(brewer.pal(10, "RdBu"))(256)
#col=redgreen(75)

table1_matrix=as.matrix(table1)
pdf("allxall_correlation.matrix.pdf",width=12, height=12)
heatmap.2(cor(table1_matrix), key=TRUE,symkey=FALSE, dendrogram = "none",na.rm = TRUE, scale="none", margins=c(10,10),symbreaks=FALSE, density.info="density", trace="none", main="Sample Correlations",cexRow=0.8, cexCol=0.8, col=hmcol, , keysize=1, key.title="Color key", key.xlab=NULL, key.ylab=NULL, key.xtickfun=NULL, key.ytickfun=NULL)
dev.off()

# Volcano plot
table2=read.table("casevscontrol_DEResults.txt", h=T)

#Choose column 6 for adj p value
#cutoff for Fc lines shd be 1, -1 logFC

logFC<-table2[,2]
pvals<-table2[,6]
names_gene<-table2[,1]

range(logFC)
#--2.842820  3.129462

summary(-log10(pvals))
#[1]   0.00000 10.7

#Change range inside volcano_plot function accordingly

volcano_plot<-function(logFC, pval) {
     x_range <- range((logFC), -5, 5)
     y_range <- range(c(-log10(pvals), 0, 11))

     plot(x_range,                                 # x-dim
          y_range,                                   # y-dim
          type="n",                                  # empty plot
          xlab="log2FC",                   # x-axis title
          ylab="-log10 adj.P-value",              # y-axis title
          main="",                       # plot title
          )
     
     ##For labels
  #textxy((logFC),(-log10(pvals)), labs=names, cex = 0.6, m = c(0, 0))

     abline(h=-log10(0.05),col="green",lty="44")# horizontal line at P=0.05
     abline(v=c(-1,1),col="violet",lty="1343")  # vertical lines at 1.5-fold
     ## Define colors based on their values:
     ## Not significant: purple
     ## Significant and smaller than half fold change: blue
     ## Significant and larger than two fold change: red
     ## Significant but between half and two fold change: orange
     color <- ifelse(-log10(pvals)>(-log10(0.05)),
                     ifelse((logFC)>(-1),
                            ifelse((logFC<1),
                                   "orange",
                                   "red"),
                            "blue"),
                     "black")
     points(
            (logFC),
            -log10(pvals),
            col=color,
            pch=20
            )

}

volcano_plot(logFC, pvals)

pdf("volcano_plot.pdf")
volcano_plot(logFC, pvals)
dev.off()

