

library('tidyverse', lib="/home/jcayford/r_libs")
library("DESeq2")
library('viridisLite')



    # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}
 





rawcounts_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"





setwd(rawcounts_wd)

# Opening up the raw counts file and setting as matrix
    counts <- read.table("raw_counts.csv", sep=",", header=T, row.names=1)
    counts_2 <- counts[, c(1:3, 8:10)]
    countdata <- as.matrix(counts_2)

# Assigning condition (first 3 are controls, next are TAC)
    (condition <- factor(c(rep("ctl", 3), rep("tac", 3))))


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
    (coldata <- data.frame(row.names=colnames(countdata), condition))
    dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
    dds





# Export wd: 
export_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"



# Run the DESeq pipeline
    dds <- DESeq(dds)



# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))



# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="healthy_hcm_deseq_result.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")




## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("healthy_hcm_p005.png", 1500, 1000, pointsize=20)
maplot(resdata, main="Heathly vs HCM MA Plot - p:0.05")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.005, main="Healthy vs HCM DESeq Genes", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=16, col="gray30", main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=16, col="firebrick", ...))
  #with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=16, col="orange", ...))
  #with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=16, col="dodgerblue", ...))
  #if (labelsig) {
  #  require(calibrate)
  #with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  #}
  #legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=16, col=c("firebrick","orange","dodgerblue"))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep="")), pch=16, col=c("firebrick"))
}



#png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
#volcanoplot(resdata, lfcthresh=1, sigthresh=0.005, textcx=.8, cex=0.6)
#dev.off()


resdata_vol <- resdata

fc_filter <- 2.5

resdata_vol$log2FoldChange[(resdata_vol$log2FoldChange > fc_filter)] <- fc_filter
resdata_vol$log2FoldChange[(resdata_vol$log2FoldChange < -fc_filter)] <- (-fc_filter)

png("healthy_HCM_volcanoplot_2.5fcthresh.png", 1200, 1000, pointsize=20)
volcanoplot(resdata_vol, lfcthresh=1, sigthresh=0.005, textcx=.8, cex=0.6, labelsig=TRUE)
dev.off()


















