

library('tidyverse', lib="/home/jcayford/r_libs")
library("DESeq2")
library('viridisLite')



    # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}
 

today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")




# Working Directories 
  rawcounts_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"
  export_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"

# Assigning condition (first 3 are controls, next are TAC)
  (condition <- factor(c(rep("ctl", 3), rep("tac", 3))))

# Data Import and DESeq

setwd(rawcounts_wd)

  # Opening up the raw counts file and setting as matrix
      counts <- read.table("raw_counts.csv", sep=",", header=T, row.names=1)
      counts_2 <- counts[, c(1:3, 8:10)]
      countdata <- as.matrix(counts_2)

  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
      (coldata <- data.frame(row.names=colnames(countdata), condition))
      dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
      dds
  # Run the DESeq pipeline
      dds <- DESeq(dds)


setwd(export_Wd)
# Plot dispersions
  png("qc-dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
# Regularized log transformation for clustering/heatmaps, etc
  rld <- rlogTransformation(dds)
  head(assay(rld))
  hist(assay(rld))

# Some QC Ploits

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
#


# Generation of the gene list with a padj < 0.05. Set up for mouse data to get the gene names
    # Get differential expression results
      res <- results(dds)
      table(res$padj<0.05)
      ## Order by adjusted p-value
      res <- res[order(res$padj), ]
    # Merge with normalized count data
      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene" 
      resdata_2 <- resdata %>% filter(padj < 0.05)

    # Getting gene names:
      library("biomaRt")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      mouse_mgi <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = mouse)

    # Combining the lists to get the gene name from the transcript and exporting the result for the transcripts and genes.
      genes_list <- data.frame("ensembl_gene_id"=as.factor(substr(resdata_2[,1], 1, 18)), resdata_2[,c(1:(ncol(resdata_2)-1))])
      all_transcripts <- left_join(genes_list, mouse_mgi, by="ensembl_gene_id")
      all_trans_ordered <- all_transcripts[c(14, 1:2, 4, 7:13, 15:18)]
      unique_genes <- all_trans_ordered %>% distinct(mgi_symbol, .keep_all=TRUE)

    # Exporting the files
      write.table(all_trans_ordered, paste0("top_transcripts_DESeq_counts_", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
      write.table(unique_genes, paste0("Top_unique_genes_DESeq_counts_", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)






## MA plot
## Could do with built-in DESeq2 function: DESeq2::plotMA(dds, ylim=c(-1,1), cex=1) or can use a custom one:
    maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
      with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
      with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=0.5))
      if (labelsig) {
        require(calibrate)
        with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
      }
    }

png(paste0("MA_Plot_healthy_hcm_p005_",today_f, ".png"), 1500, 1000, pointsize=20)
maplot(resdata, main="Heathly vs HCM MA Plot - p < 0.05")
dev.off()


# Volcano Plots: 
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
        png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
        volcanoplot(resdata, lfcthresh=1, sigthresh=0.005, textcx=.8, cex=0.6)
        dev.off()

## Volcano plot with a threshold. USE THIS ONE
resdata_vol <- resdata
fc_filter <- 2.5

resdata_vol$log2FoldChange[(resdata_vol$log2FoldChange > fc_filter)] <- fc_filter
resdata_vol$log2FoldChange[(resdata_vol$log2FoldChange < -fc_filter)] <- (-fc_filter)

png(paste0("healthy_HCM_volcanoplot_2.5fcthresh_", today_f, ".png"), 1200, 1000, pointsize=20)
volcanoplot(resdata_vol, lfcthresh=1, sigthresh=0.005, textcx=.8, cex=0.6, labelsig=TRUE)
dev.off()


















