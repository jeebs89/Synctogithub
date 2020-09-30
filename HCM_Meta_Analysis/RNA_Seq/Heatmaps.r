# Designed to make the heatmaps in R



library('tidyverse', lib="/home/jcayford/r_libs")
library('viridisLite')
library('readr')
library('pheatmap')
library('gplots')



    # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}         
            # Calculating z-scores
                zscore_sample <- function(row, column, mean_set, sd_set){
                    (gene_cluster[row,column] - mean_set)/sd_set
                }
                zscore_row <- function(row, columns){
                    scores <- NULL
                    a <- mean(as.numeric(gene_cluster[row,]))
                    b <- sd(gene_cluster[row,])

                        for(i in 1:columns){
                            scores <- cbind(scores, zscore_sample(row, i, a, b))
                        }
                        scores
                }

    #




tpm_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"
deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"


# Uploading the TPM data from the RNA-Seq pipeline
setwd(tpm_wd)
tpms_a <- read.table("TPM_counts.csv", header=T, sep=",")
tpms <- tpms_a[,c(1:4, 9:11)]
colnames(tpms)[1] <- "Gene"

# Uploading the DESeq results
setwd(deseq_wd)
deseq_a <- read.table("healthy_hcm_deseq_result.csv", header=T, sep=",")
deseq <- deseq_a[,c(2:ncol(deseq_a))]
deseq_no_counts <- deseq[,c(1:7)]


# Joining the data
joined <- left_join(deseq_no_counts, tpms, by="Gene")


# pheatmap for the top selected genes
number_of_genes <- "all" # this can be a number or "all"

sig_joined <- joined %>% filter(padj < 0.05)
sig_2 <- data.frame(sig_joined[,c(1:7)], (log2(sig_joined[,c(8:13)]+1)))


if(number_of_genes=="all"){number_of_genes <- nrow(sig_2)}
res <- pheatmap(sig_2[c(1:number_of_genes), c(8:13)], cluster_cols=F)

gene_cluster_a <- cbind(sig_2[c(1:number_of_genes), c(1, 8:13)],cluster = cutree(res$tree_row, k = 10))
gene_cluster <- gene_cluster_a %>% arrange(cluster) %>% select(2:ncol(gene_cluster_a))
gene_cluster_genes <- gene_cluster_a %>% arrange(cluster) %>% select(1)


all_scores <- NULL
    for(i in 1:nrow(gene_cluster)){
        all_scores <- rbind(all_scores, zscore_row(i, 6))
    }
    colnames(all_scores) <- colnames(gene_cluster[1:6])




meta_data <- data.frame("State"=c(rep("Healthy", 3), rep("HCM", 3)))

row.names(meta_data) <- colnames(all_scores)

# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- c("dodgerblue", "firebrick3")
names(Var1) <- c("Healthy", "HCM")
anno_cols <- list(State = Var1)

#png(paste0("top_", number_of_genes, "_genes_heatmap.png"))
pheatmap(all_scores, cluster_cols=T, 
        cluster_rows=T, color = viridis(256, option="magma", direction = -1),
        border_color=NA, clustering_method="complete", show_rownames=FALSE, show_colnames=FALSE,
        main="Healthy vs HCM Heatmap", annotation_col=meta_data, annotation_colors=anno_cols
)
#dev.off()



### YELLOW BLUE HEATMAP

setwd("/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq/")


colors <- colorRampPalette(c("#1F00EB", "black", "#FFFB0F"))




a <- data.frame(all_scores, "hmean"=rowMeans(all_scores[,1:3]), "dmean"=rowMeans(all_scores[,4:6]))

b <- a %>% filter(hmean < (-0.35) & dmean > 0.35)
c <- a %>% filter(dmean < (-0.35) & hmean > 0.35)
d <- rbind(b,c)


png(paste0("all_DESeq2_genes_heatmap.png"),
                width=400, height=1250, unit="px",  pointsize = 12)

pheatmap(d[,1:6], cluster_cols=T, 
        cluster_rows=T, color = colors(75),
        border_color=NA, clustering_method="complete", show_rownames=FALSE, show_colnames=FALSE,
        main="Healthy vs HCM All DEseq Genes", cellwidth = 35, #cellheight = 3.5,
        annotation_col=meta_data, annotation_colors=anno_cols
)

dev.off()



all_scores_flipped <- all_scores[nrow(all_scores):1,]






# Getting gene names:
library("biomaRt")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_mgi <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = mouse)


genes_list <- data.frame("ensembl_gene_id"=as.factor(substr(gene_cluster_genes[,1], 1, 18)), gene_cluster[,c(1:(ncol(gene_cluster)-1))])
all_genes <- left_join(genes_list, mouse_mgi, by="ensembl_gene_id")
unique_genes <- all_genes %>% distinct(mgi_symbol, .keep_all=TRUE)
short_unique_genes <- unique_genes[, c(1,8)]

write.table(all_genes, paste0("top_", number_of_genes, "_genes_withTPM_28SEP2020.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(short_unique_genes, paste0("top_", number_of_genes, "_unique_genes_with_TPM_28SEP2020.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(all_scores, paste0("top_", number_of_genes, "_genes_zscores_28SEP2020.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)


