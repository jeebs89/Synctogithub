
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

    #



# Importing the looping data for MEF2a
    hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    setwd(hic_data_directory)
    all_loops_a <- read.table("hic_contacts_mef2a_all_merged_q0.005_28sep2020.txt", header=T)
    all_loops <- all_loops_a %>% filter(cond=="disease" | cond=="healthy" | cond=="h_d") %>% select(c(1:8, 10:12,16))

# Importing the Genes from DESeq/RNA analysis
    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
    setwd(deseq_wd)
    all_genes_a <- read.table("top_1886_genes_withTPM_28SEP2020.txt", header=T)
    all_genes <- data.frame("chr1"=paste0("chr", all_genes_a$chromosome_name), all_genes_a[,c(1:8, 11:12)])

# Generating the gene ID and gene name with the zscores
    all_genes_zscore_a <- read.table("top_1886_genes_zscores_28SEP2020.txt", header=T)
    all_genes_unique <- all_genes %>% distinct(mgi_symbol, .keep_all=TRUE)
    all_genes_zscore <- data.frame(all_genes_unique[,c(2,9)], all_genes_zscore_a)




loop_rna_overlaps_func <- function(input_a, input_b, row){
    loops <- input_a %>% filter(fragmentMid1 < input_b$transcript_start[row] & fragmentMid2 > input_b$transcript_end[row])
    if(nrow(loops)==0){NULL}
    if(nrow(loops)>0){cbind(loops, input_b[row,c(2,9,3:8, 10:11)])}
}


gene_rna_loops <- NULL
for(j in 1:19){
# selecting chr1
    all_loops_1 <- all_loops %>% filter(chr1==paste0("chr", j))
    all_genes_1 <- all_genes %>% filter(chr1==paste0("chr", j))

        chr_genes <- NULL
        for(i in 1:nrow(all_genes_1)){
            chr_genes <- rbind(chr_genes, loop_rna_overlaps_func(all_loops_1, all_genes_1, i))
        }
    gene_rna_loops <- rbind(gene_rna_loops, chr_genes)

}

unique_gene_rna_loops <- gene_rna_loops %>% distinct(mgi_symbol, .keep_all=TRUE) %>% select(c(14:20,12))
### NEED TO GET THE DESeq PVALUES!!!!!


### 29SEP2020:

# Adding the h or d loops for the meta data
    mean_TPM <- data.frame("mgi_symbol"=unique_gene_rna_loops[,1], "cond"=unique_gene_rna_loops[,8],  "hmean"=rowMeans(unique_gene_rna_loops[,2:4]), "dmean"=rowMeans(unique_gene_rna_loops[,5:7]))
    filtered_mean_tpm <- mean_TPM %>% filter(hmean > 3.5 | dmean > 3.5)

    test <- data.frame(filtered_mean_tpm, "sub"=(filtered_mean_tpm$hmean - filtered_mean_tpm$dmean))
    test2 <- test %>% arrange(sub)
    test3 <- left_join(test2, all_genes_zscore, by="mgi_symbol")


# Meta data information
    meta_data <- data.frame("State"=c(rep("Healthy", 3), rep("HCM", 3)))
    Var1 <- c("dodgerblue", "firebrick3")
    names(Var1) <- c("Healthy", "HCM")
    anno_cols <- list(State = Var1)


# Meta data for the genes
    meta_data_genes <- data.frame("Loop"=test3$cond, row.names=test3[,1])
    Var2 <- c("#d8b365", "#5ab4ac", "grey")
    names(Var2) <- c("healthy", "disease", "h_d")
    anno_cols_row <- list(Loop = Var2)
    anno_cols_test <- list(State = Var1, Loop=Var2)


test4 <- data.frame(test3[,c(7:12)], row.names=test3[,1])



number_of_genes <- 50

test5 <- test4[c(1:number_of_genes, (nrow(test4)-number_of_genes):nrow(test4)),]


setwd("/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq/")

colors <- colorRampPalette(c("#1F00EB", "black", "#FFFB0F"))


png(paste0("top_", number_of_genes, "_updown_genes_looping_heatmap.png"),
                width=400, height=1000, unit="px")
pheatmap(test5, cluster_cols=T, 
        cluster_rows=T, color = colors(75),
        border_color=NA, clustering_method="complete", show_rownames=FALSE, show_colnames=FALSE,
        main="Healthy vs HCM Heatmap", cellwidth = 35, cellheight = 7.5,
        annotation_col=meta_data, annotation_row=meta_data_genes, annotation_colors=anno_cols_test
)
dev.off()
