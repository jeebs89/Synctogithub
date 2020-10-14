


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
            # Getting the loops with genes
                    loop_rna_overlaps_func <- function(input_a, input_b, row){
                        loops <- input_a %>% filter(fragmentMid1 < input_b$transcript_start[row] & fragmentMid2 > input_b$transcript_end[row])
                        if(nrow(loops)==0){NULL}
                        if(nrow(loops)>0){cbind(loops, input_b[row,])}
                    }
            # Calculating z-scores
                        zscore_sample <- function(df, row, column, mean_set, sd_set){
                            (df[row,column] - mean_set)/sd_set
                        }
                        zscore_row <- function(df, row, columns){
                            scores <- NULL
                            a <- mean(as.numeric(df[row,]))
                            b <- sd(df[row,])

                                for(i in 1:columns){
                                    scores <- cbind(scores, zscore_sample(df, row, i, a, b))
                                }
                                scores
                        }
            # Filtering the standard deviations
                    filtering_rna <- function(row, df){
                        data.frame("sd_h"=sd(df[row,6:8]), "sd_d"=sd(df[row,9:11]))

                    }

        #






        

    #



    hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
    tpm_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"
 
# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")

# Loop constraint, set to 0 if no constraint is wanted. This is in MB and will limit the size of a loop in the analysis 
loop_constraint <- 0

# Importing the looping data for MEF2a

    setwd(hic_data_directory)
    all_loops_a <- read.table("hic_contacts_mef2a_all_merged_q0.005_28sep2020.txt", header=T)
    all_loops_b <- all_loops_a %>% filter(cond=="disease" | cond=="healthy" | cond=="h_d", q.value < 0.005) %>% select(c(1:8, 10:12,16)) 

    if(loop_constraint==0){all_loops <- all_loops_b}
    if(loop_constraint>0){
        loop_size_limit <- loop_constraint*1000000     
        all_loops_c <- data.frame(all_loops_b, "loop_size"=(all_loops_b$fragmentMid2 - all_loops_b$fragmentMid1))
        all_loops <- all_loops_c %>% filter(loop_size < loop_size_limit) %>% select(c(1:12))

    }


# Importing the Genes from DESeq/RNA analysis

    setwd(tpm_wd)
    all_genes_a <- read.table("TPM_counts.csv", header=T, sep=",")
    colnames(all_genes_a)[1] <- "Gene"

    setwd(deseq_wd)
    deseq_result <- read.table("Top_unique_genes_DESeq_counts_05_Oct_2020.txt", header=T, sep=" ")
    deseq_genes <- left_join(deseq_result[,c(1:6, 12, 14:15)], all_genes_a[,c(1:4, 9:11)], by="Gene")


    # Filtering the TPM based on the top 85% of expression
        mean_TPM <- data.frame(deseq_genes, "h_tpm_mean"=rowMeans(deseq_genes[,10:12]), "d_tpm_mean"=rowMeans(deseq_genes[,13:15]))
        filtered_mean_tpm <- mean_TPM %>% filter(h_tpm_mean > quantile(mean_TPM$h_tpm_mean, 0.15) | d_tpm_mean > quantile(mean_TPM$d_tpm_mean, 0.15))


    zscores_a <- filtered_mean_tpm[,10:15]
    all_scores <- NULL
        for(i in 1:nrow(zscores_a)){
            all_scores <- rbind(all_scores, zscore_row(zscores_a, i, 6))
        }
        colnames(all_scores) <- colnames(zscores_a[1:6])



    all_genes <- data.frame("chr1"=paste0("chr", filtered_mean_tpm$chromosome_name), filtered_mean_tpm[,c(1:6, 8:9, 16:17)], all_scores)


    gene_rna_loops <- NULL
    for(j in 1:19){
        print(paste(Sys.time(), "..........matching the loops for chromosome: ", j,sep=""))

        all_loops_1 <- all_loops %>% filter(chr1==paste0("chr", j))
        all_genes_1 <- all_genes %>% filter(chr1==paste0("chr", j))

            chr_genes <- NULL
            for(i in 1:nrow(all_genes_1)){
                chr_genes <- rbind(chr_genes, loop_rna_overlaps_func(all_loops_1, all_genes_1, i))
            }
        gene_rna_loops <- rbind(gene_rna_loops, chr_genes)

    }

# Writing the whole table in case need to go back
write.table(gene_rna_loops, paste0("All_loop_rna_overlaps_master_", loop_constraint, "mb_limit_", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)

gene_rna_loops_b <- gene_rna_loops[c(13:19, 1:2,4, 20, 21, 5:12, 22:23,24:29)]
uniq_trans_loops <- gene_rna_loops_b %>% distinct(Gene, peak, .keep_all=TRUE)



# Filtering using the standard deviations to be less than 0.6
u_t_loops_full <- uniq_trans_loops %>% arrange(padj) %>% distinct(mgi_symbol, .keep_all=TRUE) %>% select(c(1,2,5,7,9:10,20:28))
u_t_loops_trim <- data.frame(u_t_loops_full[, c(2, 4, 7:15)])


sds <- NULL
for(i in 1:nrow(u_t_loops_trim)){
    sds <- rbind(sds, filtering_rna(i, u_t_loops_trim))
}




loop_sds <- data.frame(u_t_loops_trim, sds)
loops_sds_filtered <- loop_sds %>% filter(sd_h < 0.6 & sd_d < 0.6) %>% arrange(padj)
loops_final <- data.frame(loops_sds_filtered[,2:ncol(loops_sds_filtered)], row.names=loops_sds_filtered[,1])


write.table(loops_final, paste0("Top_loop_rna_overlaps_", loop_constraint, "mb_limit_", today_f, ".txt"), col.names=TRUE, row.names=TRUE, quote=FALSE)


# Meta data information
    meta_data <- data.frame("State"=c(rep("Healthy", 3), rep("HCM", 3)), row.names=colnames(loops_final[5:10]))
    Var1 <- c("dodgerblue", "firebrick3")
    names(Var1) <- c("Healthy", "HCM")
    anno_cols <- list(State = Var1)


# Meta data for the genes
    meta_data_genes <- data.frame("Loop"=loops_final$cond, row.names=loops_sds_filtered[,1])
    Var2 <- c("#d8b365", "#5ab4ac", "grey")
    names(Var2) <- c("healthy", "disease", "h_d")
    anno_cols_row <- list(Loop = Var2)
    anno_cols_test <- list(State = Var1, Loop=Var2)


colors <- colorRampPalette(c("#1F00EB", "black", "#FFFB0F"))



png(paste0("Heatmap_genes_looping_with_genes_top75_", loop_constraint, "mb_limit_", today_f, ".png"),
            width=500, height=1700, unit="px")


pheatmap(loops_final[,5:10], cluster_cols=T, 
        cluster_rows=T, color = colors(75),
        border_color=NA, clustering_method="complete", show_rownames=TRUE, show_colnames=FALSE,
        main="Healthy vs HCM Heatmap", #cellwidth = 35, #cellheight = 7.5,
        annotation_col=meta_data, annotation_row=meta_data_genes, annotation_colors=anno_cols_test
)

dev.off()






# TPM Export for counts data

# deseq_genes is the TPM
tpm_a <- data.frame("mgi_symbol"=row.names(loops_final), loops_final)
tpm_b <- left_join(tpm_a, deseq_genes[,c(1, 10:15)], by="mgi_symbol")
colnames(tpm_b) <- c(colnames(tpm_a), "tpm_h1", "tpm_h2", "tpm_h3", "tpm_d1", "tpm_d2", "tpm_d3")
write.table(tpm_b, paste0("Top_loop_rna_overlaps_TPM_COUNTS_", loop_constraint, "mb_limit_", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)







## GGplot for the adjusted pvalue threshold set in test3a

      # GGplot settings
                # Theme for general ggplot
                    theme<-theme_classic() +
                            theme(
                                text = element_text(color = "grey20"),
                                plot.title = element_text(hjust = 0.5),
                                axis.text.x = element_text(size=rel(1.5)), 
                                axis.text.y = element_text(size=rel(1.5), angle=90, hjust=0.5)
                            )

        #
    #

        # Plot Script
            setwd(deseq_wd)
            rowing <- c(1:nrow(test3))
            
            ggplot(test3, aes(x=rowing, y=padj)) +
                    geom_point(size=0.75) +
                    theme + 
                    geom_line(aes(y=0.02), col="grey20", linetype="dotted") +
                    labs(x="Gene",y="DESeq p.adjusted") 
            
            
            
            
            
            ggsave(paste0("padj_threshold_plot_", today_f, ".png"))

##






