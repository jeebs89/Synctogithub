


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



    
    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
    tpm_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"
    hic_data_directory <- "/home/jcayford/HCM/Analysis/Mustache/hic_files"
    chip_beds_directory <- "/home/jcayford/HCM/Analysis/ChIP_files/"
            
        
        #chip_bed <- "polr2a.bed" 
        chip_bed <- "mef2a.bed"
        tf <- "mef2"  
        qval <- 0.01

 
# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")




# Importing the looping data for MEF2a


            # FitHiC data
                 # ChIP Peaks
                setwd(chip_beds_directory)
                chip_peaks_a <- read.table(chip_bed)
                        if(chip_peaks_a[1,1]==1){chip_peaks <- data.frame(paste("chr", chip_peaks_a[,1], sep=""), chip_peaks_a[,2:5])}
                        if(chip_peaks_a[1,1]=="chr1"){chip_peaks <- chip_peaks_a[,1:5]}                       
                colnames(chip_peaks) <- c("chr1", "start", "end", "peak", "height")     



            # FitHiC data
                setwd(hic_data_directory)
                hic_data_a <- read.table(paste("mustache_hic_filtered_qval", qval, ".txt", sep=""), header=TRUE)
                hic_data <- data.frame("chr1"=hic_data_a[,1], "fragmentMid1"=rowMeans(hic_data_a[,2:3]), "fragmentMid2"=rowMeans(hic_data_a[,5:6]), 
                                    "FDR"=hic_data_a$FDR, "condition"=hic_data_a$condition)


        # Filtering FitHiC data with the ChIP-Seq Peaks

            chip_means <- round_to(rowMeans(chip_peaks[,2:3]), to=2500)
           
            # Filtering the HiC data
                hic_filtered_a <- NULL
                for(i in c((-10000), (-7500), (-5000), (-2500), 0, 2500, 5000, 7500, 10000)){
                    hic_filtered_a <- rbind(hic_filtered_a, chip_filter(i))
                }                
                hic_filtered <- hic_filtered_a %>% distinct(fragmentMid1, fragmentMid2,  .keep_all=TRUE)


                all_loops_m <- hic_filtered %>% filter(FDR < qval)


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


    gene_rna_loops_m <- NULL
    for(j in 1:19){
        print(paste(Sys.time(), "..........matching the loops for chromosome: ", j,sep=""))

        all_loops_1 <- all_loops_m %>% filter(chr1==paste0("chr", j))
        all_genes_1 <- all_genes %>% filter(chr1==paste0("chr", j))

            chr_genes <- NULL
            for(i in 1:nrow(all_genes_1)){
                chr_genes <- rbind(chr_genes, loop_rna_overlaps_func(all_loops_1, all_genes_1, i))
            }
        gene_rna_loops_m <- rbind(gene_rna_loops_m, chr_genes)

    }

# Writing the whole table in case need to go back
write.table(gene_rna_loops_m, paste0("All_loop_rna_overlaps_master_", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)

gene_rna_loops_b_m <- gene_rna_loops_m[c(1, 10:15, 2:3, 16, 17, 4, 6:8, 18:25)]
uniq_trans_loops_m <- gene_rna_loops_b_m %>% distinct(Gene, peak, .keep_all=TRUE)



# Filtering using the standard deviations to be less than 0.6
u_t_loops_full_m <- uniq_trans_loops_m %>% arrange(padj) %>% distinct(mgi_symbol, .keep_all=TRUE) %>% select(c(1, 2, 5, 7:9, 16:23))
u_t_loops_trim_m <- data.frame(u_t_loops_full[, c(2, 4)], "cond"=NA, u_t_loops_full_m[,c(7:14)])




sds_m <- NULL
for(i in 1:nrow(u_t_loops_trim_m)){
    sds_m <- rbind(sds_m, filtering_rna(i, u_t_loops_trim_m))
}




loop_sds_m <- data.frame(u_t_loops_trim_m, sds_m)
loops_sds_filtered_m <- loop_sds_m %>% filter(sd_h < 0.6 & sd_d < 0.6) %>% arrange(padj)
loops_final_m <- data.frame(loops_sds_filtered_m[,2:ncol(loops_sds_filtered_m)], row.names=loops_sds_filtered_m[,1])


write.table(loops_final_m, paste0("Top_loop_rna_overlaps_mustache_", today_f, ".txt"), col.names=TRUE, row.names=TRUE, quote=FALSE)


# Meta data information
    meta_data <- data.frame("State"=c(rep("Healthy", 3), rep("HCM", 3)), row.names=colnames(loops_final_m[5:10]))
    Var1 <- c("dodgerblue", "firebrick3")
    names(Var1) <- c("Healthy", "HCM")
    anno_cols <- list(State = Var1)


# Meta data for the genes
    meta_data_genes <- data.frame("Loop"=loops_final_m$cond, row.names=loops_sds_filtered_m[,1])
    Var2 <- c("#d8b365", "#5ab4ac", "grey")
    names(Var2) <- c("healthy", "disease", "h_d")
    anno_cols_row <- list(Loop = Var2)
    anno_cols_test <- list(State = Var1, Loop=Var2)


colors <- colorRampPalette(c("#1F00EB", "black", "#FFFB0F"))

png(paste0("Heatmap_genes_looping_with_genes_top75_", today_f, ".png"),
            width=500, height=1700, unit="px")


pheatmap(loops_final_m[,5:10], cluster_cols=T, 
        cluster_rows=T, color = colors(75),
        border_color=NA, clustering_method="complete", show_rownames=TRUE, show_colnames=FALSE,
        main="Healthy vs HCM Heatmap", #cellwidth = 35, #cellheight = 7.5,
        annotation_col=meta_data, annotation_colors=anno_cols_test
)


dev.off()






# TPM Export for counts data

# deseq_genes is the TPM
tpm_a <- data.frame("mgi_symbol"=row.names(loops_final), loops_final)
tpm_b <- left_join(tpm_a, deseq_genes[,c(1, 10:15)], by="mgi_symbol")
colnames(tpm_b) <- c(colnames(tpm_a), "tpm_h1", "tpm_h2", "tpm_h3", "tpm_d1", "tpm_d2", "tpm_d3")
write.table(tpm_b, paste0("Top_loop_rna_overlaps_TPM_COUNTS", today_f, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)

















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






