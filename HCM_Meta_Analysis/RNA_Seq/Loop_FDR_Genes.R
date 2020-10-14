
library('tidyverse', lib="/home/jcayford/r_libs")
library('viridisLite')
    




    # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}                
            # Getting the loops with genes
                loop_rna_overlaps_func <- function(input_a, input_b, row){
                        loops <- input_a %>% filter(as.numeric(fragmentMid1) < as.numeric(input_b$transcript_start[row]) & as.numeric(fragmentMid2) > as.numeric(input_b$transcript_end[row]))
                        if(nrow(loops)==0){NULL}
                        if(nrow(loops)>0){cbind(loops, input_b[row,])}
                    }           # Calculating z-scores
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





# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")




    hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"

    hic_data_directory_2 <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    deseq_const_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq/Loop_Constraints"

    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
    tpm_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/rosagarrido_rnaseq/4.Output/counts"


 

qval <- 0.01
filter_qval <- 0.01
contact_thresh <- 3
loop_constraint <- c(0, 1.5, 1, 0.5)



    # Start of the script

        # Loading and pre-processing of data

       


# Importing the filtered MEF2 HiClooping data
    setwd(hic_data_directory_2)
    all_loops_a <- read.table("hic_contacts_mef2a_all_merged_q0.005_28sep2020.txt", header=T)
    all_loops <- all_loops_a %>% filter(cond=="disease" | cond=="healthy" | cond=="h_d") %>% select(c(1:8, 10:12,16))
    loop_size <- data.frame(all_loops, "loop_size"=(all_loops[,4]-all_loops[,2]))
    all_size <- loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(loop_size)))

# Importing the filtered Loops based on RNA-Seq
    setwd(deseq_const_wd)
    filtered_loop_size_a <- read.table("All_loop_rna_overlaps_master_0mb_limit_14_Oct_2020.txt", header=TRUE)
    filtered_loop_size <- data.frame(filtered_loop_size_a, "loop_size"=(filtered_loop_size_a[,4]-filtered_loop_size_a[,2]))
    all_t_size <- filtered_loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(filtered_loop_size)))



# Importing the whole HiC Dataset
    setwd(hic_data_directory)
    hic_data_a <- read.table(paste("FitHiC_all_qval", qval, ".txt", sep=""), header=TRUE)
    hic_data <- hic_data_a %>% select(c(1:8)) 
                
                




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


#















 # FitHiC data






library('tidyverse', lib="/home/jcayford/r_libs")
library('viridisLite')
library('readr')
library('pheatmap')
library('gplots')








        

    #


 

# Importing the randomized HiC Data


          
    loop_constraint<-1

    rand_hic_loops_a <- sample_n(hic_data, 5938, replace=FALSE)


    if(loop_constraint==0){rand_hic_loops <- rand_hic_loops_a}
    if(loop_constraint>0){
        loop_size_limit <- loop_constraint*1000000     
        rand_hic_loops_b <- data.frame(rand_hic_loops_a, "loop_size"=(rand_hic_loops_a$fragmentMid2 - rand_hic_loops_a$fragmentMid1))
        rand_hic_loops <- rand_hic_loops_b %>% filter(loop_size < loop_size_limit) %>% select(c(1:8))

    }

    nrow(rand_hic_loops)

    gene_rna_loops_rand <- NULL
    for(j in 1:19){
        #print(paste(Sys.time(), "..........matching the loops for chromosome: ", j,sep=""))

        all_loops_1 <- rand_hic_loops %>% filter(chr1==paste0("chr", j)) %>% arrange(fragmentMid1, fragmentMid2)
        all_genes_1 <- all_genes %>% filter(chr1==paste0("chr", j))

            chr_genes_rand <- NULL
            for(i in 1:nrow(all_genes_1)){
                chr_genes_rand <- rbind(chr_genes_rand, loop_rna_overlaps_func(all_loops_1, all_genes_1, i))
            }
        gene_rna_loops_rand <- rbind(gene_rna_loops_rand, chr_genes_rand)

    }

    gene_rna_loops_rand_b <- gene_rna_loops_rand[c(1, 10:15, 9, 2, 4, 16:17, 18:25)]
    uniq_trans_loops_rand <- gene_rna_loops_rand_b %>% distinct(Gene, .keep_all=TRUE)

    # Filtering using the standard deviations to be less than 0.6
    u_t_loops_full_rand <- uniq_trans_loops_rand %>% arrange(padj) %>% distinct(mgi_symbol, .keep_all=TRUE) %>% select(c(1,2,5,7,9:10, 13:20))
    u_t_loops_trim_rand <- data.frame(u_t_loops_full_rand[, c(2, 4, 6:14)])

    sds_rand <- NULL
    for(i in 1:nrow(u_t_loops_trim_rand)){
        sds_rand <- rbind(sds_rand, filtering_rna(i, u_t_loops_trim_rand))
    }

    loop_sds_rand <- data.frame(u_t_loops_trim_rand, sds_rand)
    loops_sds_filtered_rand <- loop_sds_rand %>% filter(sd_h < 0.6 & sd_d < 0.6) %>% arrange(padj)
    loops_final_rand <- data.frame(loops_sds_filtered_rand[,2:ncol(loops_sds_filtered_rand)], row.names=loops_sds_filtered_rand[,1])



data.frame("input_loops"=nrow(rand_hic_loops), "trans_with_loops"=nrow(uniq_trans_loops_rand), "genes_passing_filter"=nrow(loops_final_rand))




# Meta data information
    meta_data <- data.frame("State"=c(rep("Healthy", 3), rep("HCM", 3)), row.names=colnames(loops_final_rand[5:10]))
    Var1 <- c("dodgerblue", "firebrick3")
    names(Var1) <- c("Healthy", "HCM")
    anno_cols <- list(State = Var1)


colors <- colorRampPalette(c("#1F00EB", "black", "#FFFB0F"))



pheatmap(loops_final_rand[,5:10], cluster_cols=T, 
        cluster_rows=T, color = colors(75),
        border_color=NA, clustering_method="complete", show_rownames=TRUE, show_colnames=FALSE,
        main="Healthy vs HCM Heatmap", #cellwidth = 35, #cellheight = 7.5,
        annotation_col=meta_data, annotation_colors=anno_cols
)
























    rand_hic_loops <- sample_n(hic_data, nrow(all_loops), replace=FALSE)


        all_loops_1 <- rand_hic_loops %>% filter(chr1==paste0("chr", 3)) %>% arrange(fragmentMid1, fragmentMid2)
        all_genes_1 <- all_genes %>% filter(chr1==paste0("chr", 3))


        loop_rna_overlaps_func(all_loops_1, test, 70)

            chr_genes_rand <- NULL
            for(i in 1:nrow(all_genes_1)){
                chr_genes_rand <- rbind(chr_genes_rand, loop_rna_overlaps_func(all_loops_1, all_genes_1, i))
            }
        gene_rna_loops_rand <- rbind(gene_rna_loops_rand, chr_genes_rand)

    }

