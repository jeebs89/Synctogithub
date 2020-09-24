
# Library import(s)
    library(tidyverse)

# Functions
        # Quick head() with the amount of rows
            w <- function(x){print(head(x));nrow(x)}

    
        # Importing .gtf file and selecting only chromosomes 1-19                       
                gtf_import <- function(file_name){
                    input <- read.table(file_name, header=FALSE, sep='\t')
                    pre <- input %>% mutate(V1=as.numeric(V1))
                    pre %>% filter(V1<=19)
                }


        # Splitting the last column of the .gtf import file and changing to tpms.
                gtf_cleanup<-function(input_df){
                    # Transcripts
                        filtered_df_transcript <- input_df %>% filter(V3=="transcript")
                        split_transcript <- str_replace_all(filtered_df_transcript[,9], ";", "")
                        combined_df_transcript <- data.frame(filtered_df_transcript[,1:8], "V9"=split_transcript)

                        all_split_df_transcript <- combined_df_transcript %>% separate("V9", paste("V", c(9:22), sep=""), " ")
                        col_select_transcript <- all_split_df_transcript[, c(1, 3, 4, 5, 10, 12, 14, 16, 18, 20, 22)]
                        colnames(col_select_transcript)<-c("chr", "type", "start", "end", "gene", "trans_exon", "FPKM", "frac", "conf_lo", "conf_hi", "cov")
                        
                    # Exons    
                        filtered_df_exon <- input_df %>% filter(V3=="exon")
                        split_exon <- str_replace_all(filtered_df_exon[,9], ";", "")
                        combined_df_exon <- data.frame(filtered_df_exon[,1:8], "V9"=split_exon)

                        all_split_df_exon <- combined_df_exon %>% separate("V9", paste("V", c(9:24), sep=""), " ")
                        col_select_exon <- all_split_df_exon[, c(1, 3, 4, 5, 10, 12, 16, 18, 20, 22, 24)]
                        colnames(col_select_exon)<-c("chr", "type", "start", "end", "gene", "trans_exon", "FPKM", "frac", "conf_lo", "conf_hi", "cov")  
                    
                    # Generating final df with TPMs    
                        cols <- rbind(col_select_transcript, col_select_exon)
                        fpkms <- as.numeric(as.character(cols$FPKM))       
                        tmp_df <- data.frame(cols, "TPM"=fpkms*((1e6)/sum(fpkms)))
                        tmp_df %>% arrange(chr, start, end, gene)
                }

        # Completeing a full_join on 3 data.frames. This is specific for the gtf_cleanup files. 
                double_fulljoin <- function(df_1, df_2, df_3){
                    names<-c("chr", "type", "start", "end", "gene", "trans_exon")
                    join_1 <- full_join(df_1[,c(1:6, 12)], df_2[,c(1:6, 12)], by=names)
                    join_2 <- full_join(join_1, df_3[,c(1:6, 12)], by=names)
                }

#


rna_seq_data <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/data"


# Script
    setwd(rna_seq_data)

        # Healthy RNA_Seq import and join
            h_rna_a_pre <- gtf_import('GSM2538380_Control_RNAseq_Rep1_transcripts.gtf')
            h_rna_b_pre <- gtf_import('GSM2538381_Control_RNAseq_Rep2_transcripts.gtf')
            h_rna_c_pre <- gtf_import('GSM2538382_Control_RNAseq_Rep3_transcripts.gtf');gc()

            h_rna_a <- gtf_cleanup(h_rna_a_pre);h_rna_b <- gtf_cleanup(h_rna_b_pre);h_rna_c <- gtf_cleanup(h_rna_c_pre)

            h_rna_all <- double_fulljoin(h_rna_a, h_rna_b, h_rna_c)
            colnames(h_rna_all)[7:9]<-c("h1_TPM", "h2_TPM", "h3_TPM")

        # Disease RNA_Seq import and join
            d_rna_a_pre <- gtf_import('GSM2538387_TAC_RNAseq_Rep1_transcripts.gtf')
            d_rna_b_pre <- gtf_import('GSM2538388_TAC_RNAseq_Rep2_transcripts.gtf')
            d_rna_c_pre <- gtf_import('GSM2538389_TAC_RNAseq_Rep3_transcripts.gtf');gc()

            d_rna_a <- gtf_cleanup(d_rna_a_pre);d_rna_b <- gtf_cleanup(d_rna_b_pre);d_rna_c <- gtf_cleanup(d_rna_c_pre)

            d_rna_all <- double_fulljoin(d_rna_a, d_rna_b, d_rna_c)
            colnames(d_rna_all)[7:9]<-c("d1_TPM", "d2_TPM", "d3_TPM")



        # Combining the healthy and disease data.frames
            rna_all <- full_join(h_rna_all, d_rna_all, by=c("chr", "type", "start", "end", "gene", "trans_exon"))
        # Calculating all the means and selecting rows which both samples have a read in (TPM > 0)
            rna_means<-data.frame(rna_all, "h_mean"=rowMeans(rna_all[,7:9]), "d_mean"=rowMeans(rna_all[,10:12]))
            rna_filter <- rna_means %>% filter(h_mean>0, d_mean>0)
                        
        # Calculating pvals for all rows for the healthy vs disease TPMs
            pvals <- rna_filter %>%
                        rowwise() %>%
                        mutate(pval = t.test(c(h1_TPM, h2_TPM, h3_TPM), c(d1_TPM, d2_TPM, d3_TPM))$p.value)

        # Exporting the file with all the rows 
            rna_with_pval<-data.frame(pvals)
            setwd(output_directory)
            write.table(rna_with_pval, "rna_analysis_all_windows.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)


#

