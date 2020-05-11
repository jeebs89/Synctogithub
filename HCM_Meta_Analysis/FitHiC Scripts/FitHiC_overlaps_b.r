
### Part 1 of the script - to generate the filtered HiC data with a given ChIP-Seq track
            library('tidyverse', lib="/home/jcayford/r_libs")
            library('viridisLite')
                
            # Functions
                    # Rounding data
                        round_to <- function(x, to = window_length) round(x/to)*to
                    # Quick head() with the amount of rows
                        w <- function(x){print(head(x));nrow(x)}
                    # Changing the import to a numeric value
                        df_c<- function(data_frame, col_num){
                            as.numeric(as.character(data_frame[,col_num]))
                        } 
                    # Filter for ChIP windows. The value is distance from the ChIP peak to find a HiC contact. This does both contacts
                        chip_filter <- function(sub_value){
                            df <- data.frame("chr"=chip_peaks[,1], (chip_means+(sub_value)), chip_peaks[,c(4:5)])
                            
                            colnames(df)[1:2] <- c("chr1", "fragmentMid1")
                            a <- inner_join(hic_data, df, by=(c("chr1", "fragmentMid1")))
                            a1 <- data.frame(a, "ChIP_contactMid"=1)

                            colnames(df)[1:2] <- c("chr1", "fragmentMid2")
                            b <- inner_join(hic_data, df, by=(c("chr1", "fragmentMid2")))
                            b1 <- data.frame(b, "ChIP_contactMid"=2)

                            c <- rbind(a1, b1)
                            c %>% distinct(fragmentMid1, fragmentMid2,  .keep_all=TRUE)
                        }    
                
            #

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


                # Directories for importing and exporting the data 
                        hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes"
                        chip_beds_directory <- "/home/jcayford/HCM/Analysis/ChIP_files/"            
                        export_directory <- "/home/jcayford/HCM/Analysis/HiC_R/washU_outputs"


                    #chip_bed <- "polr2a.bed" 
                    chip_bed <- "mef2a.bed"
                    tf <- "mef2a"  
                    qval <- 0.1
                    filter_qval <- 0.05


                # Start of the script

                    # Loading and pre-processing of data
            
                        # ChIP Peaks
                            setwd(chip_beds_directory)
                            chip_peaks_a <- read.table(chip_bed)
                            chip_peaks_a <- chip_peaks_a %>% arrange(V1, V2, V3)
                                    if(chip_peaks_a[1,1]==1){chip_peaks <- data.frame(paste("chr", chip_peaks_a[,1], sep=""), chip_peaks_a[,2:5])}
                                    if(chip_peaks_a[1,1]=="chr1"){chip_peaks <- chip_peaks_a[,1:5]}                       
                            colnames(chip_peaks) <- c("chr1", "start", "end", "peak", "height")     

                        # FitHiC data
                            setwd(hic_data_directory)
                            hic_data_a <- read.table(paste("FitHiC_all_qval", qval, ".txt", sep=""), header=TRUE)
                            hic_data <- data.frame(hic_data_a[,1], df_c(hic_data_a, 2), hic_data_a[,3],
                                                df_c(hic_data_a, 4), df_c(hic_data_a, 5), df_c(hic_data_a, 6),
                                                df_c(hic_data_a, 7), df_c(hic_data_a, 8), df_c(hic_data_a, 9))               
                            colnames(hic_data) <- colnames(hic_data_a)           
                            hic_data <- hic_data %>% filter(q.value < filter_qval)


                    # Filtering FitHiC data with the ChIP-Seq Peaks

                        chip_means <- round_to(rowMeans(chip_peaks[,2:3]), to=2500)
                
                        # Filtering the HiC data
                            hic_filtered_a <- NULL
                            for(i in c((-5000), (-2500), 0, 2500, 5000)){
                                hic_filtered_a <- rbind(hic_filtered_a, chip_filter(i))
                            }                
                            hic_filtered <- hic_filtered_a %>% distinct(fragmentMid1, fragmentMid2,  .keep_all=TRUE)




### Part 2 to create the BED file to be merged by bedtools

    
        bed_create <- function(input_df, frag_mid, extension){
            if(frag_mid==1){df_a <- data.frame(input_df[,1:2], (input_df[,2]+extension))}
            if(frag_mid==2){df_a <- data.frame(input_df[,c(1,3)], (input_df[,3]+extension))} 

            colnames(df_a) <- c("chr", "start", "end")
            df_a %>% arrange(chr, start, end)
        }

        export_bed <- function(input_df, export_name){
            write.table(input_df, paste0(export_name, ".bed"), col.names=FALSE, row.names=FALSE,  quote=FALSE, sep="\t")
        }



        hic_bed <- hic_filtered[,c(1, 2, 4)]
        hic <- hic_bed %>% arrange(chr1, fragmentMid1, fragmentMid2)


        mid1_bed <- bed_create(hic, 1, 10000)
        export_bed(mid1_bed, "merge_mid1")



        ### Quit R and generate the merged bed file 
            bedtools merge -d 1000 -i merge_mid1.bed > mid1_merged.bed








#### Part 3 - Merging of the files with a given bin - default in the script is 2 bins (10kb)


chr_select <- function(input_df, chromosome){
    input_df %>% filter(chr1==chromosome)
}


mid2_overlaps_func <- function(input_df, row){
    df_subs <- data.frame(input_df, "mid2"=abs(input_df[,4]-input_df[row,4]))
    df_filter <- df_subs %>% filter(mid2 ==0 | mid2 == 5000 | mid2 == 10000)

    new_val <- cbind(
                        df_filter[1,1:3],
                        median(df_filter$fragmentMid2),
                        min(df_filter$p.value),
                        min(df_filter$q.value),
                        mean(df_filter$h_contactCount),
                        mean(df_filter$d_contactCount),
                        df_filter[1, 9:12]
                    )
    colnames(new_val) <- colnames(hic_chr)
    new_val
}






frag1 <- read.table("mid1_merged.bed"); colnames(frag1) <- c("chr1", "start", "end")
hic_filtered[is.na(hic_filtered)] <- 0


all_chr_merged <- NULL
for(j in 1:19){
    print(paste0(Sys.time(), "...........working on chromosome:", j))
    hic_chr <- chr_select(hic_filtered, paste0("chr", j))
    frag1_chr <- chr_select(frag1, paste0("chr", j))


    merged_data <- NULL
    for(i in 1:nrow(frag1_chr)){
        initial_filter <- hic_chr %>% filter(fragmentMid1 >= frag1_chr[i,2] & fragmentMid1 <= frag1_chr[i,3])

            if(nrow(initial_filter)==1){new_val <- initial_filter}
            if(nrow(initial_filter)>1){

                    mid2_arrange <- initial_filter %>% arrange(fragmentMid2)
                    mid2_arrange$fragmentMid1 <- median(mid2_arrange$fragmentMid1)

                    new_vals_raw <- NULL
                                for(k in 1:nrow(mid2_arrange)){
                                    new_vals_raw <- rbind(new_vals_raw, mid2_overlaps_func(mid2_arrange, k))
                                }
                    new_val <- new_vals_raw %>% distinct(fragmentMid2, .keep_all=TRUE)
            }
            merged_data <- rbind(merged_data, new_val)
    }

    all_chr_merged <- rbind(all_chr_merged, merged_data)
}














subs <- data.frame(all_chr_merged, "sub"=(all_chr_merged$h_contactCount - all_chr_merged$d_contactCount), "cond_1"=NA)

subs[subs$sub < -3, 14] <- "disease"
subs[subs$sub > 3, 14] <- "healthy"
subs[subs$sub >= -3 & subs$sub <= 3, 14] <- "both"

disease_loops <- subs %>% filter(cond_1 == "disease" | cond_1 == "both") %>% select(c(1:6))
healthy_loops <- subs %>% filter(cond_1 == "healthy" | cond_1 == "both") %>% select(c(1:6))



setwd(export_directory)

write.table(disease_loops, "hic_contacts_mef2_disease_combo_merged_q0.05.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(healthy_loops, "hic_contacts_mef2_healthy_combo_merged_q0.05.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)




awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_mef2_healthy_combo_merged_q0.05.txt| sort -k1,1 -k2,2n > hic_contacts_mef2_healthy_combo_merged_q0.05_out_washU.bed
bgzip hic_contacts_mef2_healthy_combo_merged_q0.05_out_washU.bed
tabix -p bed hic_contacts_mef2_healthy_combo_merged_q0.05_out_washU.bed.gz




                


