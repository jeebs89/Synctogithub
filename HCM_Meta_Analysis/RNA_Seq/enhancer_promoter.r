
        library('tidyverse', lib="/home/jcayford/r_libs")
        library('viridisLite')
      
  
 # Functions
                # Quick head() with the amount of rows
                    w <- function(x){print(head(x));nrow(x)}
                # Rounding data
                    round_to <- function(x, to = window_length) round(x/to)*to
                # Importing the chip files. Rounding the peaks
                    chip_import <- function(import_df, round){
                        df <- read.table(import_df)
                        data.frame("chr1"=df[,1], "peak"=round_to(df[,2], to=round))
                    }
                # Exporting the loops for washU visualization
                    final_loop_export <- function(input_df, disease_state){
                        write.table(input_df, paste0("hic_contacts_", disease_state, "_merged_q", filter_qval, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
                    }

                # Overlap function for enhancers
                    overlaps_func_mid_enhancers <- function(input_df, df_b, mid, row){
                        if(mid == 1){mids <- 2}
                        if(mid == 2){mids <- 4}
                            df_subs <- data.frame(input_df, "mid"=abs(input_df[, mids]-df_b[row,2]), "enhancer"=df_b[row, 2], "gene"=df_b[row,4])
                            df_subs %>% filter(mid == 0 | mid == 5000)
                    }

                # Selecting the condition. This is for post overlaps. This will give a common condition which is not specific to only 1 selected condition
                    status_select <- function(input_df, condition){
                        input_df %>% filter(cond==condition | cond == "all") %>% select(c(1:6))
                    }

                # Overlap function for everything except enhancers
                    overlaps_func_mid <- function(input_df, df_b, mid, row){
                        if(mid == 1){mids <- 2}
                        if(mid == 2){mids <- 4}
                            df_subs <- data.frame(input_df, "mid"=abs(input_df[, mids]-df_b[row,2]))
                            df_subs %>% filter(mid == 0 | mid == 5000)
                    }

                # Function to complete the overlaps with an input with 2 columns (chr and position)
                    chip_all_overlaps_function <- function(comparing_df){
                        output <- NULL
                            for(j in 1:19){
                                    print(paste0(Sys.time(), ".........starting chromosome: ", j))
                                    fithic_chr <- fithic %>% filter(chr1==paste0("chr", j))
                                    compared_df <- comparing_df %>% filter(chr1==paste0("chr", j))

                                                overlaps_mid1 <- NULL
                                                overlaps_mid2 <- NULL
                                                for(i in 1:nrow(compared_df)){
                                                    overlaps_mid1 <- rbind(overlaps_mid1, overlaps_func_mid(fithic_chr, compared_df, 1, i))
                                                    overlaps_mid2 <- rbind(overlaps_mid2, overlaps_func_mid(fithic_chr, compared_df, 2, i))
                                                }
                                                overlaps_a <- rbind(overlaps_mid1, overlaps_mid2)
                                                overlaps <- overlaps_a %>% distinct(chr1, fragmentMid1, fragmentMid2, .keep_all=TRUE)
                                                output <- rbind(output, overlaps)
                            }
                        output
                    }

               
#

            
       

# Setting the directories 
    heart_directory <- "/home/jcayford/HCM/Analysis/heart_info_BingRen_2012 nature"
    fithic_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    export_directory <- "/home/jcayford/HCM/Analysis/HiC_R/washU_outputs/Enhancers.promoters"
    filter_qval <- 0.005

# Uploading the heart specific data from Bing Ren's 2012 Nature paper
    setwd(heart_directory)
        enhancers_a <- read.table("enhancer_tss_heart.csv", header=TRUE, sep=",")
        enhancers_b <- data.frame("chr1"=enhancers_a[,1], round_to(enhancers_a[,2:3], to=2500), enhancers_a[,4:5])
        enhancers <- enhancers_b %>% distinct(chr1, enhancer, TSS, .keep_all=TRUE)

        ctcf_peaks <- chip_import("heart.ctcf.peak.txt", 2500)
        pol2 <- chip_import("heart.polII.peak.txt", 2500)
        k27ac <- chip_import("heart.h3k27ac.peak.txt", 2500)  
        promoters <- chip_import("heart.h3k4me3.peak.txt", 2500)  

# Uploading the FitHiC merged file. This has a q value of 0.005
    setwd(fithic_directory)
    fithic <- read.table("hic_contacts_mef2a_all_q0.01.txt", header=TRUE)




# Creation of the enhancer overlaps    
    all_chr_enhancers <- NULL
    for(j in 1:19){
            print(paste0(Sys.time(), ".........starting chromosome: ", j))
            fithic_chr <- fithic %>% filter(chr1==paste0("chr", j))
            enhancers_chr <- enhancers %>% filter(chr1==paste0("chr", j)) %>% arrange(enhancer)

                        overlaps_mid1 <- NULL
                        overlaps_mid2 <- NULL
                        for(i in 1:nrow(enhancers_chr)){
                            overlaps_mid1 <- rbind(overlaps_mid1, overlaps_func_mid_enhancers(fithic_chr, enhancers_chr, 1, i))
                            overlaps_mid2 <- rbind(overlaps_mid2, overlaps_func_mid_enhancers(fithic_chr, enhancers_chr, 2, i))
                        }
                        overlaps_a <- rbind(overlaps_mid1, overlaps_mid2)
                        overlaps <- overlaps_a %>% distinct(chr1, fragmentMid1, fragmentMid2, .keep_all=TRUE)
                        all_chr_enhancers <- rbind(all_chr_enhancers, overlaps)
    }

# Selecting the status specific loops
    disease_loops <- status_select(all_chr_overlaps, "disease")
    healthy_loops <- status_select(all_chr_overlaps, "healthy")
    ctcf_ko_loops <- status_select(all_chr_overlaps, "ctcf_ko")

# Exporting the loops for use in washU
    setwd(export_directory)
    final_loop_export(disease_loops, "disease_enhancers")
    final_loop_export(healthy_loops, "healthy_enhancers")
    final_loop_export(ctcf_ko_loops, "ctcf_ko_enhancers")


# Checking to see which genes are specific to healthy or disease
    h_loops <- all_chr_overlaps %>% filter(cond=="healthy") %>% select(c(6:9, 19)) %>% arrange(q.value)
    d_loops <- all_chr_overlaps %>% filter(cond=="disease") %>% select(c(6:9, 19)) %>% arrange(desc(d.contactCount)) %>% distinct(gene, .keep_all=TRUE)



# Generation of the chip peaks and the promoters
    all_chr_k27 <- chip_all_overlaps_function(k27ac)
    all_chr_promoters <- chip_all_overlaps_function(promoters)
    all_chr_ctcf <- chip_all_overlaps_function(ctcf_peaks)
    all_chr_pol2 <- chip_all_overlaps_function(pol2)







all_e <- all_chr_enhancers
all_p <- data.frame(all_chr_promoters, "enhancer"="promoter", "gene"="promoter")

all_chr_merged <- rbind(all_e, all_p)

contact_thresh <- 3


   subs <- data.frame(all_chr_merged, "sub_a"=(all_chr_merged$h.contactCount - all_chr_merged$d.contactCount),
                                                    "sub_b"=(all_chr_merged$h.contactCount - all_chr_merged$c.contactCount),     
                                                        "sub_c"=(all_chr_merged$d.contactCount - all_chr_merged$c.contactCount), "cond"=NA)

            # Determining specific loops based on the subtraction of the different conditions
                subs[(subs$sub_a < -contact_thresh & subs$sub_c > contact_thresh), 19] <- "disease"
                subs[(subs$sub_a > contact_thresh & subs$sub_b > contact_thresh), 19] <- "healthy"
                subs[(subs$sub_b < -contact_thresh & subs$sub_c < -contact_thresh), 19] <- "ctcf_ko"
                subs[is.na(subs)] <- "all"




  final_loop_export <- function(input_df, prom_ench, disease_state, export){   
    if(prom_ench == "enhancer"){output_df <- input_df %>% filter(gene != prom_ench, cond==disease_state) %>% select(c(1:6))}
    if(prom_ench == "promoter"){output_df <- input_df %>% filter(gene == prom_ench, cond==disease_state) %>% select(c(1:6))}

    output_df <- output_df %>% arrange(chr1, fragmentMid1, fragmentMid2)
    
    if(export=="y" | export=="yes"){write.table(output_df, paste0("hic_contacts_", disease_state, "_", prom_ench, "_q", filter_qval, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)}
    if(export=="n" | export=="no"){output_df}
}


  
setwd(export_directory)
final_loop_export(subs, "enhancer", "healthy", "yes")
final_loop_export(subs, "promoter", "healthy", "yes")

final_loop_export(subs, "enhancer", "disease", "yes")
final_loop_export(subs, "promoter", "disease", "yes")

final_loop_export(subs, "enhancer", "ctcf_ko", "yes")
final_loop_export(subs, "promoter", "ctcf_ko", "yes")




test <- final_loop_export(subs, "enhancer", "disease", "no")





# End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
# From there you can import them into a .json for washU visualization
    awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_ctcf_ko_enhancers_merged_q0.005.txt | sort -k1,1 -k2,2n > hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
    bgzip hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
    tabix -p bed hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed.gz










# End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
# Combined values... From there you can import them into a .json for washU visualization
    awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","$NF"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","$NF"\t"(NR-1)"\t."}}}' hic_contacts_JC_TEST_merged_q0.005.txt | sort -k1,1 -k2,2n > hic_contacts_JC_TEST_merged_q0.005_out_washU.bed
    bgzip hic_contacts_JC_TEST_merged_q0.005_out_washU.bed
    tabix -p bed hic_contacts_JC_TEST_merged_q0.005_out_washU.bed.gz


cp *_out_wash* /mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/washU_filtered/



h_p <- read.table("hic_contacts_healthy_promoter_q0.005.txt", header=TRUE)
h_e <- read.table("hic_contacts_healthy_enhancer_q0.005.txt", header=TRUE)



h_p_a <- data.frame(h_p[,1:4], "score"=50)
h_p_b <- data.frame(h_e[,1:4], "score"=10)

test <- rbind(h_p_b, h_p_a)
test2 <- test %>% arrange(chr1, fragmentMid1, fragmentMid2) %>% distinct(chr1, fragmentMid1, fragmentMid2, .keep_all=TRUE)

write.table(test2, paste0("hic_contacts_", "JC_TEST_merged_q", filter_qval, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)

