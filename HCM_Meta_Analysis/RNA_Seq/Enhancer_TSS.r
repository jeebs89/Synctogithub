
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
                            df_subs %>% filter(mid == 0 | mid == 5000 | mid == 10000)
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
                            df_subs %>% filter(mid == 0 | mid == 5000 | mid == 10000)
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

                # Selecting and outputting a file. CAUTION USES 2 INDEPENDENT FUNCTIONS!
                    select_output_func <- function(df, status, out_file){
                        df_loops <- status_select(all_chr_overlaps, status)
                        setwd(export_directory)
                        final_loop_export(df_loops, out_file)
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
    fithic <- read.table("hic_contacts_mef2a_all_merged_q0.005.txt", header=TRUE)




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



select_output_func(all_chr_enhancers, "disease", "disease_enhancers")
select_output_func(all_chr_enhancers, "healthy", "healthy_enhancers")
select_output_func(all_chr_enhancers, "ctcf_ko", "ctcf_ko_enhancers")

select_output_func(all_chr_promoters, "disease", "disease_promoters")
select_output_func(all_chr_promoters, "healthy", "healthy_promoters")
select_output_func(all_chr_promoters, "ctcf_ko", "ctcf_ko_promoters")









# End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
# From there you can import them into a .json for washU visualization
    awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_ctcf_ko_enhancers_merged_q0.005.txt | sort -k1,1 -k2,2n > hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
    bgzip hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
    tabix -p bed hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed.gz


