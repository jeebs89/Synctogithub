
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



            
       


    heart_directory <- "/home/jcayford/HCM/Analysis/heart_info_BingRen_2012 nature"
    fithic_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    export_directory <- "/home/jcayford/HCM/Analysis/HiC_R/washU_outputs/Enhancers.promoters"
    filter_qval <- 0.005

setwd(heart_directory)
    enhancers_a <- read.table("enhancer_tss_heart.csv", header=TRUE, sep=",")
    enhancers_b <- data.frame("chr1"=enhancers_a[,1], round_to(enhancers_a[,2:3], to=2500), enhancers_a[,4:5])

    enhancers <- enhancers_b %>% distinct(chr1, enhancer, TSS, .keep_all=TRUE)

    ctcf_peaks <- chip_import("heart.ctcf.peak.txt", 2500)
    pol2 <- chip_import("heart.polII.peak.txt", 2500)
    k27ac <- chip_import("heart.h3k27ac.peak.txt", 2500)  
    promoters <- chip_import("heart.h3k4me3.peak.txt", 2500)  

setwd(fithic_directory)
   
   fithic <- read.table("hic_contacts_mef2a_all_merged_q0.005.txt", header=TRUE)






                    overlaps_func_mid_enhancers <- function(input_df, df_b, mid, row){
                        if(mid == 1){mids <- 2}
                        if(mid == 2){mids <- 4}
                            df_subs <- data.frame(input_df, "mid"=abs(input_df[, mids]-df_b[row,2]), "enhancer"=df_b[row, 2], "gene"=df_b[row,4])
                            df_subs %>% filter(mid == 0 | mid == 5000 | mid == 10000)
                    }



all_chr_overlaps <- NULL
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
                    all_chr_overlaps <- rbind(all_chr_overlaps, overlaps)
}


disease_loops <- all_chr_overlaps %>% filter(cond=="disease" | cond == "all") %>% select(c(1:6))
healthy_loops <- all_chr_overlaps %>% filter(cond=="healthy" | cond == "all") %>% select(c(1:6))
ctcf_ko_loops <- all_chr_overlaps %>% filter(cond=="ctcf_ko" | cond == "all") %>% select(c(1:6))




                setwd(export_directory)
                final_loop_export(disease_loops, "disease_enhancers")
                final_loop_export(healthy_loops, "healthy_enhancers")
                final_loop_export(ctcf_ko_loops, "ctcf_ko_enhancers")

    # End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
    # From there you can import them into a .json for washU visualization
        awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_ctcf_ko_enhancers_merged_q0.005.txt | sort -k1,1 -k2,2n > hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
        bgzip hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed
        tabix -p bed hic_contacts_ctcf_ko_enhancers_merged_q0.005_out_washU.bed.gz






h_loops <- all_chr_overlaps %>% filter(cond=="healthy") %>% select(c(6:9, 19)) %>% arrange(q.value)
d_loops <- all_chr_overlaps %>% filter(cond=="disease") %>% select(c(6:9, 19)) %>% arrange(d.contactCount) %>% distinct(gene, .keep_all=TRUE)



d_loops %>% filter(gene=="Zbtb16")




 #### NEED TO DO THE PROMOTERS AS WELL.... Link both of them together and see what the results are like.. 





                    overlaps_func_mid <- function(input_df, df_b, mid, row){
                        if(mid == 1){mids <- 2}
                        if(mid == 2){mids <- 4}
                            df_subs <- data.frame(input_df, "mid"=abs(input_df[, mids]-df_b[row,2]))
                            df_subs %>% filter(mid == 0 | mid == 5000 | mid == 10000)
                    }



all_chr_overlaps_promoters <- NULL
for(j in 1:19){
        print(paste0(Sys.time(), ".........starting chromosome: ", j))
        fithic_chr <- fithic %>% filter(chr1==paste0("chr", j))
        promoters_chr <- promoters %>% filter(chr1==paste0("chr", j))

                    overlaps_mid1 <- NULL
                    overlaps_mid2 <- NULL
                    for(i in 1:nrow(enhancers_chr)){
                        overlaps_mid1 <- rbind(overlaps_mid1, overlaps_func_mid(fithic_chr, promoters_chr, 1, i))
                        overlaps_mid2 <- rbind(overlaps_mid2, overlaps_func_mid(fithic_chr, promoters_chr, 2, i))
                    }
                    overlaps_a <- rbind(overlaps_mid1, overlaps_mid2)
                    overlaps <- overlaps_a %>% distinct(chr1, fragmentMid1, fragmentMid2, .keep_all=TRUE)
                    all_chr_overlaps_promoters <- rbind(all_chr_overlaps_promoters, overlaps)
}

w(all_chr_overlaps_promoters)
w(all_chr_overlaps)

test <- data.frame(all_chr_overlaps_promoters, "enhancer"=NA, "gene"=NA)
test2 <- rbind(all_chr_overlaps, test)
test3 <- test2 %>% distinct(chr1, fragmentMid1, fragmentMid2, .keep_all=TRUE) %>% arrange(chr1, fragmentMid1, fragmentMid2)