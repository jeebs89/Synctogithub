# Mustache

library('tidyverse', lib="/home/jcayford/r_libs")


# Funcations
    # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}

    # Cleaning up the data and getting ready for export into a json readable file. 
        # State is either "healthy" or "disease"
            mustache_data_prep <- function(state){
                    df <- read.table(paste0("mustache_", state, "_all_out_0.1.tsv"), header=TRUE, stringsAsFactors=F)
                    df_q_filtered <- df %>% filter(FDR < qval_threshold) %>% arrange(FDR)
                    
                    df_arranged <- df_q_filtered %>% arrange(BIN1_CHR, BIN1_START, BIN1_END, BIN2_CHROMOSOME, BIN2_START, BIN2_END)
                    df_final <- data.frame(df_arranged[,1], rowMeans(df_arranged[,2:3]), df_arranged[,4], rowMeans(df_arranged[,5:6]), df_arranged[,7], (10^(-df_arranged[,8])))
                    colnames(df_final) <- cols
                    df_final

            }
    # Exporting the loops for washU visualization
            final_loop_export <- function(input_df, disease_state, date){
                write.table(input_df, paste0("mustache_hic_contacts_", "all_", disease_state, "_non_merged_q", qval_threshold, "_", date, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
            }

    
#

# Setting the working directories 
    main_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files"
    disease_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files/disease/"
    healthy_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files/healthy/"
    export_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files/json_file_prep"

#
# q.value threshold, since the input files have a threshold of 0.1
qval_threshold <- 0.1
today <- Sys.Date()
date <- format(today, format="%d_%b_%Y")


# Column Names for export
cols <- c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "p.value", "q.value")


    #
    # Importing the files and filtering based on the given threshold for the qval
        setwd(disease_wd)
            disease_df <- mustache_data_prep("disease")
        setwd(healthy_wd)
            healthy_df <- mustache_data_prep("healthy")
            
    # Data Export
        setwd(export_wd)
            final_loop_export(disease_df, "disease", date)
            final_loop_export(healthy_df, "healthy", date)
    

  # End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
    # From there you can import them into a .json for washU visualization
    # conda activate samtools
        awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' mustache_hic_contacts_all_healthy_non_merged_q0.1_07_Oct_2020.txt | sort -k1,1 -k2,2n > mustache_hic_contacts_all_healthy_non_merged_q0.1_07_Oct_2020_out_washU.bed
        bgzip mustache_hic_contacts_all_healthy_non_merged_q0.1_07_Oct_2020_out_washU.bed
        tabix -p bed mustache_hic_contacts_all_healthy_non_merged_q0.1_07_Oct_2020_out_washU.bed.gz


    # Copy the data to the BioAdHoc folder
        cp *_07_Oct_2020* /mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/washU_filtered/


