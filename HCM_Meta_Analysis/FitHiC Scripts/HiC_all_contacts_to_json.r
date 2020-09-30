
# Generating HiC for json... for disease and healthy... all loops included




library('tidyverse', lib="/home/jcayford/r_libs")
library('viridisLite')
    
       # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}

            # Exporting the loops for washU visualization
                final_loop_export <- function(input_df, disease_state, date){
                    write.table(input_df, paste0("hic_contacts_", "all_", disease_state, "_non_merged_q", filter_qval, "_", date, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
                }

                
        #

    #

    # Directories for importing and exporting the data 
        hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"   
        bed_export_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr/merge_beds"   

        qval <- 0.01
        filter_qval <- 0.005
        contact_thresh <- 3
        date <- "30SEP2020"


    # Start of the script

        # Loading and pre-processing of data

            # FitHiC data
                setwd(hic_data_directory)
                hic_data_a <- read.table(paste("FitHiC_all_qval", qval, ".txt", sep=""), header=TRUE)       
                hic_data <- hic_data_a %>% filter(q.value < filter_qval)
               
                disease_hic_data <- hic_data %>% filter(d.contactCount > contact_thresh) %>% select(c(1:6))
                healthy_hic_data <- hic_data %>% filter(h.contactCount > contact_thresh) %>% select(c(1:6))

              setwd(hic_data_directory)
                final_loop_export(disease_hic_data, "disease", date)
                final_loop_export(healthy_hic_data, "healthy", date)


  # End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
    # From there you can import them into a .json for washU visualization
    # conda activate samtools
        awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_all_healthy_merged_q0.005_30SEP2020.txt | sort -k1,1 -k2,2n > hic_contacts_all_healthy_merged_q0.005_30SEP2020_out_washU.bed
        bgzip hic_contacts_all_healthy_merged_q0.005_30SEP2020_out_washU.bed
        tabix -p bed hic_contacts_all_healthy_merged_q0.005_30SEP2020_out_washU.bed.gz


    # Copy the data to the BioAdHoc folder
        cp *_30SEP2020* /mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/washU_filtered/

