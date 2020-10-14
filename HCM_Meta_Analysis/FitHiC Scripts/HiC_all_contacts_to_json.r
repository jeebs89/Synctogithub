
# Generating HiC for json... for disease and healthy... all loops included




library('tidyverse', lib="/home/jcayford/r_libs")
    
       # Functions
            # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}

            # Exporting the loops for washU visualization
                final_loop_export <- function(input_df, disease_state, date){
                    write.table(input_df, paste0("hic_contacts_", "all_", disease_state, "_non_merged_q", filter_qval_name, "_", date, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
                }

                
        #

    #

    # Directories for importing and exporting the data 
        hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"   
        bed_export_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr/merge_beds"   

        qval <- 0.01
        filter_qval <- 0.00005
        filter_qval_name <- "000005"
        contact_thresh <- 3
        today <- Sys.Date()
        date <- format(today, format="%d_%b_%Y")


    # Start of the script

        # Loading and pre-processing of data

            # FitHiC data
                setwd(hic_data_directory)
                hic_data <- read.table(paste("FitHiC_all_qval", qval, ".txt", sep=""), header=TRUE)       

                disease_hic_data <- hic_data %>% filter(d.contactCount > contact_thresh & q.value < filter_qval) %>% select(c(1:6))
                healthy_hic_data <- hic_data %>% filter(h.contactCount > contact_thresh & q.value < filter_qval) %>% select(c(1:6))

              setwd(hic_data_directory)
                final_loop_export(disease_hic_data, "disease", date)
                final_loop_export(healthy_hic_data, "healthy", date)


  # End of the script. Go to the export_directory and run the awk script and then use Samtools to bgzip and tabix.
    # From there you can import them into a .json for washU visualization
    # conda activate samtools
        awk '{if (NR > 1) {if ($NF > 0) {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"($2-1)"\t"($2+1)"\t"$3":"($4-1)"-"($4+1)",500\t"(NR-1)"\t."}}}' hic_contacts_all_healthy_non_merged_q000005_top50percent_07_Oct_2020.txt | sort -k1,1 -k2,2n > hic_contacts_all_healthy_non_merged_q000005_top50percent_07_Oct_2020_out_washU.bed
        bgzip hic_contacts_all_healthy_non_merged_q000005_top50percent_07_Oct_2020_out_washU.bed
        tabix -p bed hic_contacts_all_healthy_non_merged_q000005_top50percent_07_Oct_2020_out_washU.bed.gz


    # Copy the data to the BioAdHoc folder
        cp *percent_07_Oct_2020* /mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/washU_filtered/






test <- data.frame(disease_hic_data, "log"=-log(disease_hic_data$q.value))
disease2 <- test %>% filter(log > quantile(test$log, 0.5)) %>% select(c(1:6))


test <- data.frame(healthy_hic_data, "log"=-log(healthy_hic_data$q.value))
healthy2 <- test %>% filter(log > quantile(test$log, 0.5)) %>% select(c(1:6))


filter_qval_name <- "000005_top50percent"

 setwd(hic_data_directory)
                final_loop_export(disease2, "disease", date)
                final_loop_export(healthy2, "healthy", date)

