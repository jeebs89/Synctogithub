# Mustache

library('tidyverse', lib="/home/jcayford/r_libs")


# Funcations
    # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}

    # Function to filter out the regions which are 1 bin away from each other
        # df_a and df_b are the inputs of healthy or disease which have already been filtered for the chromosome
            overlaps_func<- function(df_a, df_b, row_line){     
                    df_subs <- data.frame(df_a, "bin1"=as.numeric(df_a[,2]) - as.numeric(df_b[row_line,2]), "bin2"=as.numeric(df_a[,5]) - as.numeric(df_b[row_line,5]))          
                    df_subs %>% filter(bin1 == 0 | bin1 == 5000 | bin1== -5000, bin2 == 0 | bin2 == 5000 | bin2 == -5000)
            }

    # Function to rebin the data.frame for the overlaps. 
        # Input_df is the chromosome overlaps      
            rebin_h_to_d <- function(input_df){
                df_a <- input_df %>% filter(bin1 >= 0 | bin1 <= 0 | bin2 >= 0 | bin2 <= 0 )
                df_b <- data.frame(df_a[,1], df_a[,2]-df_a[,10], df_a[,3],df_a[,4], df_a[,5]-df_a[,11], df_a[,6:8],  "condition"="both", df_a[,10:11])
                colnames(df_b) <- c("BIN1_CHR","BIN1_START","BIN1_END","BIN2_CHROMOSOME","BIN2_START","BIN2_END","FDR","DETECTION_SCALE","condition", "bin1", "bin2")
                df_b
            }
    
#

# Setting the working directories 
    main_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files"
    disease_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files/disease/"
    healthy_wd <- "/home/jcayford/HCM/Analysis/Mustache/hic_files/healthy/"

#
# q.value threshold, since the input files have a threshold of 0.1
qval_threshold <- 0.01

#
# Importing the files and filtering based on the given threshold for the qval
    setwd(disease_wd)
        disease <- read.table("mustache_disease_all_out_0.1.tsv", header=TRUE, stringsAsFactors=F)
        disease_q_filtered <- disease %>% filter(FDR < qval_threshold) %>% arrange(FDR)
    setwd(healthy_wd)
        healthy <- read.table("mustache_healthy_all_out_0.1.tsv", header=TRUE, stringsAsFactors=F)
        healthy_q_filtered <- healthy %>% filter(FDR < qval_threshold) %>% arrange(FDR)


        min_rows <- min(nrow(disease_q_filtered), nrow(healthy_q_filtered))
        disease_filtered_a <- disease_q_filtered[1:min_rows, ]
        disease_filtered <- data.frame(disease_filtered_a, "condition"="disease")
        
        healthy_filtered_a <- healthy_q_filtered[1:min_rows, ]
        healthy_filtered <- data.frame(healthy_filtered_a, "condition"="healthy")

#
# Generation of the new healthy dataframe which is binned to match with the healthy windows with the disease
    new_healthy <- NULL
    for(j in 1:19){
        print(paste(Sys.time(), "...............starting Chromosome: ", j, sep=""))
            chrom <- paste("chr", j, sep="")
                # Filtering the files based on their chromosome in order to reduce the filtering time
                    disease_chr <- disease_filtered %>% filter(BIN1_CHR == chrom)            
                    healthy_chr <- healthy_filtered %>% filter(BIN1_CHR == chrom)
                # For loop to get the overlaps between the healthy and the disease
                    healthy_chr_overlaps <- NULL
                    for(i in 1:nrow(healthy_chr)){
                        healthy_chr_overlaps <- rbind(healthy_chr_overlaps, overlaps_func(healthy_chr, disease_chr, i))
                    }
                # Re-binning the healthy data to generate a dataset which will overlap with the disease data
                # If the bin is not an overlap within 1 bin, the bin will be NA
                    healthy_chr_b <- data.frame(healthy_chr, "bin1"=NA, "bin2"=NA)
                    new_healthy_a <- rebin_h_to_d(healthy_chr_overlaps)
                # Combining the overlaps and non-overlaps data of the healthy
                    comined_a <- rbind(new_healthy_a, healthy_chr_b)
                    new_healthy_combined <- comined_a %>% distinct(BIN1_END, BIN2_END, FDR, DETECTION_SCALE, .keep_all=TRUE)
                # Generation of the new healthy data.frame
                    new_healthy <- rbind(new_healthy, new_healthy_combined)
    }

#
# Finding the overlaps with the newly binned healthy data.frame and putting it together
    # New healthy is a data.frame which should line up with the diseased ones.
        new_disease <- data.frame(disease_filtered, "bin1"=NA, "bin2"=NA)
    # Overlaps between the healthy and disease
        overlaps <- new_healthy %>% filter(condition == "both")
        disease_overlaps_a <- rbind(overlaps, new_disease)
        disease_ovalaps_b <- disease_overlaps_a %>% distinct(BIN1_CHR, BIN1_START, BIN2_START, .keep_all=TRUE)
        distinct_healthy <- new_healthy %>% filter(condition == "healthy")
    # Putting together the healthy, disease, and both overlaps together and arranging the chromosomes
        all_overlaps <- rbind(disease_ovalaps_b, distinct_healthy)
        final_overlaps <- all_overlaps %>% arrange(BIN1_CHR, BIN1_START, BIN1_END, BIN2_CHROMOSOME, BIN2_START, BIN2_END)
    # Exporting the final data.frame
        setwd(main_wd)
        write.table(final_overlaps, paste("mustache_hic_filtered", "_qval", qval_threshold, ".txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)







