## Uses the import file from the FitHiC_Overlaps.R for the hic data

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
        # Essentially the same filter for the RNA data compared to the filtered HiC data
            rna_filter <- function(distance){
                        df <- data.frame(rna_data[,c(1,2)], (rna_means+(distance)), rna_data[,c(5:ncol(rna_data))]) 
                        
                        colnames(df)[c(1,3)] <- c("chr1", "fragmentMid1")
                        df_uniq_a <- df %>% distinct(type, fragmentMid1, .keep_all=TRUE)
                        a <- inner_join(hic_filtered, df_uniq_a, by=c("chr1", "fragmentMid1"))
                        a1 <- data.frame(a, "rna_contactMid"=1)

                        colnames(df)[c(1,3)] <- c("chr1", "fragmentMid2")
                        df_uniq_b <- df %>% distinct(type, fragmentMid2, .keep_all=TRUE)
                        b <- inner_join(hic_filtered, df_uniq_b, by=c("chr1", "fragmentMid2"))
                        b1 <- data.frame(b, "rna_contactMid"=2)
                        
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
            rna_data_directory <- "/home/jcayford/HCM/Analysis/RNA"
            hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R"
            chip_beds_directory <- "/home/jcayford/HCM/Analysis/ChIP_files/"
            
            export_directory <- "/home/jcayford/HCM/Analysis/Gene_lists"


        chip_bed <- "polr2a.bed" 
        #chip_bed <- "mef2a.bed"
        tf <- "pol2"  
        qval <- 0.01


    # Start of the script

        # Loading and pre-processing of data
            # RNA
                setwd(rna_data_directory)
                rna_data_a <- read.table("rna_analysis_all_windows.txt", header=TRUE)
                rna_data <- data.frame(paste("chr", rna_data_a[,1], sep=""), rna_data_a[,2:ncol(rna_data_a)])
                colnames(rna_data) <- colnames(rna_data_a)

            # ChIP Peaks
                setwd(chip_beds_directory)
                chip_peaks_a <- read.table(chip_bed)
                        if(chip_peaks_a[1,1]==1){chip_peaks <- data.frame(paste("chr", chip_peaks_a[,1], sep=""), chip_peaks_a[,2:5])}
                        if(chip_peaks_a[1,1]=="chr1"){chip_peaks <- chip_peaks_a[,1:5]}                       
                colnames(chip_peaks) <- c("chr1", "start", "end", "peak", "height")     

            # FitHiC data
                setwd(hic_data_directory)
                hic_data_a <- read.table(paste("FitHiC_qval", qval, "_all.txt", sep=""), header=TRUE)
                hic_data <- data.frame(hic_data_a[,1], df_c(hic_data_a, 2), hic_data_a[,3],
                                    df_c(hic_data_a, 4), df_c(hic_data_a, 5), df_c(hic_data_a, 6),
                                    df_c(hic_data_a, 7), df_c(hic_data_a, 8), df_c(hic_data_a, 9))               
                colnames(hic_data) <- colnames(hic_data_a)           


        # Filtering FitHiC data with the ChIP-Seq Peaks

            chip_means <- round_to(rowMeans(chip_peaks[,2:3]), to=2500)
            rna_means <- round_to(rowMeans(rna_data[,3:4]), to=2500)      
       
            # Filtering the HiC data
                hic_filtered_a <- NULL
                for(i in c((-5000), (-2500), 0, 2500, 5000)){
                    hic_filtered_a <- rbind(hic_filtered_a, chip_filter(i))
                }                
                hic_filtered <- hic_filtered_a %>% distinct(fragmentMid1, fragmentMid2,  .keep_all=TRUE)

            # Filtering the RNA data
                rna_filtered_a <- NULL
                for(i in c((-10000), (-5000), (-2500), 0, 2500, 5000, 10000)){
                    rna_filtered_a <- rbind(rna_filtered_a, rna_filter(i))
                }                
                rna_filtered <- rna_filtered_a %>% distinct(gene, trans_exon, .keep_all=TRUE)
                genes <- rna_filtered %>% filter(pval < 0.05) %>% arrange(pval)

        ggplot(rna_filtered, aes(x=-log2(h_mean), y=-log2(d_mean), col=pval)) +
                geom_point(size=2, shape=16) +
                theme +
                scale_x_continuous(name="-log2 Healthy TPM mean") +
                scale_y_continuous(name="-log2 Disease TPM mean") +
                ggtitle(paste(tf, " | Disease vs Healthy TPMs | FitHiC q=", qval, sep="")) +          
                scale_colour_gradientn(colors = viridis(256, option = "B", direction = 1))


            setwd(export_directory)
            write.table(rna_filtered, paste("rna_analysis_", tf, "_qval", qval, "_.txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)
            write.table(genes, paste("gene_list_", tf, "_qval", qval, "_pval_05.txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)

# END