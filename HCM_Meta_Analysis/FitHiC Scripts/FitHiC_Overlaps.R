# FitHiC Overlaps. There are two outputs here. The first is a .txt file of the contacts found in FitHiC.
# The plot is the -log2 counts from the HiC file. This is completed for each chromosome and will need to be
# Cat together for the whole genome data. 



library('tidyverse', lib="/home/jcayford/r_libs")
library('viridisLite', lib="/home/jcayford/r_libs")


        # Functions
                # Quick head() with the amount of rows
                    w <- function(x){print(head(x));nrow(x)}

            
                # Function to combine the healthy and disease contact counts based on a given q_value filter
                    combining_func <- function(state, q_value){
                                    if(state == "healthy"){
                                            df_filter <- fithic_h_df
                                            df_raw <- d_raw_df 
                                            columns <- c(1:4, 6, 7, 5, 10)
                                    } 
                                    if(state == "disease"){
                                            df_filter <- fithic_d_df
                                            df_raw <- h_raw_df
                                            columns <- c(1:4, 6, 7, 10, 5)         
                                    }     

                            q_filter <- df_filter %>% filter(q.value < q_value)
                            combined_a <- left_join(q_filter, df_raw, by=c("fragmentMid1", "fragmentMid2"))
                            combined_a[is.na(combined_a)] <- 0
                            combined_df <- combined_a[, columns]
                            colnames(combined_df) <- c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "p.value", "q.value", "h.contactCount", "d.contactCount") 

                            data.frame(combined_df, "sub"=combined_df[,7] - combined_df[,8])
                    }

                # Function to import the fithic data and to rearrange the columns
                    fithic_import <- function(state){
                                    data <- paste(state, "_chr", chr_num,".spline_pass1.res5000.significances.txt" , sep="")
                                    
                                    if(state == "healthy"){setwd(healthy_wd)}
                                    if(state == "disease"){setwd(disease_wd)}
                                                            
                                    df <- read.table(data, header=TRUE)
                                    df_2 <- df %>% select(c(1:4, 6, 7, 5))  
                                    df_2 %>% filter(q.value < q_value)
                                }




                # Function to get the mid-points of the raw data and change the column names
                        raw_import_cleanup <- function(state){
                                setwd(raw_data_wd)
                                data <- paste(state, "_chr", chr_num,"_NONE.txt" , sep="")
                                df <- read.table(data)
                                df_1 <- data.frame(df[,1]+2500, df[,2]+2500, df[,3])
                                colnames(df_1) <- c("fragmentMid1", "fragmentMid2", "contactCount");gc();
                                df_1
                        }


                # Function to combined the raw counts of the other disease state
                        contact_combining <- function(state){
                                        if(state=="healthy"){
                                                df <- fithic_h
                                                df2 <- d_raw
                                                col_names <- c("h_contactCount", "d_contactCount")
                                        } else{
                                                df <- fithic_d
                                                df2 <- h_raw
                                                col_names <- c("d_contactCount", "h_contactCount")
                                        }
                                combined <- left_join(df, df2, by=c("fragmentMid1", "fragmentMid2"))
                                colnames(combined)[c(7,8)] <- col_names
                                data.frame(combined, "condition"=state)
                        
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
      
        
        # Setting the working directories and q_value threshold
                                
                main_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw"
                healthy_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/healthy/fithic_10k_2mb"
                disease_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/disease/fithic_10k_2mb"
                raw_data_wd <- "/home/jcayford/HCM/Analysis/Juicer/indiv_chr/none_normalized"
                export_wd <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes"
                plot_wd <- "/home/jcayford/HCM/Analysis/HiC_R/plots"
                
                q_value <- 0.1             


        # For loop to complete each chromosome in the analysis
            for(i in 1:19){
                    chr_num <- i
                    print(paste(Sys.time(), "..........starting Chromosome: ", chr_num, sep=""))
                            # Importing the data from Fit-HiC
                                    fithic_h <- fithic_import("healthy")
                                    fithic_d <- fithic_import("disease")           
                                    
                            # Importing the raw data
                                    d_raw <- raw_import_cleanup("disease")
                                    h_raw <- raw_import_cleanup("healthy")

                            # Generating the combined file between healthy/disease and determining the count difference          
                                    fithic_h <- contact_combining("healthy")
                                    fithic_d <- contact_combining("disease")

                            # combining the two data.frames and looking for where both windows were found with FitHiC and exporting the data                           
                                    print(paste(Sys.time(), "..........combining files",sep=""))
                                    final_a <- rbind(fithic_h, fithic_d)
                                    final_uniq <- final_a %>% distinct(fragmentMid1, fragmentMid2, .keep_all=TRUE)
                                    input_df <- data.frame(final_a[duplicated(final_a[,c(2, 4)]),1:8], "condition"="both")                                                                                            
                                # Finding the duplicates and removing them from both the healthy and diseased samples. This could
                                # likley be done in a much cleaner way... 
                                    df_healthy_a <- rbind(fithic_h, input_df)
                                    df_healthy_b <- df_healthy_a %>% distinct(fragmentMid1, fragmentMid2, .keep_all=TRUE)
                                    df_disease_a <- rbind(fithic_d, input_df)
                                    df_disease_b <- df_disease_a %>% distinct(fragmentMid1, fragmentMid2, .keep_all=TRUE)

                                    df_final <- rbind(df_healthy_b, df_disease_b, input_df)                                                       
                                    df_final2 <- df_final %>% arrange(fragmentMid1, fragmentMid2)
            
                                print(paste(Sys.time(), "..........exporting files for chromosome: ", chr_num,sep=""))

                                    setwd(export_wd)
                                    write.table(df_final2, paste("FitHiC_chr", chr_num, "_qval", q_value, ".txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)
                            
                            # Saving a plot fo the counts
                            setwd(plot_wd)
                            ggplot(df_final2, aes(x=-log2(h_contactCount), y=-log2(d_contactCount), col=(h_contactCount - d_contactCount))) +
                                    geom_point(size=0.75) +
                                    theme +
                                    scale_x_continuous(name="-log2 Healthy Contact Counts") +
                                    scale_y_continuous(name="-log2 Disease Contact Counts") +
                                    ggtitle(paste("Disease vs Healthy Counts | chr", chr_num, " | q=", q_value, sep="")) +
                                    labs(color="Count Sub", size=0.5) +                       
                                    scale_colour_gradientn(colors = viridis(256, option = "D", direction = -1))
                            ggsave(paste("FitHiC_chr", chr_num, "_qval", q_value, ".png", sep=""))
                            gc();
            }
