
################# DATA ANALYSIS TEST ############

        library('tidyverse', lib="/home/jcayford/r_libs")
        library('viridisLite')
      
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



                # Function to get the mid-points of the raw data and change the column names
                    import_cleanup <- function(raw_data_frame){
                                    df <- data.frame(raw_data_frame[,1]+2500, raw_data_frame[,2]+2500, raw_data_frame[,3])
                                    colnames(df) <- c("fragmentMid1", "fragmentMid2", "contactCount")
                                    df
                    }




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

      
      
        
        # Setting the working directories and other factors
                main_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/fithic/contacts_outputs"
                raw_data_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/juicer_files"
                export_wd <- "/home/jcayford/HCM/Analysis/HiC_R"
                q_value <- 0.01
                chr_num <- 1


        # Importing the data from Fit-HiC
                setwd(main_wd)
                fithic_h_df <- read.table("fithic_healthy_chr1_KR.txt.spline_pass1.res5000.significances.txt", header=T)
                fithic_d_df <- read.table("fithic_disease_chr1_KR.txt.spline_pass1.res5000.significances.txt", header=T)
                setwd(raw_data_wd)

        # Importing the raw data
                d_raw <- read.table("disease_chr1_KR.txt");d_raw_df <- import_cleanup(d_raw);gc()
                h_raw <- read.table("healthy_chr1_KR.txt");h_raw_df <- import_cleanup(h_raw);gc()
             
        # Generating the combined file between healthy/disease and determining the count difference          
                h_final <- combining_func("healthy", q_value)
                d_final <- combining_func("disease", q_value)
                df_final <- rbind(h_final, d_final) 

        # Exporting the data

                df_final2 <- df_final %>% distinct(fragmentMid1, fragmentMid2, .keep_all=TRUE)
                
                setwd(export_wd)
                write.table(df_final2, paste("FitHiC_chr", chr_num, "_qval", q_value, ".txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)




        ggplot(df_final2, aes(x=-log2(h.contactCount), y=-log2(d.contactCount), col=sub)) +
                geom_point(size=0.75) +
                theme +
                scale_x_continuous(name="-log2 Healthy Contact Counts") +
                scale_y_continuous(name="-log2 Disease Contact Counts") +
                ggtitle("Disease vs Healthy Counts | chr1 | q=0.01") +
                labs(color="Count Sub", size=0.5) + 
            
                scale_colour_gradientn(colors = viridis(256, option = "D", direction = -1))





#### SHOULD PROBABLY OUTPUT THESE FILES AND OPEN ANOTHER SCRIPT












