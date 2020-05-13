
################# DATA ANALYSIS TEST ############

        library('tidyverse', lib="/home/jcayford/r_libs")
        library('viridisLite')
      
        # Functions
                # Quick head() with the amount of rows
                    w <- function(x){print(head(x));nrow(x)}
                # Rounding data
                    round_to <- function(x, to = window_length) round(x/to)*to

                # Function to combine the healthy and disease contact counts based on a given q_value filter
                    combining_func <- function(state, q_value){
                                    if(state == "healthy"){
                                            df_filter <- fithic_h_df
                                            df_raw_a <- d_raw_df 
                                            df_raw_b <- c_raw_df
                                            columns <- c(1:4, 6, 7, 5, 10)
                                    } 
                                    if(state == "disease"){
                                            df_filter <- fithic_d_df
                                            df_raw_a <- h_raw_df
                                            df_raw_b <- c_raw_df
                                            columns <- c(1:4, 6, 7, 10, 5)         
                                    }  

                                    if(state == "ctcf_ko"){
                                            df_filter <- fithic_c_df
                                            df_raw_a <- h_raw_df
                                            df_raw_b <- d_raw_df
                                            columns <- c(1:4, 6, 7, 10)         
                                    }     


                            q_filter <- df_filter %>% filter(q.value < q_value)
                            combined_a <- left_join(q_filter, df_raw_a, by=c("fragmentMid1", "fragmentMid2")) 
                            combined_b <- left_join(q_filter, df_raw_b, by=c("fragmentMid1", "fragmentMid2")) 
                            
                            if(state == "disease" | state == "healthy"){combined_df <- data.frame(combined_a[, columns], combined_b[,10])}
                            if(state == "ctcf_ko"){combined_df <- data.frame(combined_a[, columns], combined_b[,10], combined_a[, 5])}                             
                            
                            combined_df[is.na(combined_df)] <- 0
                            colnames(combined_df) <- c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "p.value", "q.value", "h.contactCount", "d.contactCount", "c.contactCount") 
                            combined_df
                    }



                # Function to get the mid-points of the raw data and change the column names
                    import_cleanup <- function(raw_data_frame){
                                    df <- data.frame(raw_data_frame[,1]+2500, raw_data_frame[,2]+2500, round_to(raw_data_frame[,3], to=1))
                                    colnames(df) <- c("fragmentMid1", "fragmentMid2", "contactCount")
                                    df
                    }
        #


      
      
        
        # Setting the working directories and other factors
                main_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/fithic/contacts_outputs"
                raw_data_wd <- "/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/juicer_files"
                export_wd <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
                q_value <- 0.01

    # For loop to run all the chromosomes
        all_chromosomes <- NULL
        for(i in 2:19){
            print(paste(Sys.time(), "............ working on chromosome: ", i, sep=""))
            chr_num <- i

            # Importing the data from Fit-HiC
                    setwd(main_wd)                   
                    fithic_h_df <- read.table(paste0("fithic_healthy_chr", chr_num, "_KR.txt.spline_pass1.res5000.significances.txt"), header=TRUE)
                    fithic_d_df <- read.table(paste0("fithic_disease_chr", chr_num, "_KR.txt.spline_pass1.res5000.significances.txt"), header=TRUE)
                    fithic_c_df <- read.table(paste0("fithic_ctcf_ko_chr", chr_num, "_KR.txt.spline_pass1.res5000.significances.txt"), header=TRUE)
                    setwd(raw_data_wd)

            # Importing the raw data
                    h_raw <- read.table(paste0("healthy_chr", chr_num, "_KR.txt"));h_raw_df <- import_cleanup(h_raw);gc()
                    d_raw <- read.table(paste0("disease_chr", chr_num, "_KR.txt"));d_raw_df <- import_cleanup(d_raw);gc()
                    c_raw <- read.table(paste0("ctcf_ko_chr", chr_num, "_KR.txt"));c_raw_df <- import_cleanup(c_raw);gc()
                              
            # Generating the combined file between healthy/disease and determining the count difference          
                    h_final <- combining_func("healthy", q_value)
                    d_final <- combining_func("disease", q_value)
                    c_final <- combining_func("ctcf_ko", q_value)
                    df_final <- rbind(h_final, d_final, c_final) 

            # Exporting the data
                    df_final2 <- df_final %>% distinct(fragmentMid1, fragmentMid2, .keep_all=TRUE) 
                    all_chromosomes <- rbind(all_chromosomes, df_final2)
        }

                    setwd(export_wd)
                    write.table(all_chromosomes, paste0("FitHiC_all_qval", q_value, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)




        # # GGplot settings
        # # Theme for general ggplot
        #          theme<-theme_classic() +
        #             theme(
        #                 text = element_text(color = "grey20"),
        #                 plot.title = element_text(hjust = 0.5),
        #                 axis.text.x = element_text(size=rel(1.5)), 
        #                 axis.text.y = element_text(size=rel(1.5), angle=90, hjust=0.5)
        #             )
        # #
        # ggplot(df_final2, aes(x=-log2(h.contactCount), y=-log2(d.contactCount), col=sub)) +
        #         geom_point(size=0.75) +
        #         theme +
        #         scale_x_continuous(name="-log2 Healthy Contact Counts") +
        #         scale_y_continuous(name="-log2 Disease Contact Counts") +
        #         ggtitle("Disease vs Healthy Counts | chr1 | q=0.01") +
        #         labs(color="Count Sub", size=0.5) + 
            
        #         scale_colour_gradientn(colors = viridis(256, option = "D", direction = -1))
