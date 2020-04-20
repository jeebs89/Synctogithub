
## Chr 1 example to generate .txt from the .hic files

# java -jar /home/jcayford/tools/juicer/PBS/scripts/juicer_tools.jar dump observed KR /home/jcayford/HCM/data/GSM2544836_No_Tx.hic 1 1 BP 5000 healthy_chr1_KR.txt



# Generate all the names in R before cat-ing them together.  
                R
                awk_strings <- function(file, chr, dis_hel){
                print(paste(paste("java -jar /home/jcayford/tools/juicer/PBS/scripts/juicer_tools.jar dump observed KR /home/jcayford/HCM/data/", file, " ",chr, " ", chr, " BP 5000 ", dis_hel, "_chr", i, "_KR.txt", sep="")))
                }


                h_strings <- c()
                d_strings <- c()
                for( i in 1:19){
                        h_strings <- rbind(h_strings, awk_strings("GSM2544836_No_Tx.hic", i, "healthy"))
                        d_strings <- rbind(d_strings, awk_strings("GSM2544839_Tx_Cre-plus.hic", i, "disease"))
                }

                strings <- rbind(h_strings, d_strings)
                write.table(strings, "juicer_inputs.txt", col.names=F, row.names=F, quote=F)
                q()

#



# Run all of the names through on the cluster to generate 19 individual .txt files

                R
                library('tidyverse', lib="/home/jcayford/r_libs")

                        df<-read.table("healthy_chr1_NONE.txt")
                        mids <- df[,1]+2500

                        # Generation of the interactions data.frame
                                interactions<-data.frame("chr1", mids, "chr1", (df[,2]+2500), df[,3])
                                colnames(interactions)<-c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "contactCount")
                                interactions[is.na(interactions)]=0              
                                write.table(interactions, "healthy_interactions_chr1_NONE.txt", col.names=F, row.names=F, quote=F)

                        # Generation of the fragments data.frame
                                fragments_pre <- data.frame(df[,1]+2500, df[,3])
                                colnames(fragments_pre) <- c("fragmentMid", "ContactCount")
                                fragments_mid <- fragments_pre %>% group_by(fragmentMid) %>% summarise(ContactCount=sum(ContactCount))
                                fragments <- data.frame("chr"="chr1", "extra"=1, fragments_mid, "extra.1"=1)
                                fragments[is.na(fragments)]=0          
                                write.table(fragments, "healthy_fragments_chr1_NONE.txt", col.names=F, row.names=F, quote=F)

                        q()

#     


# Zipping and moving the files to the correct place for Fit-HiC
        gzip -ck healthy_interactions_chr1_NONE.txt > healthy_interactions_chr1_NONE.gz
        mv healthy_interactions_chr1_NONE.gz fithic_test

        gzip -ck healthy_fragments_chr1_NONE.txt > healthy_fragments_chr1_NONE.gz
        mv healthy_fragments_chr1_NONE.gz fithic_test
#



## Need to set lower bound to about 10-20 kb and the upper bound to 1-2

# Running the Fit-Hic code
        fithic -r 5000 -l test_fithic_healthy_chr1 -f healthy_fragments_chr1_NONE.gz -i healthy_interactions_chr1_NONE.gz -o test_fithic_healthy_chr1 -U 25000 -x intraOnly
#












################# DATA ANALYSIS TEST ############



R
        library('tidyverse', lib="/home/jcayford/r_libs")
        library('viridisLite')
      
        # Functions
                # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}


        # Setting the working directories
                main_wd <- "/home/jcayford/HCM/Juicer/fithic_test"
                raw_data_wd <- "/home/jcayford/HCM/Juicer/"
                fithic_h_wd <- "/home/jcayford/HCM/Juicer/fithic_test/test_fithic_healthy_chr1"
                fithic_d_wd <- "/home/jcayford/HCM/Juicer/fithic_test/test_fithic_disease_chr1"

        # Importing the data from Fit-HiC
                setwd(fithic_h_wd)
                fithic_h_df <- read.table("test_fithic_healthy_chr1.spline_pass1.res5000.significances.txt", header=T)

                setwd(fithic_d_wd)
                fithic_d_df <- read.table("test_fithic_disease_chr1.spline_pass1.res5000.significances.txt", header=T)
                setwd(raw_data_wd)


        # Joining and getting rid of all the windows which are not found in both comparisons
                join_data <- full_join(fithic_h_df, fithic_d_df, by=c("fragmentMid1", "fragmentMid2"))
                join_data_b <- data.frame(join_data[,c(1, 2, 4:7, 12:14)])
                colnames(join_data_b) <- c("chr", "fragmentMid1", "fragmentMid2", "h.contact", "h.p.value", "h.q.value", "d.contact", "d.p.value", "d.q.value")             
                sub <- na.omit(data.frame(join_data_b, "sub"=(join_data_b$h.contact - join_data_b$d.contact)))
                
        # Ranking and getting rid of one large outliar
                ranked <- sub %>% arrange(desc(sub))
                ranked_b <- ranked %>% filter(sub<quantile(sub, 0.01))
                ranked_c <- ranked %>% filter(sub>quantile(sub, 0.99))
              


                ggplot(ranked_b, aes(x=(h.contact), y=(d.contact), col=sub)) +
                        geom_point() +
                        geom_point(data=ranked_c[2:nrow(ranked_c),], aes(x=(h.contact), y=(d.contact), col=sub)) +
                        theme + 
                        scale_colour_gradientn(colors = viridis(256, option = "D", direction = -1))



        # GGplot settings
        # Theme for general ggplot
                 theme<-theme_classic() +
                    theme(
                        text = element_text(color = "grey20"),
                        plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(size=rel(1.5)), 
                        axis.text.y = element_text(size=rel(1.5), angle=90, hjust=0.5)
                    )




                ggplot(ranked_b, aes(x=-log2(h.contact), y=-log2(d.contact), col=sub)) +
                        geom_point(size=1) +
                        theme + 
                        scale_colour_gradientn(colors = viridis(256, option = "D", direction = 1))




                ggplot(ranked_b, aes(x=(h.contact), y=(d.contact), col=sub)) +
                        geom_point(size=1) +
                        theme + 
                        scale_colour_gradientn(colors = viridis(256, option = "D", direction = -1))

                        








