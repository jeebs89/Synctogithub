# Generation of the import files for FitHiC, they still need to be zipped after this script is run.

library('tidyverse', lib="/home/jcayford/r_libs")


# Functions
        # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}

        # Function to generate the necessary interactions and fragments files for input into FitHiC
                # Chr_num is coded to be a chr number in a for loop but it can be done separately
                # Status is either "healthy" or "disease". No other inputs will work
                        frags_inters_func <- function(chr_num, status, bias){
                        # Clearing variables for the loop
                                interactions <- NULL
                                fragments <- NULL
                                normalization_caps <- toupper(normalization)                       
                        # Setting up the directories, the chromosome, and the file to get
                                chromosome <- paste("chr", chr_num, sep="")
                                export_dir_final <- paste(export_dir, "/", status, sep="")
                                df_file <- paste(status, "_", chromosome, "_",normalization_caps, ".txt", sep="")
                        # Reading the file and setting the export directory based on whether the sample is healty or disease
                                setwd(import_dir)
                                df <- read.table(df_file)
                                mids <- df[,1]+2500
                                setwd(export_dir_final)
                        # Generation of the interactions data.frame
                                interactions<-data.frame(chromosome, mids, chromosome, (df[,2]+2500), df[,3])
                                colnames(interactions)<-c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "contactCount")
                                interactions[is.na(interactions)]=0     
                                write.table(interactions, paste(status, "_interactions_", chromosome, "_", normalization_caps, ".txt", sep=""), col.names=F, row.names=F, quote=F)         
                        # Generation of the fragments data.frame
                                fragments_pre <- data.frame(df[,1]+2500, df[,3])
                                colnames(fragments_pre) <- c("fragmentMid", "ContactCount")
                                fragments_mid <- fragments_pre %>% group_by(fragmentMid) %>% summarise(ContactCount=sum(ContactCount))
                                fragments <- data.frame("chr1"=chromosome, "extra"=1, fragments_mid, "extra.1"=1)
                                fragments[is.na(fragments)]=0          
                                write.table(fragments, paste0(status, "_fragments_", chromosome, "_", normalization_caps, ".txt"), col.names=F, row.names=F, quote=F)
                        # Generation of a bias file, need to set the value to True
                                if(bias=="TRUE"){
                                        print(paste0(Sys.time(), "........working on the bias file"))
                                        setwd(import_dir)
                                        df_bias <- read.table(paste0(status, "_", chromosome, "_KR.txt"))
                                        mids_bias <- df_bias[,1]+2500
                                        
                                        setwd(export_dir_final)
                                        fragments_bias_pre <- data.frame(df_bias[,1]+2500, df_bias[,3])
                                        colnames(fragments_bias_pre) <- c("fragmentMid", "ContactCount")
                                        fragments_bias_mid <- fragments_bias_pre %>% group_by(fragmentMid) %>% summarise(ContactCount=sum(ContactCount))
                                        bias_df <- data.frame(chromosome, fragments_bias_mid)
                                        colnames(bias_df) <- c("chr1", "midpoint", "bias")
                                        bias_df[is.na(bias_df)]=0
                              
                                        write.table(bias_df, paste0(status, "_bias_", chromosome, "_KR.txt"), col.names=F, row.names=F, quote=F)
                                }
                        }

#
          
# Normalization. Generally will be either kr or none.
normalization <- "none"

# Setting the import and export directories

import_dir <- paste0("/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/juicer_files")
export_dir <- paste0("/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/")

for(j in 1:19){
    print(paste(Sys.time(), "............ working on healthy chromosome: ", j, sep=""))
        frags_inters_func(j, "healthy", TRUE)
    print(paste(Sys.time(), "............ working on disease chromosome: ", j, sep=""))
        frags_inters_func(j, "disease", TRUE)
    print(paste(Sys.time(), "............ working on ctcf_ko chromosome: ", j, sep=""))
        frags_inters_func(j, "ctcf_ko", TRUE)

}



for(j in 1:19){
    print(paste(Sys.time(), "............ working on healthy chromosome: ", j, sep=""))
        frags_inters_func(j, "healthy", TRUE)
    print(paste(Sys.time(), "............ working on disease chromosome: ", j, sep=""))
        frags_inters_func(j, "disease", TRUE)
    print(paste(Sys.time(), "............ working on ctcf_ko chromosome: ", j, sep=""))
        frags_inters_func(j, "ctcf_ko", TRUE)

}


for(j in 1:19){
    print(paste(Sys.time(), "............ working on healthy chromosome: ", j, sep=""))
        frags_inters_func(j, "healthy", TRUE)
    print(paste(Sys.time(), "............ working on disease chromosome: ", j, sep=""))
        frags_inters_func(j, "disease", TRUE)
    print(paste(Sys.time(), "............ working on ctcf_ko chromosome: ", j, sep=""))
        frags_inters_func(j, "ctcf_ko", TRUE)

}
