# Generation of the import files for FitHiC, they still need to be zipped after this script is run.

library('tidyverse', lib="/home/jcayford/r_libs")


# Functions
        # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}
        # Rounding data
                round_to <- function(x, to = window_length) round(x/to)*to

        # Function to generate the necessary interactions and fragments files for input into FitHiC
                # Chr_num is coded to be a chr number in a for loop but it can be done separately
                # Status is either "healthy" or "disease". No other inputs will work
                        frags_inters_func <- function(chr_num, status){
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
                                interactions<-data.frame(chromosome, mids, chromosome, (df[,2]+2500), round_to(df[,3], to=1))
                                colnames(interactions)<-c("chr1", "fragmentMid1", "chr2", "fragmentMid2", "contactCount")
                                interactions[is.na(interactions)]=0 
                                write.table(interactions, paste(status, "_interactions_", chromosome, "_", normalization_caps, ".txt", sep=""), col.names=F, row.names=F, quote=F)         
                        # Generation of the fragments data.frame
                                fragments_pre <- data.frame(df[,1]+2500, df[,3])
                                colnames(fragments_pre) <- c("fragmentMid", "ContactCount")
                                fragments_pre[is.na(fragments_pre)]=0
                                fragments_mid <- fragments_pre %>% group_by(fragmentMid) %>% summarise(ContactCount=sum(ContactCount))
                                fragments <- data.frame("chr1"=chromosome, "extra"=1, round_to(fragments_mid, to=1), "extra.1"=1)
                                fragments[is.na(fragments)]=0          
                                write.table(fragments, paste0(status, "_fragments_", chromosome, "_", normalization_caps, ".txt"), col.names=F, row.names=F, quote=F)                   
                                
                        }

#
          
# Normalization. Generally will be either kr or none.
normalization <- "kr"

# Setting the import and export directories

import_dir <- paste0("/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/juicer_files")
export_dir <- paste0("/home/jcayford/HCM/Analysis/Fithic_Raw/kr_bias_10may2020/")

for(j in 1:19){
    print(paste(Sys.time(), "............ working on healthy chromosome: ", j, sep=""))
        frags_inters_func(j, "healthy")
}



for(j in 1:19){
    print(paste(Sys.time(), "............ working on disease chromosome: ", j, sep=""))
        frags_inters_func(j, "disease")
}


for(j in 1:19){
    print(paste(Sys.time(), "............ working on ctcf_ko chromosome: ", j, sep=""))
        frags_inters_func(j, "ctcf_ko")

}
