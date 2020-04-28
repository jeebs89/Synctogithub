# Script to run the MEDIPS for the JC_ChIP. This is not a clean copy/paste script.


library(MEDIPS, , lib="/home/jcayford/r_libs")
library(BSgenome.Hsapiens.UCSC.hg38, lib="/home/jcayford/r_libs")
library('tidyverse', lib="/home/jcayford/r_libs")
library(reshape2)


# Functions
        # Rounding data
            round_to <- function(x, to = window_length) round(x/to)*to
        # Quick head() with the amount of rows
            w <- function(x){print(head(x));nrow(x)}


        # Function to merge 5 samples together and export the file
            merge_5_samples <- function(input_df, out_name){
                    MSets_group <- input_df

                    merged_set=MEDIPS.mergeSets(MSet1=MSets_group[[1]], MSet2=MSets_group[[2]], name="Merged Set")
                    merged_set1=MEDIPS.mergeSets(MSet1=MSets_group[[3]], MSet2=MSets_group[[4]], name="Merged Set_1")
                    merged_set2=MEDIPS.mergeSets(MSet1=merged_set, MSet2=MSets_group[[5]], name="Merged Set_2")
                    merged_set3=MEDIPS.mergeSets(MSet1=merged_set1, MSet2=merged_set2, name="Merged Set_3")                     
                    MEDIPS.exportWIG(Set=merged_set3, file=paste0(out_name, ".wig"), format="rpkm", descr=out_name)

            }

        # Function to generate the correlations. This generates the plots if needed and exports the values
            correlation_func <- function(input_df, output, plots, methods){
                corr <- MEDIPS.correlation(MSets=input_df, plot = plots, method=methods)
                write.table(corr, paste0(output, "_correlation_", methods, ".txt"), row.names=T, col.names=T, quote=FALSE)
                corr
            }

        # Generation of heatmaps 
            heatmap_func <- function(input_df, out_name, method){
                melted_df <- melt(input_df, na.rm=TRUE)
                ggplot(melted_df, aes(x=Var2, y=Var1, fill=value)) +
                    geom_tile(color="white") +
                    scale_color_brewer(type = 'qual', palette = 'Dark2') +
                    theme +
                    geom_text(aes(Var2, Var1, label = round_to(value, 0.01)), color = "black", size = 8)
                    #ggsave(paste0(out_name, "_", method, "_heatmap.png"))
            }

        # ggplot theme
            theme<-theme_classic() +
            theme(
                text = element_text(color = "grey20"),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(size=rel(1.5)), 
                axis.text.y = element_text(size=rel(1.5), angle=90, hjust=0.5)
            )


#

# Input conditions
    extend = 300
    shift = 0
    window_size = 500
    BSgenome = "BSgenome.Hsapiens.UCSC.hg38"
    uniq = 1
    chr.select = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22")
    paired = F
    samples <- read.csv("medips_input_27apr2020.csv", header=TRUE, stringsAsFactors=F); colnames(samples) <- c("Sample", "File", "C1")


# Selection of samples from the sample sheet
    cd4 <- samples %>% filter(C1=="cd4") %>% select(Sample, File)
    mono <- samples %>% filter(C1=="mono") %>% select(Sample, File)
    nk <- samples %>% filter(C1=="nk") %>% select(Sample, File)


# Creation of the sets for MEDIPS
    MSets_groupA = NULL
    for(i in 1:nrow(cd4)){
    MSets_groupA = c(MSets_groupA, MEDIPS.createSet(file=cd4[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=cd4[i,1]))
    gc()
    }

    MSets_groupB = NULL  
    for(i in 1:nrow(mono)){
    MSets_groupB = c(MSets_groupB, MEDIPS.createSet(file=mono[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=mono[i,1]))
    gc()
    }

    MSets_groupC = NULL  
    for(i in 1:nrow(nk)){
    MSets_groupC = c(MSets_groupC, MEDIPS.createSet(file=nk[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=nk[i,1]))
    gc()
    }


# Merging all of the samples together
    merge_5_samples(MSets_groupA, "cd4_merge")
    merge_5_samples(MSets_groupB, "monocytes_merge")
    merge_5_samples(MSets_groupC, "nk_merge")

# Pearson Correlations
    setA_correlation <- MEDIPS.correlation(MSets=MSets_groupA, plot = FALSE, method="pearson")
    setB_correlation <- MEDIPS.correlation(MSets=MSets_groupB, plot = FALSE, method="pearson")
    setC_correlation <- MEDIPS.correlation(MSets=MSets_groupC, plot = FALSE, method="pearson")

    setA_correlation_s <- MEDIPS.correlation(MSets=MSets_groupA, plot = FALSE, method="spearman")
    setB_correlation_s <- MEDIPS.correlation(MSets=MSets_groupB, plot = FALSE, method="spearman")
    setC_correlation_s <- MEDIPS.correlation(MSets=MSets_groupC, plot = FALSE, method="spearman")






 
                melted_df <- melt(setC_correlation, na.rm=TRUE)
                
                ggplot(melted_df, aes(x=Var2, y=Var1, fill=value)) +
                    geom_tile(color="white") +
                    scale_fill_gradient(
                                        name = "Cor", # changes legend title
                                        low = "grey",
                                        high = "gold3",
                                        limit = c(0.9, 1),
                                        space = "Lab",
                                        guide = "colourbar"
                                    ) +
                    theme +
                    geom_text(aes(Var2, Var1, label = round_to(value, 0.01)), color = "black", size = 8)
                    












            heatmap_func <- function(input_df, out_name, method, color){
                melted_df <- melt(setA_correlation, na.rm=TRUE)
                ggplot(melted_df, aes(x=Var2, y=Var1, fill=value)) +
                    geom_tile(color="white") +
                    scale_fill_gradient(
                                        name = "Cor", # changes legend title
                                        low = "grey",
                                        high = color,
                                        limit = c(0.9, 1),
                                        space = "Lab",
                                        guide = "colourbar"
                                    ) +
                    theme +
                    geom_text(aes(Var2, Var1, label = round_to(value, 0.01)), color = "black", size = 8)
                    ggsave(paste0(out_name, "_", method, "_heatmap_new.png"))
            }




    heatmap_func(setA_correlation, "cd4", "pearson", "dodgerblue3")
    heatmap_func(setB_correlation, "monocytes", "pearson", "firebrick3")
    heatmap_func(setC_correlation, "nk", "pearson", "gold3")






# Spearman Correlations
    correlation_func(MSets_groupA, "cd4", FALSE, "spearman")
    correlation_func(MSets_groupB, "monocytes", FALSE, "spearman")
    correlation_func(MSets_groupC, "nk", FALSE, "spearman")
