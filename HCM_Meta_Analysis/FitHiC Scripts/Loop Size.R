library('tidyverse', lib="/home/jcayford/r_libs")
library(reshape2)
    # Functions
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}                

    #

    hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
 
# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")


# Importing the looping data for MEF2a
    setwd(hic_data_directory)
    all_loops_a <- read.table("hic_contacts_mef2a_all_merged_q0.005_28sep2020.txt", header=T)
    all_loops <- all_loops_a %>% filter(cond=="disease" | cond=="healthy" | cond=="h_d") %>% select(c(1:8, 10:12,16))
    loop_size <- data.frame(all_loops, "loop_size"=(all_loops[,4]-all_loops[,2]))
    all_size <- loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(loop_size)))

# Importing the filtered Loops based on RNA-Seq
    setwd(deseq_wd)
    filtered_loop_size_a <- read.table("All_loop_rna_overlaps_master_7_OCT_20.txt", header=TRUE, rown.names=TRUE)
    filtered_loop_size <- data.frame(filtered_loop_size_a, "loop_size"=(filtered_loop_size_a[,4]-filtered_loop_size_a[,2]))
    all_t_size <- filtered_loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(filtered_loop_size)))

# Graphing the density of the loops
a <- data.frame("State"="Pre_Filter", "loop_size"=all_size[,3])
b <- data.frame("State"="Post_Filter", "loop_size"=all_t_size[,3])
c <- rbind(a, b)
colnames(c)[3] <- "HiC Loops"

both <- melt(c, id=c("loop_size"))

ggplot(both, aes(x=(loop_size/1000), col=value)) +
    geom_density() + 
    theme +
    scale_color_manual(values=c("#808080", "#FF8F1E")) +
    labs(x="Loop Size (kb)", y="Loop Size Frequency", main="Density of Pre/Post Filtering of Loop Size")
ggsave(paste0("pre_post_loop_frequency_plot_", today_f, ".png"))
