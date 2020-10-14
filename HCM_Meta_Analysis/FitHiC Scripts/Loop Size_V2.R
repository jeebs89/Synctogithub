library('tidyverse', lib="/home/jcayford/r_libs")
library(reshape2)
    # Functions
            # Quick head() with the amount of rows
                w <- function(x){print(head(x));nrow(x)}                


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
    #
    #

    hic_data_directory <- "/home/jcayford/HCM/Analysis/HiC_R/chromosomes_all_kr"
    deseq_const_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq/Loop_Constraints"
 
# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")

loop_constraint <- 0.5

# Importing the looping data for MEF2a
    setwd(hic_data_directory)
    all_loops_a <- read.table("hic_contacts_mef2a_all_merged_q0.005_28sep2020.txt", header=T)
    all_loops <- all_loops_a %>% filter(cond=="disease" | cond=="healthy" | cond=="h_d") %>% select(c(1:8, 10:12,16))
    loop_size <- data.frame(all_loops, "loop_size"=(all_loops[,4]-all_loops[,2]))
    all_size <- loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(loop_size)))

# Importing the filtered Loops based on RNA-Seq
    setwd(deseq_wd)
    filtered_loop_size_a <- read.table("All_loop_rna_overlaps_master_0mb_limit_14_Oct_2020.txt", header=TRUE)
    filtered_loop_size <- data.frame(filtered_loop_size_a, "loop_size"=(filtered_loop_size_a[,4]-filtered_loop_size_a[,2]))
    all_t_size <- filtered_loop_size %>% arrange(loop_size) %>% select(c(2, 4, ncol(filtered_loop_size)))



loop_constraint <- c(0, 1.5, 1, 0.5)

looping_df <- NULL
for(i in 1:length(loop_constraint)){
    loops <- NULL;loops_a <- NULL
        if(loop_constraint[i] == 0){loops <- data.frame(all_t_size, "State"="Post_Filter_0")}
        if(loop_constraint[i] > 0){
                    loops_a <- all_t_size %>% filter(loop_size < (loop_constraint[i]*1000000))
                    loops <- data.frame(loops_a, "State"=paste0("Post_Filter_", loop_constraint[i]))
        }

        looping_df <- rbind(looping_df, loops)
}




# Graphing the density of the loops
a <- data.frame("State"="Pre_Filter", "loop_size"=all_size[,3])
b <- data.frame("State"=looping_df[,4], "loop_size"=looping_df[,3])
c <- rbind(a, b)


both <- melt(c, id=c("loop_size"))

ggplot(both, aes(x=(loop_size/1000), fill=value, col=value)) +
    geom_density(alpha=0.1) + 
    theme +
    #scale_color_manual(values=c("#808080", "#FF8F1E")) +
    labs(x="Loop Size (kb)", y="Loop Size Frequency", main="Density of Pre/Post Filtering of Loop Size")





ggsave(paste0("pre_post_loop_frequency_plot_all_mb_limit_", today_f, ".png"))





