

        library('tidyverse', lib="/home/jcayford/r_libs")
        library('reshape')
        library('viridisLite')


# Today's date 
today <- Sys.Date()
today_f <- format(today, format="%d_%b_%Y")

  
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





    deseq_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq"
    export_wd <- "/mnt/BioAdHoc/Groups/vd-vijay/justin/HCM_Meta/rna-seq/DESeq/Gene_Count_Plots"



    setwd(deseq_wd)

    genes <- read.table("Top_loop_rna_overlaps_TPM_COUNTS06_Oct_2020.txt", header=TRUE, sep=" ")
    
    dots_a <- genes[c(1, 14:19)]
    dots_melt <- melt(dots_a, id.vars="mgi_symbol")
      
    dots_melt_healthy <- dots_melt %>% filter(variable=="tpm_h1" | variable=="tpm_h2" | variable=="tpm_h3")
    dots_h <- data.frame("mgi_symbol"=dots_melt_healthy[,1], "State"="Healthy", "TPM"=log2(dots_melt_healthy[,3]+1))
 
    dots_melt_disease <- dots_melt %>% filter(variable=="tpm_d1" | variable=="tpm_d2" | variable=="tpm_d3")
    dots_d <- data.frame("mgi_symbol"=dots_melt_disease[,1], "State"="HCM", "TPM"=log2(dots_melt_disease[,3]+1))

    dots_melted <- rbind(dots_h, dots_d)


   

    setwd(export_wd)
    for(i in 1:nrow(genes)){
        
     gene_test <- dots_melted %>% filter(mgi_symbol==genes[i,1])

    plot <- ggplot(gene_test, aes(x=State, y=TPM, fill=State)) + 
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.1), cex=2) + 
        theme +
        scale_fill_manual(values=c("dodgerblue", "firebrick3")) +
        labs(y="Log2(TPM+1)", title=paste0(genes[i,1], " TPM Counts Healthy vs Disease")) +
        expand_limits(y=0)
    ggsave(plot, filename=paste0(genes[i,1], "_count_plot_tpm_", today_f, ".png"))

    }

    
