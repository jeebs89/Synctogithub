## Script to find genes which are nearby to the regions picked out from MEDIPs analysis


library('tidyverse', lib="/home/jcayford/r_libs")



setwd("/home/jcayford/HCM/jc_chip/MEDIPS/outputs")
  df <- read.table("merged_serum_noserum_all_1kb_bin_10kb_dist_Counts_Full_genome.csv", header=T, sep=",")
  chip_peaks <- read.table("merged_serum_noserum_all_1kb_bin_10kb_dist_pvalue_001.csv", header=T, sep=",")

setwd("/home/jcayford/HCM/jc_chip/rat_genes")


unordered_genes <- read.table("refGene.txt", header=FALSE)
rat_genes <- read.table("refGene_ordered.txt", header=TRUE)

rat_genes <- data.frame(rat_genes, "strand"=unordered_genes[,4])


distance <- 20000


genes_all <- data.frame("chr" = rat_genes[,1], "start" = rat_genes[,2]-distance, "end" = rat_genes[,3]+distance, "gene" = rat_genes[,4], "strand"=rat_genes[,5])
chip_mid <- data.frame(chip_peaks, round_to(rowMeans(chip_peaks[,2:3]), to=500));colnames(chip_mid)[28]<-"mid"


full_genes <- NULL
for(j in 1:20){
    chromosome <- paste("chr", j,sep="")
    genes_chr <- genes_all %>% filter(chr == chromosome) %>% arrange(start)
    chip_chr <- chip_mid %>% filter(chr == chromosome) 

    print(paste(Sys.time(), "..........starting chromosome ", j, sep=""))
        genes <- NULL
        for(i in 1:(nrow(chip_chr))){
            low_range <- NULL
            
            if(i==1){
                low_range <- genes_chr %>% filter(start < chip_chr[i, 28])
                genes_r <- low_range %>% filter(end > chip_chr[i, 28] & start < chip_chr[i, 28])     
                    if(nrow(genes_r)==0){
                         genes <- NULL
                    } else{genes <- data.frame(genes_r, "mid"=chip_chr[i, 28])}
                gene_filter <- genes_chr %>% filter(start > chip_chr[i, 28])
            }

            if(i > 1){
                low_range <- gene_filter %>% filter(start < chip_chr[i, 28])
                gene_r2 <- low_range %>% filter(end > chip_chr[i, 28] & start < chip_chr[i, 28])

                    if(nrow(gene_r2)==0){
                        gene_a <- 0
                    } else{gene_a <- data.frame(gene_r2, "mid"=chip_chr[i, 28])}

                gene_filter <- genes_chr %>% filter(start > chip_chr[i, 28])
            }

            if(length(gene_a)==1){
                genes <- genes
            } else{genes <- rbind(genes, gene_a)}
        }

    full_genes <- rbind(full_genes, genes)
}


test <- full_genes %>% distinct(gene, .keep_all=TRUE)
test2 <- inner_join(test, chip_mid, by=c("chr", "mid"))
test3 <- test2 %>% arrange(edgeR.adj.p.value) 


pos_strand <- test3 %>% filter(strand=="+")

minus <- data.frame(pos_strand, "sub"=(pos_strand$start.x - pos_strand$mid))
minus_ordered <- minus %>% arrange(desc(sub))



write.table(test3, "DER_genes_10kb_window_pval_0.01.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)




