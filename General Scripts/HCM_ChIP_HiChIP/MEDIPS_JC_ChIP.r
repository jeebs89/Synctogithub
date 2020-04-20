# Script to run the MEDIPS for the JC_ChIP. This is not a clean copy/paste script.

library(MEDIPS, , lib="/home/jcayford/r_libs")
library(BSgenome.Rnorvegicus.UCSC.rn6, lib="/home/jcayford/r_libs")
library('tidyverse', lib="/home/jcayford/r_libs")


extend = 0
shift = 0
window_size = 500
BSgenome = "BSgenome.Rnorvegicus.UCSC.rn6"
uniq = 1
#chr.select = c("chr11")
chr.select = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20")
paired = T




samples<-read.csv("medips_hcm_30mar20.csv", header=TRUE, stringsAsFactors=F); colnames(samples) <- c("Sample", "File", "C1", "C2")


# Rounding data
            round_to <- function(x, to = window_length) round(x/to)*to
        # Quick head() with the amount of rows
            w <- function(x){print(head(x));nrow(x)}
  



df_subset <- function(condition_a, condition_b){
    samples %>% filter(C1==condition_a, C2==condition_b)
}


serum_ac <- samples %>% filter(C1=="serum" | C1=="serum_dmso", C2=="ac") %>% select(Sample, File)
no_serum_ac <- samples %>% filter(C1=="no_serum" | C1=="no_serum_dmso", C2=="ac") %>% select(Sample, File)

serum_me <- samples %>% filter(C1=="serum" | C1=="serum_dmso", C2=="me3") %>% select(Sample, File)
no_serum_me <- samples %>% filter(C1=="no_serum" | C1=="no_serum_dmso", C2=="me3") %>% select(Sample, File)


MSets_groupA = NULL

for(i in 1:nrow(no_serum_ac)){
  MSets_groupA = c(MSets_groupA, MEDIPS.createSet(file=no_serum_ac[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=no_serum_ac[i,1]))
  gc()
  }

      #merged_setA=MEDIPS.mergeSets(MSet1=MSets_groupA[[1]], MSet2=MSets_groupA[[2]], name="Merged Set_A")
      #merged_setA1=MEDIPS.mergeSets(MSet1=MSets_groupA[[3]], MSet2=MSets_groupA[[4]], name="Merged Set_A1")
      #merged_setA3=MEDIPS.mergeSets(MSet1=merged_setA, MSet2=merged_setA1, name="Merged Set_A2")
      
      #MEDIPS.exportWIG(Set=merged_setA3, file="all_no_serum_merge_me3.wig", format="rpkm", descr="all_no_serum_merge_me3")


MSets_groupB = NULL  
for(i in 1:nrow(serum_ac)){
  MSets_groupB = c(MSets_groupB, MEDIPS.createSet(file=serum_ac[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=serum_ac[i,1]))
  gc()
}

      #merged_setB=MEDIPS.mergeSets(MSet1=MSets_groupB[[1]], MSet2=MSets_groupB[[2]], name="Merged Set_B")
      #merged_setB1=MEDIPS.mergeSets(MSet1=MSets_groupB[[3]], MSet2=MSets_groupB[[4]], name="Merged Set_B1")
      #merged_setB3=MEDIPS.mergeSets(MSet1=merged_setB, MSet2=merged_setB1, name="Merged Set_B2")
      
      #MEDIPS.exportWIG(Set=merged_setB3, file="all_serum_merge_me3.wig", format="rpkm", descr="all_serum_merge_me3")


MSets_groupC = NULL  
for(i in 1:nrow(no_serum_me)){
  MSets_groupC = c(MSets_groupC, MEDIPS.createSet(file=no_serum_me[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=no_serum_me[i,1]))
  gc()
}

MSets_groupD = NULL  
for(i in 1:nrow(serum_me)){
  MSets_groupD = c(MSets_groupD, MEDIPS.createSet(file=serum_me[i,2], extend=extend, window_size=window_size, BSgenome=BSgenome, uniq=uniq, chr.select=chr.select, paired=paired, sample_name=serum_me[i,1]))
  gc()
}






MEDIPS.correlation(MSets=MSets_groupA, plot = T, method="spearman")









#MinRowSum_Cond1_Cond2 = (nrow(bam.files.groupA) + nrow(bam.files.groupA.input) + nrow(bam.files.groupB) + nrow(bam.files.groupB.input)) * 8
#MinRowSum_Cond1_Cond2

res_1_T = MEDIPS.meth(MSet1 = MSets_groupA, MSet2 = MSets_groupB, chr = chr.select, p.adj = "BH", diff.method = "edgeR", CNV = FALSE, MeDIP = F)
res.s_1 = MEDIPS.selectSig(res_1_T, p.value=0.01, adj=T)
res.s_filtered<-MEDIPS.mergeFrames(frames=res.s_1, distance=10000)


setwd("/home/jcayford/HCM/jc_chip/MEDIPS")
write.table(res.s_1, "merged_serum_noserum_ac_all_500bp_bin_10kb_dist_pvalue_001.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=',')
#write.table(res_1_T, "merged_serum_noserum_ac_all_500bp_bin_10kb_dist_Counts_Full_genome.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=',')


write.table(res.s_filtered, "merged_serum_noserum_ac_all_500bp_bin_10kb_dist_Counts_Full_genome_pval001.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')






res_2_T = MEDIPS.meth(MSet1 = MSets_groupC, MSet2 = MSets_groupD, chr = chr.select, p.adj = "BH", diff.method = "edgeR", CNV = FALSE, MeDIP = F)
res.s_2 = MEDIPS.selectSig(res_2_T, p.value=0.1, adj=T)
res.s_filtered_2<-MEDIPS.mergeFrames(frames=res.s_2, distance=10000)


setwd("/home/jcayford/HCM/jc_chip/MEDIPS")
write.table(res.s_2, "merged_serum_noserum_me3_all_500bp_bin_10kb_dist_pvalue_001.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=',')
#write.table(res_2_T, "merged_serum_noserum_me3_all_500bp_bin_10kb_dist_Counts_Full_genome.csv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep=',')


write.table(res.s_filtered_2, "merged_serum_noserum_me3_all_500bp_bin_10kb_dist_Counts_Full_genome_pval001.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')









no_serum_peaks <- res.s_1 %>% filter(MSets1.rpkm.mean>MSets2.rpkm.mean)
res.s_filtered<-MEDIPS.mergeFrames(frames=no_serum_peaks, distance=1000)
write.table(res.s_filtered, "merged_noserum_ac_peaks_all_500bp_bin_1kb_dist_Counts_Full_genome_pval0.01.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')


serum_peaks <- res.s_1 %>% filter(MSets1.rpkm.mean<MSets2.rpkm.mean)
res.s_filtered<-MEDIPS.mergeFrames(frames=serum_peaks, distance=1000)
write.table(res.s_filtered, "merged_serum_ac_peaks_all_500bp_bin_1kb_dist_Counts_Full_genome_pval0.01.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')






res.s_2 = MEDIPS.selectSig(res_2_T, p.value=0.2, adj=T)

no_serum_peaks <- res.s_2 %>% filter(MSets1.rpkm.mean>MSets2.rpkm.mean)
res.s_filtered<-MEDIPS.mergeFrames(frames=no_serum_peaks, distance=1000)
write.table(res.s_filtered, "merged_noserum_me_peaks_all_500bp_bin_1kb_dist_Counts_Full_genome_pval0.01.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')


serum_peaks <- res.s_2 %>% filter(MSets1.rpkm.mean<MSets2.rpkm.mean)
res.s_filtered<-MEDIPS.mergeFrames(frames=serum_peaks, distance=1000)
write.table(res.s_filtered, "merged_serum_me_peaks_all_500bp_bin_1kb_dist_Counts_Full_genome_pval0.01.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')


