



##### TEST2 is from the DESeq_Looping_analysis_V2 and the subs is from the FitHiC_overlaps_3cond.

post_a <- data.frame(test2[,c(1, 9:10, 17)], "pre_post"="post")
post <- data.frame(table(post_a$peak), "pre_post"="post")
post_multi_contacts <- post %>% filter(Freq>1)

pre_a <-  data.frame(subs[,c(1:2, 4, 10)], "pre_post"="pre")
pre <- data.frame(table(pre_a$peak), "pre_post"="pre")
pre_multi_contacts <- pre %>% filter(Freq>1)






##### NEED TO GET THE MEF2:MEF2 INTERACTIONS BY LOOKING AT THE CHIP-SEQ PEAK TO FIND IT IN EITHER CONTACT 1 OR 2. 
##### THIS FEATURE IS REMOVED FROM THE ANALYSIS IN FitHiC_overlaps_3cond LINE 177. CAN LIKELY REMOVE THE DISTINCT
##### FEATURE IN ORDER TO FIND OUT WHERE THERE ARE MULTIPLE CONTACTS ONLY AND THEN USE THAT TO FIND THE ENRICHMENT.
##### THERE SHOULD BE A REDUCTION OF MEF2:MEF2 CONTACTS IN HCM.




both <- rbind(post, pre)



ggplot(both, aes(x=Var1, col=pre_post)) +
    geom_density() + 
    theme +
    scale_color_manual(values=c("#808080", "#FF8F1E")) +
    labs(x="Mustache Q-value", y="Q Value Frequency")



ggsave(paste0("mustache_qvalue_density_plot_", today_f, ".png"))






