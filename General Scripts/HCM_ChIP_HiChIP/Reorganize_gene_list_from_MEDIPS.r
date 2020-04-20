
### Reorganizing the gene list do the +/- strands
            re_org <- function(row, df, distance){
                if(df[i,2]<df[i,3]){
                    
                    a <- data.frame("chr"=df[i, 1], "start"=(df[i, 3]-distance), "end"=(df[i, 2]+distance), "gene"=df[i,4], "strand"="-")
                } else{
                    a <- data.frame("chr"=df[i, 1], "start"=(df[i, 2]-distance), "end"=(df[i, 3]+distance), "gene"=df[i,4], "strand"="+")
                }
                a
            }

            rat_genes <- read.table("refGene.txt")
            gene_list <- 

            distance <- 0
            ordered_genes <- NULL
            for(i in 1:nrow(gene_list)){
                ordered_genes <- rbind(ordered_genes, re_org(i, gene_list, distance))
            }
            genes <- data.frame("chr"=ordered_genes[,1], "start"=ordered_genes[,2]+10000, "end"=ordered_genes[,3]-10000, "gene"=ordered_genes[,4])
            genes2 <- data.frame(genes[,1], genes[,3], genes[,2], genes[,4])
            colnames(genes2) <- c("chr", "start", "end", "gene")
            write.table(genes2, paste("refGene_ordered.txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)
###

