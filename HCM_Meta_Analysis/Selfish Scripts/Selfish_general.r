# Selfish


./selfish/selfish/selfish.py -f1 /path/to/contact/map1.txt \
                             -f2 /path/to/contact/map2.txt \
	                     -ch 2 \
                             -r 100kb -o ./output.npy






# Differences between chr1 in healthy and disease

./selfish/selfish/selfish.py -f1 /home/jcayford/HCM/Analysis/Fithic_Raw/disease/disease_interactions_chr1_NONE.txt -f2 /home/jcayford/HCM/Analysis/Fithic_Raw/healthy/healthy_interactions_chr1_NONE.txt -ch 1 -r 20kb -o test.npy

./selfish/selfish/selfish.py -f1 /home/jcayford/HCM/Analysis/Fithic_Raw/disease/disease_interactions_chr1_NONE.txt -f2 /home/jcayford/HCM/Analysis/Fithic_Raw/healthy/healthy_interactions_chr1_NONE.txt -ch 1 -r 20kb -o test.npy




/home/jcayford/selfish/selfish/selfish.py -f1 /home/jcayford/HCM/Analysis/Fithic_Raw/disease/disease_interactions_chr1_NONE.txt -f2 /home/jcayford/HCM/Analysis/Fithic_Raw/healthy/healthy_interactions_chr1_NONE.txt -ch 1 -r 5kb -t 0.001 -plot=True -d 5000000 -o chr1_out_5kb_0.001.tsv





library('tidyverse', lib="/home/jcayford/r_libs")


# Functions
        # Rounding data
            round_to <- function(x, to = window_length) round(x/to)*to
        # Quick head() with the amount of rows
            w <- function(x){print(head(x));nrow(x)}
  


selfish_data <- read.table("selfish_all_chr_085_0.01.tsv", header=T)

selfish_data_chr1 <- read.table("chr1_out_0.1.tsv", header=T, stringsAsFactors=F)

test <- selfish_data_chr1 %>% arrange(P_VAL)

test2 <- data.frame(test, "distance"=(test$LOC2 - test$LOC1))

test3 <- test2 %>% filter(distance < 15000000)
test3 <- test3 %>% arrange(distance)
