#Declare the path of the pruned Q ADMIXTURE output file containing ancestry estimates for the SCCS cohort compared with 1000 Genomes CEU and YRI individuals
pruned_q_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.anno.pruned.Q"

#Open the Q file and read it in as a data frame
pruned_q_file <- read.table(pruned_q_path, header = FALSE)
pruned_q_df <- as.data.frame(pruned_q_file)

#Create a new data frame ordering the individuals from highest to lowest proportion of African ancestry
ordered_pruned_q_df <- pruned_q_df[order(-pruned_q_df[,3]),]

#Write the ordered data frame out to a new file
out_name <- "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.anno.pruned.ordered.Q"
write.table(ordered_pruned_q_df, file = paste("/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri", out_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)