#Declare the path of the Q file, which was generated from the initial ADMIXTURE analysis
q_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.Q"

#Read in the Q file as a data frame
q_file <- read.table(q_path, header = FALSE)
q_df <- as.data.frame(q_file)

#Append the individual IDs from the fam file to the Q file
fam_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.fam"
fam_file <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
fam_df <- as.data.frame(fam_file)

q_df <- data.frame(fam_df[,1], q_df[,1:2])

#Find the number of 1000 Genomes reference individuals in the .pop file used in the ADMIXTURE analysis
pop_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.pop"
pop_file <- read.table(pop_path, header = FALSE)
pop_df <- as.data.frame(pop_file)
pop_df_nonref_count <- length(pop_df[(pop_df[,1] == "-"),])

#This number turns out to be 601-807. So we will exclude these individuals, leaving only the unlabeled SCCS cohort individuals
cohort_only_q_df <- q_df[(1:pop_df_nonref_count),]

#Write the new data frame out to a new .Q file
out_name <- "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.anno.pruned.Q"
write.table(cohort_only_q_df, file = paste("/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri", out_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)