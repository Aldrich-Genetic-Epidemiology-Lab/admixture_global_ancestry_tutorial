#Declare the file path for the merged fam file with select founder populations and cohort individuals labeled
fam_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.fam"

#Open the fam file as a data frame and slice out the first column containing founder populations
fam_file <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
fam_df <- as.data.frame(fam_file, stringsAsFactors = FALSE)

pops <- fam_df[,1]

#For any individual that does not have a population listed as ASW or CEU, replace the value in pops with "-"
is_not_founder <- !(pops == "CEU" | pops == "YRI")

pops[is_not_founder] <- "-"

#Write out the pops data frame to a .pop file
out_name <- "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.pop"
write.table(pops, file = paste("/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri", out_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)