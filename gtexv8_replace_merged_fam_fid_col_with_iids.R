#Declare the path of the fam file (for plink set of SCCS cohort merged with CEU and YRI 1000 Genomes individuals)
fam_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.fam"

#Open the fam file as a data frame
fam_file <- read.table(fam_path, header = FALSE)
fam_df <- as.data.frame(fam_file)

#Make a new data frame in which the FID column of all 0s is replaced with the IID column (column 2)
new_df <- data.frame(fam_df[,2], fam_df[,2], fam_df[,3], fam_df[,4], fam_df[,5], fam_df[,6])

pops <- new_df[,1]
iids <- new_df[,2]

#Replace the family IDs for founder populations with their 3-letter code
pops[startsWith(pops, "CEU")] <- "CEU"
pops[startsWith(pops, "YRI")] <- "YRI"

#Relace the IIDs for founder populations with the IID code only (removing population prefix)
iids[startsWith(pops, "CEU") | startsWith(pops, "YRI")] <- substring(iids[startsWith(pops, "CEU") | startsWith(pops, "YRI")], 5)

final_df <- data.frame(pops, iids, new_df[,3], new_df[,4], new_df[,5], new_df[,6])

#Rewrite the existing fam file with this new data frame
write.table(final_df, file = fam_path, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE) 