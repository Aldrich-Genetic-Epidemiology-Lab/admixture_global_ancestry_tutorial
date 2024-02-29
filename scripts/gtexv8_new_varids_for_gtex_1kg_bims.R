library("data.table")

#Declare the paths of the GTEx v8 lifted bim file, 1000 Genomes bim file, and the output directory
bim_path_gtex <- "/home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.bim"
bim_path_1kg <- "/home/bettimj/gamazon_rotation/1000_genomes/b37/all_phase3.bim"
out_dir <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids"

out_name_gtex <- "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.bim"
out_name_1kg <- "all_phase3.bim"

#Open the bim files as data frames
bim_file_gtex <- fread(bim_path_gtex, header = FALSE, quote = "", sep = "\t")
bim_df_gtex <- as.data.frame(bim_file_gtex)

bim_file_1kg <- fread(bim_path_1kg, header = FALSE, quote = "", sep = "\t")
bim_df_1kg <- as.data.frame(bim_file_1kg)

#Make new varIDs for each of the bims
bim_df_gtex[,2] <- paste(paste0("chr", bim_df_gtex[,1]), bim_df_gtex[,4], bim_df_gtex[,6], bim_df_gtex[,5], "b37", sep = "_")

bim_df_1kg[,2] <- paste(paste0("chr", bim_df_1kg[,1]), bim_df_1kg[,4], bim_df_1kg[,6], bim_df_1kg[,5], "b37", sep = "_")

#Write out each of the new bims to the specified output directory
write.table(bim_df_gtex, file = paste(out_dir, out_name_gtex, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(bim_df_1kg, file = paste(out_dir, out_name_1kg, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)