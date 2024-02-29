#Declare the path of the input bed. Also within this directory should be the accompanying bim, fam, and pop files
in_bed=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.bed

#Run ADMIXTURE using the --supervised flag and a K of 2, as we are working with two known reference populations from 1000 Genomes in this analysis
admixture --supervised $in_bed 2