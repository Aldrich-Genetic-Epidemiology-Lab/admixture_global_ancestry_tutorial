#The tutorial found on https://www.gnxp.com/WordPress/2018/07/13/tutorial-to-run-supervised-admixture-analyses/?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+RazibKhansTotalFeed+%28Razib+Khan%27s+total+feed%29 was helpful in generating this workflow.

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
SOURCE_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids
TARGET_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids/1kg_filt_ceu_yri
FAM_NAME=all_phase3

#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from either the ASW or CEU 
grep "CEU\|YRI" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/keep.keep

#Use plink to generate an entire new set of binary files containing just the CEU and YRI 1000 Genomes individuals
plink --bfile $SOURCE_DIR\/$FAM_NAME --keep $TARGET_DIR\/keep.keep --make-bed --allow-extra-chr --out $TARGET_DIR\/ceu_yri_phase3
