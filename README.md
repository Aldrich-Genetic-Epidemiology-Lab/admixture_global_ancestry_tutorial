# Estimating global genetic ancestry using ADMIXTURE
## Installing dependencies
The following dependencies are required to run this entire workflow:
*  plink
*  plink2
*  Picard
*  R
*  the data.table R package
*  ADMIXTURE

All of these dependencies can easily be installed via miniconda using the included ```admixture_env.yml``` file. If you do not already have miniconda3 installed on your ACCRE account, first download the install script to your home directory on ACCRE (more information at https://docs.conda.io/en/latest/miniconda.html):
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Next modify the script's permissions, giving it all read/write access:
```
chmod 777 Miniconda3-latest-Linux-x86_64.sh
```

Run the install script and follow the terminal prompts:
```
./Miniconda3-latest-Linux-x86_64.sh
```

There should be an indicator next to your command line prompt saying ```(base)```, which indicates that your installation of miniconda is active. It should look something like the following:
```
(base) Michaels-Mini:admixture_global_ancestry michaelbetti$
```

If you do not see this indicator, then run the following command to refresh your bashrc profile, which should activate miniconda:
```
source ~/.bashrc
```

Once miniconda is actively running, you can run the following command to create a miniconda environment for running this ADMIXTURE workflow:
```
conda env create -f admixture_env.yml
```

Once the environment has finished compiling, you can activate it using the following command:
```
conda activate admixture_env
```

You should now see that the ```(base)``` indicator on your command line has changed to ```(admixture_env)```, looking something like the following:
```
(admixture_env) Michaels-Mini:admixture_global_ancestry michaelbetti$ 
```

If you would like to deactivate this environment, run the following command:
```
conda deactivate
```

## Required datasets
In order to initiate a supervised run of ADMIXTURE, the following two core datasets will be required:
*  A set of variants (in plink format) from reference individuals on which the global ancestry predictions will be based (such as 1000 Genomes populations)
*  The variants (in plink format) from the target individual, for whom you aim to calculate global ancestry

## Compiling reference dataset
The full download instructions for compiling a set of 1000 Genomes variant data in plink format can be found at https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf. However, these commands should work:
```
pgen=https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1

wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

wget $sample
mv 'phase3_corrected.psam?dl=1' all_phase3.psam

plink2 --pfile all_phase3 --max-alleles 2 --make-bed --out all_phase3
```

The resulting plink binary files (bed/bim/fam) can be saved to a directory of your choosing.

*Note that this 1000 Genomes dataset is based in human genome build GRCh37/hg19. The lifted-over hg38 dataset has been pulled from the 1000 Genomes FTP for quality issues, and while possible, it would be a more laborious process to manually download all of the 1000 Genomes WGS datasets based on hg38 to compile a new dataset on this genome build. So if your target dataset is not on the hg19 genome build, the simplest solution is to lift over the coordinates to hg19. An easy-to-use tool to directly lift over variants in VCF format is Picard's LiftoverVcf (https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-). Picard is already included in this repository's Anaconda environment, if this step is necessary.

## Merging reference and target datasets
Once the reference dataset has been compiled, the next step is to merge the reference individuals representing the global ancestries you want to estimate with the individuals from your target dataset. This workflow is detailed in downloaded from https://cran.r-project.org/web/packages/plinkQC/vignettes/AncestryCheck.pdf, but the steps described below should function in most cases.

While the 1000 Genomes dataset is not directly labeled with global population for each individual, the International Genome Sample Resource provides these data in a separate spreadsheet (https://www.internationalgenome.org/faq/can-i-get-phenotype-gender-and-family-relationship-information-samples/). A copy of these annotations in TSV format (```Sample_Info-Table_1.tsv```) is included in this repository.

The family IDs in the 1000 Genomes fam file we have compiled (in column 1) are coded as zeros, so you will next want to recode these with the corresponding global population for each individual instead. Modify and then run the ```annotate_fam_with_pops_1kg_b37.R``` script to perform this recoding:
```
#Declare the path of the 1KG fam file, the reference TSV file with the corresponding populations of individuals in the cohort, and the target output directory
fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/b37/all_phase3.fam"
pop_key_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/Sample_Info-Table_1.tsv"
out_dir <- "/home/bettimj/gamazon_rotation/1000_genomes/b37/pop_anno"

#Open each of the files as a data frame
fam_file <- read.table(fam_path, header = FALSE)
fam_df <- as.data.frame(fam_file)

pop_key_file <- read.delim(pop_key_path, header = TRUE)
pop_key_df <- as.data.frame(pop_key_file)

#For each match between sample IDs in the fam and key files, replace the family ID in the fam file with the 1000 Genomes population identifier from the key file
fam_ind_ids <- fam_df[,2]
pop_key_ind_ids <- pop_key_df[,1]

fam_df[,1] <- pop_key_df$Population[match(fam_ind_ids, pop_key_ind_ids)]

#Write this new fam file out to a fam file
out_file_name <- "all_phase3.fam"
write.table(fam_df, file = paste(out_dir, out_file_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The new fam file should be output to a different directory from your original 1000 Genomes dataset:
```
Rscript annotate_fam_with_pops_1kg_b37.R
``` 

Once the new fam file has been generated, create symbolic links from the original bed and bim files to the same directory as your new fam, such as in the example below:
```
ln -s /home/bettimj/gamazon_rotation/1000_genomes/b37/all_phase3.bed /home/bettimj/gamazon_rotation/1000_genomes/b37/pop_anno

ln -s /home/bettimj/gamazon_rotation/1000_genomes/b37/all_phase3.bim /home/bettimj/gamazon_rotation/1000_genomes/b37/pop_anno
```
At this point in the workflow, variant merging and pruning will fail if your reference and target datasets do not have the same naming scheme for their variable names (column 2 of the bim files). For example, the 1000 Genomes variants are annotated with rsIDs, but GTEx variants are annotated with varIDs. So if there is discordance between variant naming, the easiest solution is to create completely new varIDs for the bim files of both the reference and target dataset. An example of how this can be accomplished is in the ```gtexv8_new_varids_for_gtex_1kg_bims.R``` script:
```
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
```

Wherever these new bims are saved to, symbolic links from the corresponding bed and fam files should also be linked to the same directory:
```
ln -s /home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.bed \
/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids

ln -s /home/bettimj/gamazon_rotation/gtex8_geno/liftover_hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.fam \
/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids

ln -s /home/bettimj/gamazon_rotation/1000_genomes/b37/all_phase3.bed \
/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids

ln -s /home/bettimj/gamazon_rotation/1000_genomes/b37/pop_anno/all_phase3.fam \
/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids
```

Next, depending on how many reference populations you would like to include in your ADMIXTURE estimate, subset the 1000 Genomes plink files so that they retain only individuals of your global ancestries of interest. The ```mod_bim_plink_subset_1000genomes_ceu_yri.sh``` script details what pruning for only individuals of European (CEU) and African (YRI) ancestry would look like:
```
#The tutorial found on https://www.gnxp.com/WordPress/2018/07/13/tutorial-to-run-supervised-admixture-analyses/?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+RazibKhansTotalFeed+%28Razib+Khan%27s+total+feed%29 was helpful in generating this workflow.

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
SOURCE_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids
TARGET_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids/1kg_filt_ceu_yri
FAM_NAME=all_phase3

#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from either the ASW or CEU 
grep "CEU\|YRI" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/keep.keep

#Use plink to generate an entire new set of binary files containing just the CEU and YRI 1000 Genomes individuals
plink --bfile $SOURCE_DIR\/$FAM_NAME --keep $TARGET_DIR\/keep.keep --make-bed --allow-extra-chr --out $TARGET_DIR\/ceu_yri_phase3
```

The reference dataset can now be pruned and merged with the target data. The LD cutoff for pruning was an r2 threshold of 0.1. The ```gtexv8_merge_1kg_pop_ceu_yri.sh``` script details this step. The input files for the ```HIGHLD```variable are included in this repository and are borrowed from the plinkQC R package:
```
#Import plink 1.9, which is installed as a module on ACCRE
module load PLINK/1.9b_5.2

#Declare the paths of the processed AALC cohort and 1000 Genomes data within the Aldrich Lab folder on ACCRE
BED_IN_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids
BED_IN_NAME=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19
BED_1KG_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/bim_convert_rsids/1kg_filt_ceu_yri
BED_1KG_NAME=ceu_yri_phase3

#Declare the path of file containing genomic ranges of high-LD structure. This file is included with the plinkQC R package.
HIGHLD=/gpfs52/home/bettimj/miniconda3/lib/R/library/plinkQC/extdata/high-LD-regions-hg19-GRCh37.txt

#Create the necessary output directories for outputs
ROOT_DIR=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri

QCDIR=$ROOT_DIR\/qcdir

mkdir $QCDIR
mkdir $QCDIR\/plink_log

#Merge the cohort genotype data with the 1000 Genomes reference data using plink --merge
	#Filter reference and study data for non A-T ot G-C SNPs. This is because these SNPs are more difficult to align, and only a subset of SNPs is required for this type of analysis.
		#Input cohort
awk 'BEGIN {OFS="\t"}($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA"){print $2}' \
$BED_IN_DIR\/$BED_IN_NAME\.bim > \
$QCDIR\/$BED_IN_NAME\.ac_gt_snps

plink --bfile $BED_IN_DIR\/$BED_IN_NAME \
	--exclude $QCDIR/$BED_IN_NAME\.ac_gt_snps \
	--make-bed \
	--out $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
	--allow-extra-chr

mv $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.log $QCDIR\/plink_log/$BED_IN_NAME\.no_ac_gt_snps.log

		#1000 Genomes
awk 'BEGIN {OFS="\t"}($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA"){print $2}' \
$BED_1KG_DIR\/$BED_1KG_NAME\.bim > \
$QCDIR\/$BED_1KG_NAME\.ac_gt_snps

plink --bfile $BED_1KG_DIR\/$BED_1KG_NAME \
	--exclude $QCDIR\/$BED_1KG_NAME\.ac_gt_snps \
	--make-bed \
	--allow-extra-chr \
	--out $QCDIR\/$BED_1KG_NAME\.no_ac_gt_snps

mv $QCDIR\/$BED_1KG_NAME\.no_ac_gt_snps.log $QCDIR\/plink_log/$BED_1KG_NAME\.no_ac_gt_snps.log

	#Prune study data. Conduct PCA on genetic variants that are pruned for variants in LD with r2 > 0.1 in a 50 kb window.
		#Input cohort
plink --bfile $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
	--exclude range $HIGHLD \
	--indep-pairwise 50 5 0.1 \
	--out $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
	--allow-extra-chr
mv $QCDIR\/$BED_IN_NAME\.prune.log $QCDIR\/plink_log/$BED_IN_NAME\.prune.log

plink --bfile $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
	--extract $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.prune.in \
	--make-bed \
	--out $QCDIR\/$BED_IN_NAME\.pruned \
	--allow-extra-chr
mv $QCDIR\/$BED_IN_NAME\.pruned.log $QCDIR\/plink_log/$BED_IN_NAME\.pruned.log

	#Filter reference data for the same SNP set as in study. Use the list of pruned variants from the study sample to reduce the reference dataset to the size of the study samples
plink --bfile $BED_1KG_DIR\/$BED_1KG_NAME \
	--extract $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.prune.in \
	--make-bed \
	--allow-extra-chr \
	--out $QCDIR\/$BED_1KG_NAME\.pruned

mv $QCDIR\/$BED_1KG_NAME\.pruned.log $QCDIR\/plink_log/$BED_1KG_NAME\.pruned.log

	#Check and correct chromosome mismatch. Check that the variant IDs of the reference data have the same chromosome ID as the study data. Merging the files via plink will only work for variants with perfectly matching attributes. Because sex chromosomes and are often encoded differently and might make the matching more difficult, we will ignore sex chromosomes.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
	($2 in a && a[$2] != $1)	{print a[$2],$2}' \
	$QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim | \
	sed -n '/^[XY]/!p' > $QCDIR\/$BED_1KG_NAME\.toUpdateChr

plink --bfile $QCDIR\/$BED_1KG_NAME\.pruned \
	--update-chr $QCDIR\/$BED_1KG_NAME\.toUpdateChr 1 2 \
	--make-bed \
	--out $QCDIR\/$BED_1KG_NAME\.updateChr \
	--allow-extra-chr
mv $QCDIR\/$BED_1KG_NAME\.updateChr.log $QCDIR\/plink_log/$BED_1KG_NAME\.updateChr.log

	#Position mismatch - find variants with mis-matching chromosome positions
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
	($2 in a && a[$2] != $4)	{print a[$2],$2}' \
	$QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim > \
	$QCDIR\/$BED_1KG_NAME\.toUpdatePos

	#Possible allele flips - check if non-matching allele codes are a simple case of allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
	($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
	$QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim > \
    $QCDIR\/$BED_1KG_NAME\.toFlip
    
	#Update positions and flip alleles
plink --bfile $QCDIR\/$BED_1KG_NAME\.updateChr \
      --update-map $QCDIR\/$BED_1KG_NAME\.toUpdatePos 1 2 \
      --flip $QCDIR\/$BED_1KG_NAME\.toFlip \
      --make-bed \
      --out $QCDIR\/$BED_1KG_NAME\.flipped \
      --allow-extra-chr
      
mv $QCDIR\/$BED_1KG_NAME\.flipped.log $QCDIR\/plink_log/$BED_1KG_NAME\.flipped.log

	#Remove mismatches - any alleles that do not match after allele flipping are identified and removed from the reference dataset.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.flipped.bim > \
    $QCDIR\/$BED_1KG_NAME\.mismatch

plink --bfile $QCDIR\/$BED_1KG_NAME\.flipped \
      --exclude $QCDIR\/$BED_1KG_NAME\.mismatch \
      --make-bed \
      --out $QCDIR\/$BED_1KG_NAME\.clean \
      --allow-extra-chr
      
mv $QCDIR\/$BED_1KG_NAME\.clean.log $QCDIR\/plink_log/$BED_1KG_NAME\.clean.log

	#Merge study genotypes and reference data - merge the AALC study data and 1000 Genomes reference dataset into a combined dataset
plink --bfile $QCDIR\/$BED_IN_NAME\.pruned  \
      --bmerge $QCDIR\/$BED_1KG_NAME\.clean.bed $QCDIR\/$BED_1KG_NAME\.clean.bim \
         $QCDIR\/$BED_1KG_NAME\.clean.fam  \
      --make-bed \
      --out $ROOT_DIR\/$BED_IN_NAME\.merge.$BED_1KG_NAME \
      --allow-extra-chr
      
mv $ROOT_DIR\/$BED_IN_NAME\.merge.$BED_1KG_NAME\.log $QCDIR/plink_log
```

It is possible that the resulting set of merged plink files might contain some variants on non-autosomal chromosomes, such as chrX. Because ADMIXTURE only works with the canonical autosomal chromosomes (chr1-chr22), it will throw an error if any non-autosomal chromosomes are present. Plink can be used to filter these out, retaining only variants in the 22 autosomes:
```
plink --bfile \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3 \
--chr 1-22 \
--make-bed \
--out GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal \
--allow-extra-chr
```

## Running ADMIXTURE
ADMIXTURE requires the following two sets of files to run:
* A merged set of plink files (bed/bim/fam)
* A single-column pop file containing the global population of each individual, with the order corresponding to the order of individuals in the merged fam file. Target individuals should be denoted as "-".

All of these files should be together in the same directory.

A pop file can ge generated using the ```make_pop_file_gtexv8_1000genomes_ceu_yri.R``` script included in this repository:
```
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
```

After generating the pop file, the fam file should finally be modified one last time to replace the family IDs (FIDs) with the individual IDs (IIDs), so that the first and second column have identical values. This can be performed using the code in the ```gtexv8_replace_merged_fam_fid_col_with_iids.R``` script:
```
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
```

Finally, a supervised run of ADMIXTURE can be run using your merged set of plink binary files. The K value represents the number of reference populations you are using, so in the case of working with CEU and YRI, you would use a K value of 2. The ```gtexv8_admixture_supervised_1000genomes_ceu_yri.sh``` script provides an example of this ADMIXTURE run:
```
#Declare the path of the input bed. Also within this directory should be the accompanying bim, fam, and pop files
in_bed=/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.bed

#Run ADMIXTURE using the --supervised flag and a K of 2, as we are working with two known reference populations from 1000 Genomes in this analysis
admixture --supervised $in_bed 2
```

Once ADMIXTURE has successfully completed a run, its two main outputs are files ending in .P and .Q. The Q file is the one containing global ancestry estimates for each individual based on the two reference populations. This file contains only unlabeled estimates, and all of the reference individuals are also still present with annotations of 100% ancestry for their reference population.

So a script like ```gtexv8_q_output_ceu_yri_prune.R``` can be used to add the individual identifiers back to the respective estimates and prune out the reference individuals so that only those from your target population remain:
```
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
```

After generating this pruned, annotated file of global ancestry estimates, you might also want to sort based on ancestry proportion. The ```order_gtexv8_q_ceu_yri_pruned.R``` script contains code that orders results from highese proportion of African ancestry to lowest:
```
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
```

Your final results should look something like the ```GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.anno.pruned.ordered.Q``` file.

If you would like to plot these global ancestry calculations visually, the code in the ```gtexv8_ceu_yri_barplot_q_admixture.R``` script can be used as a template:
```
#Declare the directory of the pruned .Q file containing the Q ancestry estimates of the GTEx v8 individuals
pruned_q_path <- "/home/bettimj/aldrich_rotation/admixture/gtexv8_supervised_ceu_yri/plink_merge_gtexv8_1000Genomes_ceu_yri/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hg19.merge.ceu_yri_phase3_autosomal.2.anno.pruned.ordered.Q"

#Read in and plot the Q file, saving it as a pdf
pruned_q_file <- read.table(pruned_q_path, header = FALSE)

pdf("gtexv8_ceu_yri_admixture_q.pdf")
barplot(t(as.matrix(pruned_q_file[,2:3])), main = "GTEx v8 African and European Ancestry Estimates", col = rainbow(2),
               xlab = "GTEx v8 Individuals", ylab = "Ancestry", border = NA, legend = c("CEU", "YRI"), args.legend=list(
      x=700,
      y=0,
      bty = "n"
    ))
dev.off()
```

An example of what this resulting image would look like is depicted in the ```gtexv8_ceu_yri_admixture_q.pdf``` PDF image.
