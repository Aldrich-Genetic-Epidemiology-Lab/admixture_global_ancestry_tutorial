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