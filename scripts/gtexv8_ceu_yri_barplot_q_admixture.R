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