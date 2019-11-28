setwd("/mnt/2EA01BDBA01BA7FB/Working issues/BI_Science/Systems Genetics/Pleiotropy/Paper/Revision/Figures/expression")

library(ggplot2)
library(reshape2)


# Incrementing degree
#############gtex########################################
#########################################################
#########################################################

all_exp = read.table('../../new_eqtl_tpm/gtex.tsv', sep='\t', header = T)

strat1 = read.table('../../final_loci/1_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat2 = read.table('../../final_loci/2_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat3 = read.table('../../final_loci/3_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat4 = read.table('../../final_loci/4_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat5 = read.table('../../final_loci/5_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat6 = read.table('../../final_loci/6_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat7 = read.table('../../final_loci/7_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat8 = read.table('../../final_loci/8_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat9 = read.table('../../final_loci/9_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)
strat10 = read.table('../../final_loci/10_noHLA_index_genes.GTEx.medians.tsv', sep='\t', header=F)


summaries_med = data.frame(tpm = c(apply(all_exp[, 3:56], 1, median),
                                   apply(strat1, 1, median),
                                   apply(strat2, 1, median),
                                   apply(strat3, 1, median),
                                   apply(strat4, 1, median),
                                   apply(strat5, 1, median),
                                   apply(strat6, 1, median),
                                   apply(strat7, 1, median),
                                   apply(strat8, 1, median),
                                   apply(strat9, 1, median),
                                   apply(strat10, 1, median)),
                           pleio = c(rep(0, nrow(all_exp)),
                                     rep(1, nrow(strat1)),
                                     rep(2, nrow(strat2)),
                                     rep(3, nrow(strat3)),
                                     rep(4, nrow(strat4)),
                                     rep(5, nrow(strat5)),
                                     rep(6, nrow(strat6)),
                                     rep(7, nrow(strat7)),
                                     rep(8, nrow(strat8)),
                                     rep(9, nrow(strat9)),
                                     rep(10, nrow(strat10))))

ggplot(summaries_med, aes(x=as.factor(pleio), y=log2(tpm+1), 
                          fill=as.factor(pleio))) + 
  geom_violin(scale='width') +
  geom_boxplot(width=0.4, fill='white', outlier.shape=NA) +
  scale_y_continuous(limits=c(0, 10)) + theme_bw() + guides(fill=F)
