library(dplyr)
library(ggplot2)
require(grid)

manhattan_file = read.csv('herit_Clustered_mahattan_pval_jaccard.tsv', header = T, sep = "\t") 
colnames(manhattan_file) = c('CHR', "BP", "SNP", "Assos_list", "Assos_num")


don <- manhattan_file %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(manhattan_file, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

dev.new()
pdf("_pval_clust_manhattan_2.pdf")
ggplot(don, aes(x=BPcum, y=Assos_num)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0, size=6) +
  scale_color_manual(values = rep(c("#1DAEC5","#AF77BD"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 30)) +
  coord_fixed(1 * max(don$BPcum) / (max(Assos_num)+100))+
  xlab('Genomic coordinate') + ylab('Number of associations')+
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=8,angle = 45),
    axis.line = element_line(colour = "black")
  )
dev.off()

manhattan_file %>% select(SNP, Assos_list, Assos_num) %>% write.table(file = "herit_clust_diags_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

#Для манхэттена пвалов:

colnames(manhattan_file) = c('CHR', "BP", "SNP", "Assos_list", "Assos_num", 'logp')
dev.new()
pdf("_pval_clust_manhattan_2.pdf")
ggplot(don, aes(x=BPcum, y=-log10(logp))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0, size=6) +
  scale_color_manual(values = rep(c("#1DAEC5","#AF77BD"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 30)) +
  coord_fixed(1 * max(don$BPcum) / (max(-log10(don$logp))+100))+
  xlab('Genomic coordinate') + ylab('Number of associations')+
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=8,angle = 45),
    axis.line = element_line(colour = "black")
  )
dev.off()


