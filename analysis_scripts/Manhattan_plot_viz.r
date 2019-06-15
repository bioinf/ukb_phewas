library(qqman)
library(dplyr)
library(ggplot2)
manhattan_file = read.csv('herit_Clustered_mahattan_pval_jaccard.tsv', header = T, sep = "\t") 
#signif_names = read.csv('herit_valuable_SNP_jaccard.tsv', header = FALSE, sep = "\t")
#str(manhattan_file)

#manhattan_file[which(manhattan_file$V3 %in% signif_names$V1),] %>% write.table(file = "herit_valuable_SNP_jaccard_for_func", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

#manhattan(manhattan_file, chr="V1", bp="V2", snp="V3", p="V5", logp = FALSE)
colnames(manhattan_file) = c('CHR', "BP", "SNP", "Assos_list", "Assos_num", 'logp')

don <- manhattan_file %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(manhattan_file, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


#require(grid)
#dev.new()
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



ggplot(don, aes(x=BPcum, y=-log10(logp))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.9, size=2.8) +
  scale_color_manual(values = rep(c("#808B96", "#A4EAE8"), 22)) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 35)) +
  coord_fixed((1 * max(don$BPcum) / (max(-log10(don$logp))))/2.5)+
  xlab('Genomic coordinate') + ylab('Number of associations')+theme_void()+
  theme( legend.position="none")
  



c("#808B96", "#A4EAE8")
c("#1DAEC5","#AF77BD")
manhattan_file %>% select(SNP, Assos_list, Assos_num) %>% write.table(file = "herit_clust_diags_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


require(grid)
# Make the plots
dev.new()  # Reducing width and height of this device will produce smaller raster files
plot(rnorm(2e4), rnorm(2e4), pch="+", cex=0.6)
cap1 <- grid.cap()
plot(rnorm(2e4), rnorm(2e4), pch="+", cex=0.6, col="red")
cap2 <- grid.cap()
dev.off()

# Write them to a pdf
pdf("test.pdf", onefile=TRUE)
grid.raster(cap1)
plot.new()
grid.raster(cap2)
dev.off()

#c("blue","yellow", "green", "grey", "skyblue", "purple","black","orange", "red", "cyan", "pink"), 2